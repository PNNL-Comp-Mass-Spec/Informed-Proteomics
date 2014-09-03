using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public class RafLcMsRun: ILcMsRun, IMassSpecDataReader
    {
        public const int FileFormatId = 150603;

        public static ILcMsRun GetLcMsRun(string specFilePath, MassSpecDataType dataType)
        {
            return GetLcMsRun(specFilePath, dataType, 0.0, 0.0);
        }

        public static ILcMsRun GetLcMsRun(string specFilePath, MassSpecDataType dataType,
            double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold)
        {
            var rafFilePath = Path.ChangeExtension(specFilePath, ".raf");

            if (!File.Exists(rafFilePath) || !CheckFileFormatVersion(rafFilePath))
            {
                LcMsRun run;
                if (dataType == MassSpecDataType.XCaliburRun)
                {
                    run = new LcMsRun(new XCaliburReader(specFilePath),
                        0, 0);
                }
                else run = null;
                if (run != null) run.WriteAsRaf(rafFilePath);
                else throw new Exception("Unsupported raw file format!");
            }

            return new RafLcMsRun(rafFilePath, precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold);
        }

        public RafLcMsRun(string specFileName, double precursorSignalToNoiseRatioThreshold = 0.0, double productSignalToNoiseRatioThreshold = 0.0)
        {
            _reader = new BinaryReader(new BufferedStream(File.Open(specFileName, FileMode.Open, FileAccess.Read, FileShare.Read)));
            if (ReadMetaInfo() == false)
            {
                throw new FormatException("Illegal raf file format!");
            }
            _precursorSignalToNoiseRatioThreshold = precursorSignalToNoiseRatioThreshold;
            _productSignalToNoiseRatioThreshold = productSignalToNoiseRatioThreshold;
        }

        public int MinLcScan { get; private set; }
        public int MaxLcScan { get; private set; }

        public Spectrum GetSpectrum(int scanNum)
        {
            long offset;
            if (!_scanNumToSpecOffset.TryGetValue(scanNum, out offset)) return null;
            var spec = ReadSpectrum(offset);
            if (spec.MsLevel == 1 && _precursorSignalToNoiseRatioThreshold > 0.0) spec.FilterNoise(_precursorSignalToNoiseRatioThreshold);
            else if (_productSignalToNoiseRatioThreshold > 0.0) spec.FilterNoise(_productSignalToNoiseRatioThreshold);
            return spec;
        }

        public double GetElutionTime(int scanNum)
        {
            double elutionTime;
            return _scanNumElutionTimeMap.TryGetValue(scanNum, out elutionTime) ? elutionTime : 0.0;
        }

        public int GetMsLevel(int scanNum)
        {
            int msLevel;
            return _scanNumToMsLevel.TryGetValue(scanNum, out msLevel) ? msLevel : 0;
        }

        public int GetPrevScanNum(int scanNum, int msLevel)
        {
            for (var curScanNum = scanNum - 1; curScanNum >= MinLcScan; curScanNum--)
            {
                if (GetMsLevel(curScanNum) == msLevel) return curScanNum;
            }
            return MinLcScan - 1;
        }

        public int GetNextScanNum(int scanNum, int msLevel)
        {
            for (var curScanNum = scanNum + 1; curScanNum <= MaxLcScan; curScanNum++)
            {
                if (GetMsLevel(curScanNum) == msLevel)
                {
                    return curScanNum;
                }
            }
            return MaxLcScan + 1;
        }

        public IList<int> GetScanNumbers(int msLevel)
        {
            var scanNumbers = new List<int>();
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                if (GetMsLevel(scanNum) == msLevel) scanNumbers.Add(scanNum);
            }
            return scanNumbers;
        }

        public Xic GetFullPrecursorIonExtractedIonChromatogram(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetFullPrecursorIonExtractedIonChromatogram(minMz, maxMz);
        }

        public Xic GetFullPrecursorIonExtractedIonChromatogram(double minMz, double maxMz)
        {
            var xic = GetPrecursorExtractedIonChromatogram(minMz, maxMz);

            var hasXicPoint = new bool[MaxLcScan - MinLcScan + 1];
            foreach (var xicPoint in xic) hasXicPoint[xicPoint.ScanNum - MinLcScan] = true;

            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                if (GetMsLevel(scanNum) == 1 && !hasXicPoint[scanNum - MinLcScan]) xic.Add(new XicPoint(scanNum, 0, 0));
            }
            xic.Sort();
            return xic;
        }

        public Xic GetFullProductExtractedIonChromatogram(double mz, Tolerance tolerance, double precursorIonMz)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetFullProductExtractedIonChromatogram(minMz, maxMz, precursorIonMz);
        }

        public Xic GetFullProductExtractedIonChromatogram(double minMz, double maxMz, double precursorMz)
        {
            var targetOffset = GetOffset(minMz, maxMz, _offsetProductChromatogramBebin, _offsetProductChromatogramEnd);
            if (targetOffset < _offsetProductChromatogramBebin) return new Xic();
            var xic = GetXicPointsWithin(minMz, maxMz, _offsetProductChromatogramBebin, _offsetProductChromatogramEnd, targetOffset);
            if (!xic.Any()) return xic;

            var scanToXicPoint = new XicPoint[MaxLcScan - MinLcScan + 1];
            foreach (var xicPoint in xic)
            {
                var prev = scanToXicPoint[xicPoint.ScanNum - MinLcScan];
                if (prev == null || xicPoint.Intensity > prev.Intensity) scanToXicPoint[xicPoint.ScanNum - MinLcScan] = xicPoint;
            }

            var newXic = new Xic();
            
            newXic.AddRange(GetFragmentationSpectraScanNums(precursorMz).Select(scanNum => scanToXicPoint[scanNum - MinLcScan] ?? new XicPoint(scanNum, 0, 0)));
            return newXic;
        }

        /// <summary>
        /// Gets scan numbers of the fragmentation spectra whose isolation window contains the precursor ion specified
        /// </summary>
        /// <param name="mostAbundantIsotopeMz"></param>
        /// <returns>scan numbers of fragmentation spectra</returns>
        public int[] GetFragmentationSpectraScanNums(double mostAbundantIsotopeMz)
        {
            var targetIsoBin = (int)Math.Round(mostAbundantIsotopeMz * LcMsRun.IsolationWindowBinningFactor);
            int[] scanNums;
            return _isolationMzBinToScanNums.TryGetValue(targetIsoBin, out scanNums) ? scanNums : new int[0];
        }

        private readonly BinaryReader _reader;
        private readonly double _precursorSignalToNoiseRatioThreshold;
        private readonly double _productSignalToNoiseRatioThreshold;

        private long _offsetChromatogramStart;
        private long _offsetProductChromatogramBebin;
        private long _offsetProductChromatogramEnd;
        private long _offsetChromatogramEnd;
        //private long _offsetMetaInfo;

        private Dictionary<int, int> _scanNumToMsLevel;
        private Dictionary<int, double> _scanNumElutionTimeMap;
        private Dictionary<int, int[]> _isolationMzBinToScanNums;
        private Dictionary<int, long> _scanNumToSpecOffset;
        private int _minMzIndex;
        private int _maxMzIndex;

        private long[] _chromMzIndexToOffset;
        private const double MzBinSize = 1;
        private const int NumBytePeak = 16;

        public bool ReadMetaInfo()
        {
            _reader.BaseStream.Seek(-1 * sizeof(int), SeekOrigin.End);
            var fileFormatId = _reader.ReadInt32();
            if (fileFormatId != FileFormatId) return false;

            _reader.BaseStream.Seek(-3*sizeof (long) - 1*sizeof (int), SeekOrigin.End);

            _offsetChromatogramStart = _reader.ReadInt64();
            _offsetChromatogramEnd = _offsetProductChromatogramBebin = _reader.ReadInt64();
            _offsetProductChromatogramEnd = _reader.ReadInt64();
            // Read meta information
            var offsetMetaInfo = _offsetProductChromatogramEnd;

            _reader.BaseStream.Seek(offsetMetaInfo, SeekOrigin.Begin);
            MinLcScan = _reader.ReadInt32();
            MaxLcScan = _reader.ReadInt32();

            _scanNumToMsLevel = new Dictionary<int, int>();
            _scanNumElutionTimeMap = new Dictionary<int, double>();
            _scanNumToSpecOffset = new Dictionary<int, long>();
            //_scanNumToIsolationWindow = new Dictionary<int, IsolationWindow>();
            //_isolationMzBinToScanNums = new Dictionary<int, List<int>>();
            var isolationMzBinToScanNums = new Dictionary<int, List<int>>();
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                var msLevel = _reader.ReadInt32();
                _scanNumToMsLevel[scanNum] = msLevel;
                _scanNumElutionTimeMap[scanNum] = _reader.ReadDouble();
                if (msLevel == 2)
                {
                    var minMz = _reader.ReadSingle();
                    var maxMz = _reader.ReadSingle();
                    var minBinNum = (int)Math.Round(minMz * LcMsRun.IsolationWindowBinningFactor);
                    var maxBinNum = (int)Math.Round(maxMz * LcMsRun.IsolationWindowBinningFactor);
                    for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
                    {
                        List<int> scanNumList;
                        if (!isolationMzBinToScanNums.TryGetValue(binNum, out scanNumList))
                        {
                            scanNumList = new List<int>();
                            isolationMzBinToScanNums[binNum] = scanNumList;
                        }
                        scanNumList.Add(scanNum);
                    }
                }
                _scanNumToSpecOffset[scanNum] = _reader.ReadInt64();
            }

            _isolationMzBinToScanNums = new Dictionary<int, int[]>();
            foreach (var entry in isolationMzBinToScanNums)
            {
                var binNum = entry.Key;
                entry.Value.Sort();
                var scanNumList = entry.Value.ToArray();
                _isolationMzBinToScanNums[binNum] = scanNumList;
            }

            var minIndex = _reader.ReadInt32();
            var maxIndex = _reader.ReadInt32();
            _chromMzIndexToOffset = new long[maxIndex - minIndex + 1];

            for (var i = 0; i < _chromMzIndexToOffset.Length; i++)
            {
                _chromMzIndexToOffset[i] = _reader.ReadInt64();
            }
            _minMzIndex = minIndex;
            _maxMzIndex = maxIndex;

            return true;
        }

        public static int GetMzBinIndex(double mz)
        {
            return (int) (mz / MzBinSize);
        }

        public Xic GetPrecursorExtractedIonChromatogram(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetPrecursorExtractedIonChromatogram(minMz, maxMz);
        }

        public Xic GetPrecursorExtractedIonChromatogram(double minMz, double maxMz)
        {
            var minBinIndex = GetMzBinIndex(minMz);
            var maxBinIndex = GetMzBinIndex(maxMz);

            long targetOffset;
            if (minBinIndex == maxBinIndex)
            {
                if (maxBinIndex < _minMzIndex || maxBinIndex > _maxMzIndex) return new Xic();
                var offset = _chromMzIndexToOffset[maxBinIndex - _minMzIndex];
                if(offset < _offsetChromatogramStart) return new Xic();

                // binary search
                var beginOffset = offset;
                var endOffset = _chromMzIndexToOffset[maxBinIndex - _minMzIndex + 1];
                targetOffset = GetOffset(minMz, maxMz, beginOffset, endOffset);
            }
            else
            {
                if(maxBinIndex < _minMzIndex || minBinIndex > _maxMzIndex) return new Xic();
                targetOffset = maxBinIndex > _maxMzIndex ? _offsetChromatogramEnd : _chromMzIndexToOffset[maxBinIndex - _minMzIndex];
            }

            if(targetOffset < _offsetChromatogramStart) return new Xic();
            var xic = GetXic(minMz, maxMz, _offsetChromatogramStart, _offsetChromatogramEnd, targetOffset);
            if (!xic.Any()) return xic;
            xic.Sort();
            return Xic.GetSelectedXic(xic);
        }

        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            return _scanNumToSpecOffset.OrderBy(e => e.Key).Select(e => ReadSpectrum(e.Value));
        }

        public static bool CheckFileFormatVersion(string filePath)
        {
            var fs = File.OpenRead(filePath);
            using (var reader = new BinaryReader(fs))
            {
                fs.Seek(-1 * sizeof(int), SeekOrigin.End);

                var fileFormatId = reader.ReadInt32();
                if (fileFormatId != FileFormatId) return false;
            }
            return true;
        }

        public void Close()
        {
            _reader.Close();
        }

        private Spectrum ReadSpectrum(long offset)
        {
            _reader.BaseStream.Seek(offset, SeekOrigin.Begin);
            while (_reader.BaseStream.Position != (_reader.BaseStream.Length - sizeof(int)))
            {
                var spec = ReadSpectrum(_reader);
                return spec;
            }
            return null;
        }

        private static Spectrum ReadSpectrum(BinaryReader reader)
        {
            var scanNum = reader.ReadInt32();
            var msLevel = reader.ReadByte();
            var elutionTime = reader.ReadDouble();

            if (msLevel > 1)
            {
                double? precursorMass = reader.ReadDouble();
                if (precursorMass == 0.0) precursorMass = null;
                int? precursorCharge = reader.ReadInt32();
                if (precursorCharge == 0) precursorCharge = null;
                var activationMethod = (ActivationMethod)reader.ReadByte();
                var isolationWindowTargetMz = reader.ReadDouble();
                var isolationWindowLowerOffset = reader.ReadDouble();
                var isolationWindowUpperOffset = reader.ReadDouble();
                var peakList = ReadPeakList(reader);
                return new ProductSpectrum(peakList, scanNum)
                {
                    MsLevel = msLevel,
                    ElutionTime = elutionTime,
                    ActivationMethod = activationMethod,
                    IsolationWindow = new IsolationWindow(
                        isolationWindowTargetMz,
                        isolationWindowLowerOffset,
                        isolationWindowUpperOffset,
                        precursorMass,
                        precursorCharge
                        )
                };
            }
            else
            {
                var peakList = ReadPeakList(reader);
                return new Spectrum(peakList, scanNum)
                {
                    ElutionTime = elutionTime
                };
            }
        }

        private static List<Peak> ReadPeakList(BinaryReader reader)
        {
            var peakList = new List<Peak>();
            var numPeaks = reader.ReadInt32();
            for (var i = 0; i < numPeaks; i++)
            {
                var mz = reader.ReadDouble();
                var intensity = reader.ReadSingle();
                peakList.Add(new Peak(mz, intensity));
            }
            return peakList;
        }

        private long GetOffset(double minMz, double maxMz, long beginOffset, long endOffset)
        {
            var minOffset = beginOffset;
            var maxOffset = endOffset;
            var curOffset = -1L;
            while (minOffset <= maxOffset)
            {
                curOffset = minOffset + (maxOffset - minOffset) / NumBytePeak / 2 * NumBytePeak;
                _reader.BaseStream.Seek(curOffset, SeekOrigin.Begin);
                var curMz = _reader.ReadDouble();
                if (curMz < minMz) minOffset = curOffset + NumBytePeak;
                else if (curMz > maxMz) maxOffset = curOffset - NumBytePeak;
                else
                {
                    return curOffset;
                }
            }
            return curOffset;
        }

        // beginOffset: inclusive, endOffset: exclusive
        private Xic GetXic(double minMz, double maxMz, long beginOffset, long endOffset, long targetOffset)
        {
            var xic = GetXicPointsWithin(minMz, maxMz, beginOffset, endOffset, targetOffset);
            if (!xic.Any()) return xic;

            xic.Sort();
            return Xic.GetSelectedXic(xic);
        }

        private Xic GetXicPointsWithin(double minMz, double maxMz, long beginOffset, long endOffset,
            long targetOffset)
        {
            var xic = new Xic();
            var curOffset = targetOffset - NumBytePeak;
            // go down
            while (curOffset >= beginOffset)
            {
                _reader.BaseStream.Seek(curOffset, SeekOrigin.Begin);
                var mz = _reader.ReadDouble();
                var intensity = _reader.ReadSingle();
                var scanNum = _reader.ReadInt32();
                if (mz < minMz) break;
                xic.Add(new XicPoint(scanNum, mz, intensity));
                curOffset -= NumBytePeak;
            }

            // go up
            curOffset = targetOffset;
            while (curOffset < endOffset)
            {
                _reader.BaseStream.Seek(curOffset, SeekOrigin.Begin);
                var mz = _reader.ReadDouble();
                var intensity = _reader.ReadSingle();
                var scanNum = _reader.ReadInt32();
                if (mz > maxMz) break;
                xic.Add(new XicPoint(scanNum, mz, intensity));
                _reader.BaseStream.Seek(NumBytePeak, SeekOrigin.Current);
                curOffset += NumBytePeak;
            }

            return xic;
        }
    }
}
