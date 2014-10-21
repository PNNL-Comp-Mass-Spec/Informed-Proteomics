using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public class PbfLcMsRun: LcMsRun, IMassSpecDataReader
    {
        public const int FileFormatId = 150604;
        public const string FileExtension = ".pbf";
        
        public static LcMsRun GetLcMsRun(string specFilePath, MassSpecDataType dataType = MassSpecDataType.XCaliburRun)
        {
            return GetLcMsRun(specFilePath, dataType, 0.0, 0.0);
        }

        public static LcMsRun GetLcMsRun(string specFilePath, MassSpecDataType dataType,
            double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold)
        {
            var pbfFilePath = Path.ChangeExtension(specFilePath, FileExtension);

            if (!File.Exists(pbfFilePath) || !CheckFileFormatVersion(pbfFilePath))
            {
                InMemoryLcMsRun run;
                if (dataType == MassSpecDataType.XCaliburRun)
                {
                    run = new InMemoryLcMsRun(new XCaliburReader(specFilePath),
                        0, 0);
                }
	            else if (dataType == MassSpecDataType.MzMLFile)
	            {
		            run = new InMemoryLcMsRun(new MzMLReader(specFilePath), 0, 0);
	            }
                else run = null;
                if (run == null) throw new Exception("Unsupported raw file format!");
                try
                {
                    run.WriteAsPbf(pbfFilePath);
                }
                catch (UnauthorizedAccessException) // Cannot write to same directory, attemp to write to temp directory
                {
                    var fileName = Path.GetFileName(pbfFilePath);
                    if (String.IsNullOrEmpty(fileName)) throw;  // invalid path?
                    var tempPath = Path.Combine(Path.GetTempPath(), fileName);
                    if (!File.Exists(tempPath) || !CheckFileFormatVersion(tempPath)) run.WriteAsPbf(tempPath);
                    pbfFilePath = tempPath;
                }
            }

            return new PbfLcMsRun(pbfFilePath, precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold);
        }

        public PbfLcMsRun(string specFileName, double precursorSignalToNoiseRatioThreshold = 0.0, double productSignalToNoiseRatioThreshold = 0.0)
        {
            _reader = new BinaryReader(new BufferedStream(File.Open(specFileName, FileMode.Open, FileAccess.Read, FileShare.Read)));
            if (ReadMetaInfo() == false)
            {
                throw new FormatException("Illegal pbf file format!");
            }
            CreatePrecursorNextScanMap();
            _precursorSignalToNoiseRatioThreshold = precursorSignalToNoiseRatioThreshold;
            _productSignalToNoiseRatioThreshold = productSignalToNoiseRatioThreshold;


            _reader.BaseStream.Seek(_offsetPrecursorChromatogramStart, SeekOrigin.Begin);
            _minMs1Mz = _reader.ReadDouble();

            _reader.BaseStream.Seek(_offsetPrecursorChromatogramEnd - NumBytePeak, SeekOrigin.Begin);
            _maxMs1Mz = _reader.ReadDouble();
        }

        public override double MinMs1Mz
        {
            get { return _minMs1Mz; }
        }

        public override double MaxMs1Mz
        {
            get { return _maxMs1Mz; }
        }

        public override Spectrum GetSpectrum(int scanNum)
        {
            long offset;
            if (!_scanNumToSpecOffset.TryGetValue(scanNum, out offset)) return null;
            var spec = ReadSpectrum(offset);
            if (spec.MsLevel == 1 && _precursorSignalToNoiseRatioThreshold > 0.0) spec.FilterNoise(_precursorSignalToNoiseRatioThreshold);
            else if (_productSignalToNoiseRatioThreshold > 0.0) spec.FilterNoise(_productSignalToNoiseRatioThreshold);
            return spec;
        }


        public override Xic GetFullProductExtractedIonChromatogram(double minMz, double maxMz, double precursorMz)
        {
            var targetOffset = GetOffset(minMz, maxMz, _offsetProductChromatogramBegin, _offsetProductChromatogramEnd);
            if (targetOffset < _offsetProductChromatogramBegin) return new Xic();
            var xic = GetXicPointsWithin(minMz, maxMz, _offsetProductChromatogramBegin, _offsetProductChromatogramEnd, targetOffset);
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

        public override Xic GetPrecursorExtractedIonChromatogram(double minMz, double maxMz)
        {
            var minBinIndex = GetMzBinIndex(minMz);
            var maxBinIndex = GetMzBinIndex(maxMz);

            long targetOffset;
            if (minBinIndex == maxBinIndex)
            {
                if (maxBinIndex < _minMzIndex || maxBinIndex > _maxMzIndex) return new Xic();
                var offset = _chromMzIndexToOffset[maxBinIndex - _minMzIndex];
                if (offset < _offsetPrecursorChromatogramStart) return new Xic();

                // binary search
                var beginOffset = offset;
                var endOffset = _chromMzIndexToOffset[maxBinIndex - _minMzIndex + 1];
                targetOffset = GetOffset(minMz, maxMz, beginOffset, endOffset);
            }
            else
            {
                if (maxBinIndex < _minMzIndex || minBinIndex > _maxMzIndex) return new Xic();
                targetOffset = maxBinIndex > _maxMzIndex ? _offsetPrecursorChromatogramEnd : _chromMzIndexToOffset[maxBinIndex - _minMzIndex];
            }

            if (targetOffset < _offsetPrecursorChromatogramStart) return new Xic();
            var xic = GetXic(minMz, maxMz, _offsetPrecursorChromatogramStart, _offsetPrecursorChromatogramEnd, targetOffset);
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



        private readonly BinaryReader _reader;
        private readonly double _precursorSignalToNoiseRatioThreshold;
        private readonly double _productSignalToNoiseRatioThreshold;

        private readonly double _minMs1Mz;
        private readonly double _maxMs1Mz;

        private readonly Object _lock = new Object();
        private long _offsetPrecursorChromatogramStart;
        private long _offsetPrecursorChromatogramEnd;   // exclusive
        private long _offsetProductChromatogramBegin;
        private long _offsetProductChromatogramEnd;
        //private long _offsetMetaInfo;

//        private Dictionary<int, int> _scanNumToMsLevel;
//        private Dictionary<int, double> _scanNumElutionTimeMap;
//        private Dictionary<int, int[]> _isolationMzBinToScanNums;

        private Dictionary<int, long> _scanNumToSpecOffset;
        private int _minMzIndex;
        private int _maxMzIndex;

        private long[] _chromMzIndexToOffset;
        private const double MzBinSize = 1;
        private const int NumBytePeak = 16;

        private bool ReadMetaInfo()
        {
            ScanNumToMsLevel = new Dictionary<int, int>();
            ScanNumElutionTimeMap = new Dictionary<int, double>();
            IsolationMzBinToScanNums = new Dictionary<int, int[]>();

            _reader.BaseStream.Seek(-1 * sizeof(int), SeekOrigin.End);
            var fileFormatId = _reader.ReadInt32();
            if (fileFormatId != FileFormatId) return false;

            _reader.BaseStream.Seek(-3*sizeof (long) - 1*sizeof (int), SeekOrigin.End);

            _offsetPrecursorChromatogramStart = _reader.ReadInt64();
            _offsetPrecursorChromatogramEnd = _offsetProductChromatogramBegin = _reader.ReadInt64();
            _offsetProductChromatogramEnd = _reader.ReadInt64();
            // Read meta information
            var offsetMetaInfo = _offsetProductChromatogramEnd;

            _reader.BaseStream.Seek(offsetMetaInfo, SeekOrigin.Begin);
            MinLcScan = _reader.ReadInt32();
            MaxLcScan = _reader.ReadInt32();

            _scanNumToSpecOffset = new Dictionary<int, long>();
            //_scanNumToIsolationWindow = new Dictionary<int, IsolationWindow>();
            //_isolationMzBinToScanNums = new Dictionary<int, List<int>>();
            var isoWindowSet = new HashSet<IsolationWindow>();
            var isDda = false;
            var isolationMzBinToScanNums = new Dictionary<int, List<int>>();
            var minMsLevel = int.MaxValue;
            var maxMsLevel = int.MinValue;
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                var msLevel = _reader.ReadInt32();
                if (msLevel < minMsLevel) minMsLevel = msLevel;
                if (msLevel > maxMsLevel) maxMsLevel = msLevel;

                ScanNumToMsLevel[scanNum] = msLevel;
                ScanNumElutionTimeMap[scanNum] = _reader.ReadDouble();
                if (msLevel == 2)
                {
                    var minMz = _reader.ReadSingle();
                    var maxMz = _reader.ReadSingle();
                    if (!isDda)
                    {
                        var isoWindow = new IsolationWindow((minMz + maxMz) / 2, (maxMz - minMz) / 2, (maxMz - minMz) / 2);
                        isoWindowSet.Add(isoWindow);
                        if (isoWindowSet.Count >= NumUniqueIsolationWindowThresholdForDia) isDda = true;
                    }
                    var minBinNum = (int)Math.Round(minMz * IsolationWindowBinningFactor);
                    var maxBinNum = (int)Math.Round(maxMz * IsolationWindowBinningFactor);
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

            MinMsLevel = minMsLevel;
            MaxMsLevel = maxMsLevel;

            IsDiaOrNull = !isDda;

            foreach (var entry in isolationMzBinToScanNums)
            {
                var binNum = entry.Key;
                entry.Value.Sort();
                var scanNumList = entry.Value.ToArray();
                IsolationMzBinToScanNums[binNum] = scanNumList;
            }

            var minIndex = _reader.ReadInt32();
            var maxIndex = _reader.ReadInt32();
            _chromMzIndexToOffset = new long[maxIndex - minIndex + 2];

            for (var i = 0; i < _chromMzIndexToOffset.Length; i++)
            {
                _chromMzIndexToOffset[i] = _reader.ReadInt64();
            }
            _chromMzIndexToOffset[_chromMzIndexToOffset.Length - 1] = _offsetPrecursorChromatogramEnd;
            _minMzIndex = minIndex;
            _maxMzIndex = maxIndex;

            return true;
        }

        public static int GetMzBinIndex(double mz)
        {
            return (int) (mz / MzBinSize);
        }

        private Spectrum ReadSpectrum(long offset)
        {
            lock (_lock)
            {
                _reader.BaseStream.Seek(offset, SeekOrigin.Begin);
                while (_reader.BaseStream.Position != (_reader.BaseStream.Length - sizeof(int)))
                {
                    var spec = ReadSpectrum(_reader);
                    return spec;
                }
                return null;
            }
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
            // binary search
            lock (_lock)
            {
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
            lock (_lock)
            {
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
            }

            return xic;
        }
    }
}
