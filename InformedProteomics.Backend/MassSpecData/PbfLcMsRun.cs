using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.MassSpecData
{
    public class PbfLcMsRun: LcMsRun, IMassSpecDataReader
    {
        // This constant should be incremented by 1 if the binary file format is changed
        public const int FileFormatId = 150604;

        public const string FileExtension = ".pbf";

        /// <summary>
        /// Function to convert a spectra file name/path to a *.pbf name, even when it has multiple extensions (i.e., .mzML.gz)
        /// </summary>
        /// <param name="specFileName"></param>
        /// <returns></returns>
        /// <remarks>It is recommended that "MassSpecDataReaderFactory.NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path</remarks>
        public static string GetPbfFileName(string specFileName)
        {
            return MassSpecDataReaderFactory.ChangeExtension(specFileName, FileExtension);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="specFilePath"></param>
        /// <param name="dataType"></param>
        /// <param name="progress"></param>
        /// <returns></returns>
        /// <remarks>It is recommended that "MassSpecDataReaderFactory.NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path</remarks>
        [ObsoleteAttribute("Remove MassSpecDataType -> now uses MassSpecDataReaderFactory", true)]
        public static LcMsRun GetLcMsRun(string specFilePath, MassSpecDataType dataType, IProgress<ProgressData> progress = null)
        {
            return GetLcMsRun(specFilePath, dataType, 0.0, 0.0, progress);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="specFilePath"></param>
        /// <param name="progress"></param>
        /// <returns></returns>
        /// <remarks>It is recommended that "MassSpecDataReaderFactory.NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path</remarks>
        public static LcMsRun GetLcMsRun(string specFilePath, IProgress<ProgressData> progress = null)
        {
            return GetLcMsRun(specFilePath, 0.0, 0.0, progress);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="specFilePath"></param>
        /// <param name="dataType"></param>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        /// <param name="progress"></param>
        /// <returns></returns>
        /// <remarks>It is recommended that "MassSpecDataReaderFactory.NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path</remarks>
        [ObsoleteAttribute("Remove MassSpecDataType -> now uses MassSpecDataReaderFactory", true)]
        public static LcMsRun GetLcMsRun(string specFilePath, MassSpecDataType dataType,
            double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold, IProgress<ProgressData> progress = null)
        {
            //var pbfFilePath = InMemoryLcMsRun.ConvertToPbf(specFilePath, dataType, precursorSignalToNoiseRatioThreshold,
            //    productSignalToNoiseRatioThreshold, null, progress);
            //
            //return new PbfLcMsRun(pbfFilePath, precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold);
            return GetLcMsRun(specFilePath, MassSpecDataReaderFactory.GetMassSpecDataReader(specFilePath), precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold, progress);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="specFilePath"></param>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        /// <param name="progress"></param>
        /// <returns></returns>
        /// <remarks>It is recommended that "MassSpecDataReaderFactory.NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path</remarks>
        public static LcMsRun GetLcMsRun(
            string specFilePath, 
            double precursorSignalToNoiseRatioThreshold, 
            double productSignalToNoiseRatioThreshold, 
            IProgress<ProgressData> progress = null)
        {
            return GetLcMsRun(specFilePath, MassSpecDataReaderFactory.GetMassSpecDataReader(specFilePath), precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold, progress);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="specFilePath"></param>
        /// <param name="specReader"></param>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        /// <param name="progress"></param>
        /// <returns></returns>
        /// <remarks>It is recommended that "MassSpecDataReaderFactory.NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path</remarks>
        public static LcMsRun GetLcMsRun(string specFilePath, IMassSpecDataReader specReader, double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold,
            IProgress<ProgressData> progress = null)
        {
            var pbfFilePath = ConvertToPbf(specFilePath, specReader, precursorSignalToNoiseRatioThreshold,
                productSignalToNoiseRatioThreshold, null, progress);

            return new PbfLcMsRun(pbfFilePath, precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="specFilePath"></param>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        /// <param name="pbfFilePath"></param>
        /// <param name="progress"></param>
        /// <returns></returns>
        /// <remarks>It is recommended that "MassSpecDataReaderFactory.NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path</remarks>
        public static string ConvertToPbf(string specFilePath, double precursorSignalToNoiseRatioThreshold,
            double productSignalToNoiseRatioThreshold, string pbfFilePath = null, IProgress<ProgressData> progress = null)
        {
            return ConvertToPbf(specFilePath, MassSpecDataReaderFactory.GetMassSpecDataReader(specFilePath),
                precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold, pbfFilePath, progress);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="specFilePath"></param>
        /// <param name="specReader"></param>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        /// <param name="pbfFilePath">If supplied, file will be written to this path; otherwise the file will be written to the same directory as specFilePath, or to the temp directory if the user does not have write permissions</param>
        /// <param name="progress">Progress data, as a percentage</param>
        /// <returns></returns>
        /// <remarks>It is recommended that "MassSpecDataReaderFactory.NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path</remarks>
        public static string ConvertToPbf(string specFilePath, IMassSpecDataReader specReader,
            double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold, string pbfFilePath = null,
            IProgress<ProgressData> progress = null)
        {
            if (specFilePath.ToLower().EndsWith(PbfLcMsRun.FileExtension))
            {
                return specFilePath;
            }

            string pbfPath = pbfFilePath;
            string fileName = String.Empty;
            string tempPath = String.Empty;

            Progress<ProgressData> prog = new Progress<ProgressData>();
            var progData = new ProgressData();
            progData.IsPartialRange = true;
            progData.MaxPercentage = 75.0;
            if (progress != null)
            {
                prog = new Progress<ProgressData>(p =>
                {
                    progData.Status = p.Status;
                    progress.Report(progData.UpdatePercent(p.Percent));
                });
            }

            if (string.IsNullOrWhiteSpace(pbfFilePath))
            {
                // Calls "NormalizeDatasetPath" to make sure we save the file to the containing directory
                pbfPath = PbfLcMsRun.GetPbfFileName(MassSpecDataReaderFactory.NormalizeDatasetPath(specFilePath));
                fileName = Path.GetFileName(pbfPath);
                if (String.IsNullOrEmpty(fileName))
                {
                    throw new ArgumentException("Cannot create .pbf cache file", "specFilePath");
                }

                tempPath = Path.Combine(Path.GetTempPath(), fileName);
                // Return the temp path if the pbf file of proper format already exists in the temp directory
                if (File.Exists(tempPath) && PbfLcMsRun.CheckFileFormatVersion(tempPath))
                {
                    return tempPath;
                }
            }

            if (!File.Exists(pbfPath) || !PbfLcMsRun.CheckFileFormatVersion(pbfPath))
            {
                if (specReader == null)
                {
                    throw new Exception("Unsupported file format!");
                }
                InMemoryLcMsRun run = new InMemoryLcMsRun(specReader, 0, 0, prog);
                try
                {
                    progData.StepRange(100.0);
                    PbfLcMsRun.WriteAsPbf(run, pbfPath, prog);
                }
                catch (UnauthorizedAccessException) // Cannot write to same directory, attempt to write to temp directory
                {
                    if (string.IsNullOrWhiteSpace(pbfFilePath))
                    {
                        throw;
                    }
                    //var fileName = Path.GetFileName(pbfFilePath);
                    if (String.IsNullOrEmpty(fileName)) throw; // invalid path?
                    //var tempPath = Path.Combine(Path.GetTempPath(), fileName);
                    if (!File.Exists(tempPath) || !PbfLcMsRun.CheckFileFormatVersion(tempPath))
                        PbfLcMsRun.WriteAsPbf(run, tempPath, prog);
                    pbfPath = tempPath;
                }
            }
            return pbfPath;
        }

        public PbfLcMsRun(string specFileName, double precursorSignalToNoiseRatioThreshold = 0.0, double productSignalToNoiseRatioThreshold = 0.0)
        {
            var specFile = new FileInfo(specFileName);
            if (!specFile.Exists)
            {
                throw new FileNotFoundException("File not found by PbfLcMsRun", specFile.FullName);
            }

            if (specFile.Length < 28)
            {
                throw new FormatException("Illegal pbf file (too small)!");
            }

            _reader = new BinaryReader(new BufferedStream(File.Open(specFileName, FileMode.Open, FileAccess.Read, FileShare.Read)));
            _fileLock = new Object();

            _precursorSignalToNoiseRatioThreshold = precursorSignalToNoiseRatioThreshold;
            _productSignalToNoiseRatioThreshold = productSignalToNoiseRatioThreshold;

            lock(_fileLock)
            {
                if (ReadMetaInfo() == false)
                {
                    throw new FormatException("Illegal pbf file format!");
                }
                _reader.BaseStream.Seek(_offsetPrecursorChromatogramStart, SeekOrigin.Begin);
                _minMs1Mz = _reader.ReadDouble();

                if (_offsetPrecursorChromatogramEnd - NumBytePeak < 0)
                {
                    throw new FormatException("Corrupt pbf file (_offsetPrecursorChromatogramEnd is < 0)");
                }

                if (_offsetPrecursorChromatogramEnd - NumBytePeak >= specFile.Length)
                {
                    throw new FormatException("Corrupt pbf file (_offsetPrecursorChromatogramEnd is past the end of the file)");
                }

                _reader.BaseStream.Seek(_offsetPrecursorChromatogramEnd - NumBytePeak, SeekOrigin.Begin);
                _maxMs1Mz = _reader.ReadDouble();
            }

            NumSpectra = _scanNumToSpecOffset.Count;

            CreatePrecursorNextScanMap();
        }

        public override double MinMs1Mz
        {
            get { return _minMs1Mz; }
        }

        public override double MaxMs1Mz
        {
            get { return _maxMs1Mz; }
        }

        public override Spectrum GetSpectrum(int scanNum, bool includePeaks = true)
        {
            long offset;
            if (!_scanNumToSpecOffset.TryGetValue(scanNum, out offset)) return null;
            var spec = ReadSpectrum(offset, includePeaks);
            if (spec.MsLevel == 1 && _precursorSignalToNoiseRatioThreshold > 0.0) spec.FilterNoise(_precursorSignalToNoiseRatioThreshold);
            else if (_productSignalToNoiseRatioThreshold > 0.0) spec.FilterNoise(_productSignalToNoiseRatioThreshold);
            return spec;
        }

        public Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true)
        {
            return GetSpectrum(scanNum, includePeaks);
        }

        public bool TryMakeRandomAccessCapable()
        {
            return true;
        }

        public override IsolationWindow GetIsolationWindow(int scanNum)
        {
            long offset;
            return !_scanNumToSpecOffset.TryGetValue(scanNum, out offset) ? null : ReadIsolationWindow(offset);
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
            var pbfFile = new FileInfo(filePath);
            if (!pbfFile.Exists || pbfFile.Length < sizeof(int))
                return false;

            var fs = new FileStream(pbfFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            using (var reader = new BinaryReader(fs))
            {
                fs.Seek(-1 * sizeof(int), SeekOrigin.End);

                var fileFormatId = reader.ReadInt32();
                if (fileFormatId != FileFormatId) 
                    return false;
            }
            return true;
        }

        public void Close()
        {
            _reader.Close();
        }

        private readonly Object _fileLock;
        private readonly BinaryReader _reader;

        private readonly double _precursorSignalToNoiseRatioThreshold;
        private readonly double _productSignalToNoiseRatioThreshold;

        private readonly double _minMs1Mz;
        private readonly double _maxMs1Mz;

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

            _reader.BaseStream.Seek(-3 * sizeof (long) - 1 * sizeof (int), SeekOrigin.End);

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

        public override Ms1Spectrum GetMs1Spectrum(int scanNum)
        {
            long offset;
            if (!_scanNumToSpecOffset.TryGetValue(scanNum, out offset)) return null;

            var ms1ScanNums = GetMs1ScanVector();
            var ms1ScanIndex = Array.BinarySearch(ms1ScanNums, scanNum);
            if (ms1ScanIndex < 0) return null;

            lock (_fileLock)
            {
                _reader.BaseStream.Seek(offset, SeekOrigin.Begin);
                while (_reader.BaseStream.Position != (_reader.BaseStream.Length - sizeof(int)))
                {
                    var ms1ScanNum  = _reader.ReadInt32();
                    var msLevel     = _reader.ReadByte();
                    var elutionTime = _reader.ReadDouble();
                    var numPeaks    = _reader.ReadInt32();

                    // load peaks
                    var peaks = new Ms1Peak[numPeaks];
                    for (var i = 0; i < numPeaks; i++)
                    {
                        var mz = _reader.ReadDouble();
                        var intensity = _reader.ReadSingle();
                        peaks[i] = new Ms1Peak(mz, intensity, i) { Ms1SpecIndex = (ushort) ms1ScanIndex };
                    }

                    // Create Ms1Spectrum
                    var ms1Spec = new Ms1Spectrum(scanNum, ms1ScanIndex, peaks)
                    {
                        ElutionTime = elutionTime,
                        MsLevel = msLevel
                    };
                    
                    return ms1Spec;
                }
                return null;
            }
        }

        private Spectrum ReadSpectrum(long offset, bool includePeaks = true)
        {
            lock(_fileLock)
            {
                _reader.BaseStream.Seek(offset, SeekOrigin.Begin);
                while (_reader.BaseStream.Position != (_reader.BaseStream.Length - sizeof (int)))
                {
                    var spec = ReadSpectrum(_reader, includePeaks);
                    return spec;
                }
                return null;
            }
        }

        private IsolationWindow ReadIsolationWindow(long offset)
        {
            lock(_fileLock)
            {
                _reader.BaseStream.Seek(offset, SeekOrigin.Begin);
                while (_reader.BaseStream.Position != (_reader.BaseStream.Length - sizeof(int)))
                {
                    var isolationWindow = ReadIsolationWindow(_reader);
                    return isolationWindow;
                }
                return null;
            }
        }

        private static IsolationWindow ReadIsolationWindow(BinaryReader reader)
        {
            reader.ReadInt32(); // ScanNum
            var msLevel = reader.ReadByte();
            reader.ReadDouble();  // Elution time

            if (msLevel <= 1) return null;

            double? precursorMass = reader.ReadDouble();
            if (precursorMass == 0.0) precursorMass = null;
            int? precursorCharge = reader.ReadInt32();
            if (precursorCharge == 0) precursorCharge = null;
            reader.ReadByte();  // Activation Method
            var isolationWindowTargetMz = reader.ReadDouble();
            var isolationWindowLowerOffset = reader.ReadDouble();
            var isolationWindowUpperOffset = reader.ReadDouble();

            return new IsolationWindow(isolationWindowTargetMz, isolationWindowLowerOffset,
                isolationWindowUpperOffset, precursorMass, precursorCharge);
        }

        private static Spectrum ReadSpectrum(BinaryReader reader, bool includePeaks = true)
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
                var peakList = ReadPeakList(reader, includePeaks);
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
                var peakList = ReadPeakList(reader, includePeaks);
                return new Spectrum(peakList, scanNum)
                {
                    ElutionTime = elutionTime
                };
            }
        }

        public static new void WriteSpectrum(Spectrum spec, BinaryWriter writer)
        {
            // scan number: 4
            writer.Write(spec.ScanNum);

            // ms level: 1
            writer.Write(Convert.ToByte(spec.MsLevel));

            // elution time: 4
            writer.Write(spec.ElutionTime);

            var productSpec = spec as ProductSpectrum;
            if (productSpec != null)    // product spectrum
            {
                var isolationWindow = productSpec.IsolationWindow;
                // precursor mass: 8
                writer.Write(isolationWindow.MonoisotopicMz ?? 0.0);
                // precursor charge: 4
                writer.Write(isolationWindow.Charge ?? 0);
                // Activation method: 1
                writer.Write((byte)productSpec.ActivationMethod);
                // Isolation window target m/z: 8
                writer.Write(productSpec.IsolationWindow.IsolationWindowTargetMz);
                // Isolation window lower offset: 8
                writer.Write(productSpec.IsolationWindow.IsolationWindowLowerOffset);
                // Isolation window uppoer offset: 8
                writer.Write(productSpec.IsolationWindow.IsolationWindowUpperOffset);
            }
            // Number of peaks: 4
            writer.Write(spec.Peaks.Length);
            foreach (var peak in spec.Peaks)
            {
                // m/z: 8
                writer.Write(peak.Mz);
                // intensity: 4
                writer.Write(Convert.ToSingle(peak.Intensity));
            }
        }

        private static List<Peak> ReadPeakList(BinaryReader reader, bool includePeaks = true)
        {
            var peakList = new List<Peak>();
            var numPeaks = reader.ReadInt32();

            // Skip the read if peaks aren't requested
            if (!includePeaks)
            {
                // first the number of peaks, then 12 bytes per peak (mz, double, 8 bytes, then intensity, single, 4 bytes)
                reader.BaseStream.Seek(numPeaks * 12, SeekOrigin.Current);
                return peakList;
            }

            for (var i = 0; i < numPeaks; i++)
            {
                var mz = reader.ReadDouble();
                var intensity = reader.ReadSingle();
                peakList.Add(new Peak(mz, intensity));
            }
            return peakList;
        }


        public static void WriteAsPbf(InMemoryLcMsRun imlr, string outputFilePath, IProgress<ProgressData> progress = null)
        {
            using (var writer = new BinaryWriter(File.Open(outputFilePath, FileMode.Create)))
            {
                WriteAsPbf(imlr, writer, progress);
            }
        }

        public static void WriteAsPbf(InMemoryLcMsRun imlr, BinaryWriter writer, IProgress<ProgressData> progress = null)
        {
            long countTotal = 1;
            long counter = 0;
            if (progress == null)
            {
                progress = new Progress<ProgressData>();
            }
            var progressData = new ProgressData();
            progressData.IsPartialRange = true;
            progressData.Status = "Writing spectra data";

            var scanNumToSpecOffset = new long[imlr.MaxLcScan - imlr.MinLcScan + 1];
            var scanNumToIsolationWindow = new IsolationWindow[imlr.MaxLcScan - imlr.MinLcScan + 1];

            // Spectra
            countTotal = imlr.MaxLcScan - imlr.MinLcScan;
            counter = 0;
            progressData.MaxPercentage = 42.9; // SpecData: Approximately 43% of total file size
            long countMS2Spec = 0;
            for (var scanNum = imlr.MinLcScan; scanNum <= imlr.MaxLcScan; scanNum++)
            {
                progress.Report(progressData.UpdatePercent((double)counter / countTotal * 100.0));
                counter++;
                scanNumToSpecOffset[scanNum - imlr.MinLcScan] = writer.BaseStream.Position;
                var spec = imlr.GetSpectrum(scanNum);
                if (spec == null) continue;
                var productSpec = spec as ProductSpectrum;
                scanNumToIsolationWindow[scanNum - imlr.MinLcScan] = null;
                if (productSpec != null)
                {
                    scanNumToIsolationWindow[scanNum - imlr.MinLcScan] = productSpec.IsolationWindow;
                    countMS2Spec++;
                }
                PbfLcMsRun.WriteSpectrum(spec, writer);
            }

            // Precursor ion chromatograms
            var offsetBeginPrecursorChromatogram = writer.BaseStream.Position;

            var minMzIndex = imlr.Ms1PeakList.Any() ? PbfLcMsRun.GetMzBinIndex(imlr.Ms1PeakList[0].Mz) : 0;
            var maxMzIndex = imlr.Ms1PeakList.Any() ? PbfLcMsRun.GetMzBinIndex(imlr.Ms1PeakList[imlr.Ms1PeakList.Count - 1].Mz) : -1;

            var chromMzIndexToOffset = new long[maxMzIndex - minMzIndex + 1];
            var prevMzIndex = -1;
            counter = 0;
            countTotal = imlr.Ms1PeakList.Count;
            progressData.Status = "Writing precursor chromatograms";
            progressData.StepRange(42.9 + 15.7); // Approximately 16% of total file size
            foreach (var peak in imlr.Ms1PeakList)
            {
                progress.Report(progressData.UpdatePercent((double)counter / countTotal * 100.0));
                counter++;
                var mz = peak.Mz;
                var mzIndex = GetMzBinIndex(mz);
                if (mzIndex > prevMzIndex)
                {
                    chromMzIndexToOffset[mzIndex - minMzIndex] = writer.BaseStream.Position;
                    prevMzIndex = mzIndex;
                }
                writer.Write(peak.Mz);
                writer.Write((float)peak.Intensity);
                writer.Write(peak.ScanNum);
            }

            // Product ion chromatograms
            var ms2PeakList = new List<LcMsPeak>();
            counter = 0;
            countTotal = countMS2Spec;
            progressData.Status = "Processing product chromatograms";
            progressData.StepRange(42.9 + 15.7 + (41.2 / 2)); // Approximately 41% of total file size
            foreach (var ms2ScanNum in imlr.GetScanNumbers(2))
            {
                progress.Report(progressData.UpdatePercent((double)counter / countTotal * 100.0));
                counter++;
                var productSpec = imlr.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                if (productSpec == null) continue;
                foreach (var peak in productSpec.Peaks)
                {
                    ms2PeakList.Add(new LcMsPeak(peak.Mz, peak.Intensity, ms2ScanNum));
                }
            }
            ms2PeakList.Sort();

            var offsetBeginProductChromatogram = writer.BaseStream.Position;
            counter = 0;
            countTotal = ms2PeakList.Count;
            progressData.Status = "Writing product chromatograms";
            progressData.StepRange(42.9 + 15.7 + 41.2); // Approximately 41% of total file size
            foreach (var peak in ms2PeakList)
            {
                progress.Report(progressData.UpdatePercent((double)counter / countTotal * 100.0));
                counter++;
                writer.Write(peak.Mz);
                writer.Write((float)peak.Intensity);
                writer.Write(peak.ScanNum);
            }

            // Meta information
            var offsetBeginMetaInformation = writer.BaseStream.Position;
            progressData.Status = "Writing metadata";
            progressData.IsPartialRange = false;
            progress.Report(progressData.UpdatePercent(99.8)); // Metadata: Approximately 0.2% of total file size

            var warnedInvalidScanNum = false;
            var warnedNullScanToIsolationWindow = false;

            writer.Write(imlr.MinLcScan);
            writer.Write(imlr.MaxLcScan);
            for (var scanNum = imlr.MinLcScan; scanNum <= imlr.MaxLcScan; scanNum++)
            {
                var msLevel = imlr.GetMsLevel(scanNum);
                writer.Write(imlr.GetMsLevel(scanNum));
                writer.Write(imlr.GetElutionTime(scanNum));

                if (msLevel == 2)
                {
                    float minMz = 0;
                    float maxMz = 0;

                    if (scanNum - imlr.MinLcScan < 0 || scanNum - imlr.MinLcScan >= scanNumToIsolationWindow.Length)
                    {
                        if (!warnedInvalidScanNum)
                        {
                            Console.WriteLine("\nWriteAsPbf encountered an invalid scan number: " + scanNum + "; " +
                                              "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown");
                            warnedInvalidScanNum = true;
                        }
                    }
                    else
                    {
                        if (scanNumToIsolationWindow[scanNum - imlr.MinLcScan] == null)
                        {
                            if (!warnedNullScanToIsolationWindow)
                            {
                                Console.WriteLine("\nWriteAsPbf encountered a Null entry in scanNumToIsolationWindow for scan " + scanNum + "; " +
                                                  "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown");
                                warnedNullScanToIsolationWindow = true;
                            }
                        }
                        else
                        {
                            minMz = (float)scanNumToIsolationWindow[scanNum - imlr.MinLcScan].MinMz;
                            maxMz = (float)scanNumToIsolationWindow[scanNum - imlr.MinLcScan].MaxMz;
                        }
                    }

                    writer.Write(minMz);
                    writer.Write(maxMz);
                }
                writer.Write(scanNumToSpecOffset[scanNum - imlr.MinLcScan]);
            }

            // Precursor chromatogram index
            writer.Write(minMzIndex);   // min index
            writer.Write(maxMzIndex);
            progress.Report(progressData.UpdatePercent(99.9)); // Metadata: Approximately 0.2% of total file size

            var prevOffset = offsetBeginMetaInformation;
            for (var i = chromMzIndexToOffset.Length - 1; i >= 0; i--)
            {
                if (chromMzIndexToOffset[i] < offsetBeginPrecursorChromatogram) chromMzIndexToOffset[i] = prevOffset;
                else prevOffset = chromMzIndexToOffset[i];
            }

            foreach (var offset in chromMzIndexToOffset)
            {
                writer.Write(offset);
            }

            writer.Write(offsetBeginPrecursorChromatogram); // 8
            writer.Write(offsetBeginProductChromatogram); // 8
            writer.Write(offsetBeginMetaInformation); // 8
            progress.Report(progressData.UpdatePercent(100.0));
            writer.Write(FileFormatId); // 4
        }

        private long GetOffset(double minMz, double maxMz, long beginOffset, long endOffset)
        {
            var minOffset = beginOffset;
            var maxOffset = endOffset;
            var curOffset = -1L;
            lock(_fileLock)
            {
                // binary search
                while (minOffset <= maxOffset)
                {
                    curOffset = minOffset + (maxOffset - minOffset)/NumBytePeak/2*NumBytePeak;
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
            lock(_fileLock)
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
