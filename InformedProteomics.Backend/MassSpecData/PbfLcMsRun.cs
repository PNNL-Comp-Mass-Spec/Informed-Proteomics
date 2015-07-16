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

        #region Public static functions

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
                    // Fail out if the output path was specified, and we cannot write to it.
                    if (!string.IsNullOrWhiteSpace(pbfFilePath))
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

        /// <summary>
        /// Gets valid possible pbf file paths
        /// </summary>
        /// <param name="specFilePath">Path to the spectra file</param>
        /// <param name="pbfPath">Path to the default pbf file (in the same folder as the spectra file dataset)</param>
        /// <param name="fileName"></param>
        /// <param name="tempPath"></param>
        /// <returns>The default path to the pbf file, unless a valid pbf file exists at the temp path</returns>
        public static string GetCheckPbfFilePath(string specFilePath, out string pbfPath, out string fileName, out string tempPath)
        {
            // Calls "NormalizeDatasetPath" to make sure we save the file to the containing directory
            pbfPath = PbfLcMsRun.GetPbfFileName(MassSpecDataReaderFactory.NormalizeDatasetPath(specFilePath));
            fileName = Path.GetFileName(pbfPath);
            if (String.IsNullOrEmpty(fileName))
            {
                throw new ArgumentException("Cannot create .pbf cache file", "specFilePath");
            }

            tempPath = Path.Combine(Path.GetTempPath(), fileName);
            if (File.Exists(pbfPath) && PbfLcMsRun.CheckFileFormatVersion(pbfPath))
            {
                return pbfPath;
            }

            // Return the temp path if the pbf file of proper format already exists in the temp directory
            if (File.Exists(tempPath) && PbfLcMsRun.CheckFileFormatVersion(tempPath))
            {
                return tempPath;
            }
            return pbfPath;
        }

        #endregion

        #region Constructors

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

            _precursorSignalToNoiseRatioThreshold = precursorSignalToNoiseRatioThreshold;
            _productSignalToNoiseRatioThreshold = productSignalToNoiseRatioThreshold;

            lock (_filelock)
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

        public PbfLcMsRun(string specFileName, IMassSpecDataReader msdr, string pbfFileName = null, 
            double precursorSignalToNoiseRatioThreshold = 0.0, double productSignalToNoiseRatioThreshold = 0.0, IProgress<ProgressData> progress = null)
        {
            string pbfPath2, fileName, tempPath;
            string pbfPath = GetCheckPbfFilePath(specFileName, out pbfPath2, out fileName, out tempPath);
            if (!string.IsNullOrWhiteSpace(pbfFileName))
            {
                pbfPath = pbfFileName;
            }
            if (specFileName.EndsWith(FileExtension) || File.Exists(pbfPath) && CheckFileFormatVersion(pbfPath))
            {
                // Use regular construction
                if (!specFileName.EndsWith(FileExtension))
                {
                    specFileName = pbfPath;
                }
                PbfFilePath = specFileName;
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

                _precursorSignalToNoiseRatioThreshold = precursorSignalToNoiseRatioThreshold;
                _productSignalToNoiseRatioThreshold = productSignalToNoiseRatioThreshold;

                lock (_filelock)
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
                return;
            }

            if (msdr == null)
            {
                msdr = MassSpecDataReaderFactory.GetMassSpecDataReader(specFileName);
            }
            NumSpectra = msdr.NumSpectra;
            _precursorSignalToNoiseRatioThreshold = precursorSignalToNoiseRatioThreshold;
            _productSignalToNoiseRatioThreshold = productSignalToNoiseRatioThreshold;

            try
            {
                _rawFilePath = specFileName;
                PbfFilePath = pbfPath;
                using (var writer =
                        new BinaryWriter(File.Open(pbfPath, FileMode.Create, FileAccess.ReadWrite, FileShare.Read)))
                {
                    WriteToPbf(msdr, writer, progress);
                }
            }
            catch (UnauthorizedAccessException)
            {

                PbfFilePath = tempPath;
                using (var writer =
                        new BinaryWriter(File.Open(tempPath, FileMode.Create, FileAccess.ReadWrite, FileShare.Read)))
                {
                    WriteToPbf(msdr, writer, progress);
                }
            }
            _reader = new BinaryReader(new BufferedStream(File.Open(PbfFilePath, FileMode.Open, FileAccess.Read, FileShare.Read)));

            CreatePrecursorNextScanMap();
        }

        #endregion

        public string PbfFilePath { get; private set; }
        private string _rawFilePath;

        #region IMassSpecDataReader implementation

        public Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true)
        {
            return GetSpectrum(scanNum, includePeaks);
        }

        public bool TryMakeRandomAccessCapable()
        {
            return true;
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

        #endregion

        #region LcMsRun Public functions

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

        public override Ms1Spectrum GetMs1Spectrum(int scanNum)
        {
            long offset;
            if (!_scanNumToSpecOffset.TryGetValue(scanNum, out offset)) return null;

            var ms1ScanNums = GetMs1ScanVector();
            var ms1ScanIndex = Array.BinarySearch(ms1ScanNums, scanNum);
            if (ms1ScanIndex < 0) return null;

            lock (_filelock)
            {
                _reader.BaseStream.Seek(offset, SeekOrigin.Begin);
                while (_reader.BaseStream.Position != (_reader.BaseStream.Length - sizeof(int)))
                {
                    var ms1ScanNum = _reader.ReadInt32();
                    var msLevel = _reader.ReadByte();
                    var elutionTime = _reader.ReadDouble();
                    var numPeaks = _reader.ReadInt32();

                    // load peaks
                    var peaks = new Ms1Peak[numPeaks];
                    for (var i = 0; i < numPeaks; i++)
                    {
                        var mz = _reader.ReadDouble();
                        var intensity = _reader.ReadSingle();
                        peaks[i] = new Ms1Peak(mz, intensity, i) { Ms1SpecIndex = (ushort)ms1ScanIndex };
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

        #endregion

        private readonly object _filelock = new object();
        private BinaryReader _reader;

        private readonly double _precursorSignalToNoiseRatioThreshold;
        private readonly double _productSignalToNoiseRatioThreshold;

        private double _minMs1Mz;
        private double _maxMs1Mz;

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

        private Spectrum ReadSpectrum(long offset, bool includePeaks = true)
        {
            lock (_filelock)
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
            lock (_filelock)
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
                // Isolation window upper offset: 8
                writer.Write(productSpec.IsolationWindow.IsolationWindowUpperOffset);
            }

            // Guarantee sorted peaks.
            var peaks = spec.Peaks.ToList();
            peaks.Sort();
            // Number of peaks: 4
            //writer.Write(spec.Peaks.Length);
            writer.Write(peaks.Count);
            //foreach (var peak in spec.Peaks)
            foreach (var peak in peaks)
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

            var scanNumToSpecOffset = new long[imlr.NumSpectra + 1];
            var scanNumToIsolationWindow = new IsolationWindow[imlr.NumSpectra + 1];

            // Spectra
            countTotal = imlr.NumSpectra;
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

        public static void WriteAsPbf(IMassSpecDataReader msdr, string outputFilePath, IProgress<ProgressData> progress = null)
        {
            using (var writer = new BinaryWriter(File.Open(outputFilePath, FileMode.Create)))
            {
                WriteAsPbf(msdr, writer, progress);
            }
        }

        public static void WriteAsPbf(IMassSpecDataReader msdr, BinaryWriter writer, IProgress<ProgressData> progress = null)
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

            //var scanNumToSpecOffset = new long[msdr.NumSpectra + 1];
            var scanNumToSpecOffset = new Dictionary<int, long>(msdr.NumSpectra + 1);
            //var scanNumToIsolationWindow = new IsolationWindow[msdr.NumSpectra + 1];
            var scanNumToIsolationWindow = new Dictionary<int, IsolationWindow>(msdr.NumSpectra + 1);
            var ms1Scans = new List<int>();
            var ms2Scans = new List<int>();

            // Spectra
            var ms1PeakList = new List<LcMsPeak>();
            var ms2PeakList = new List<LcMsPeak>();
            var scanMetadata = new List<ScanMetadata>(msdr.NumSpectra);
            var minLcScan = int.MaxValue;
            var maxLcScan = int.MinValue;
            countTotal = msdr.NumSpectra;
            counter = 0;
            progressData.MaxPercentage = 42.9; // SpecData: Approximately 43% of total file size
            //long countMS2Spec = 0;
            //for (var scanNum = msdr.MinLcScan; scanNum <= msdr.MaxLcScan; scanNum++)
            foreach (var spec in msdr.ReadAllSpectra())
            {
                progress.Report(progressData.UpdatePercent((double)counter / countTotal * 100.0));
                counter++;
                //scanNumToSpecOffset[scanNum - msdr.MinLcScan] = writer.BaseStream.Position;
                scanNumToSpecOffset.Add(spec.ScanNum, writer.BaseStream.Position);
                //var spec = msdr.GetSpectrum(scanNum);
                //if (spec == null) continue;
                var productSpec = spec as ProductSpectrum;
                //scanNumToIsolationWindow[scanNum - msdr.MinLcScan] = null;
                scanNumToIsolationWindow[spec.ScanNum] = null;
                if (productSpec != null)
                {
                    //scanNumToIsolationWindow[scanNum - msdr.MinLcScan] = productSpec.IsolationWindow;
                    scanNumToIsolationWindow[spec.ScanNum] = productSpec.IsolationWindow;
                    //countMS2Spec++;
                    ms2Scans.Add(productSpec.ScanNum);
                    ms2PeakList.AddRange(productSpec.Peaks.Select(peak => new LcMsPeak(peak.Mz, peak.Intensity, productSpec.ScanNum)));
                }
                else
                {
                    ms1Scans.Add(spec.ScanNum);
                    ms1PeakList.AddRange(spec.Peaks.Select(peak => new LcMsPeak(peak.Mz, peak.Intensity, spec.ScanNum)));
                }
                if (spec.ScanNum < minLcScan)
                {
                    minLcScan = spec.ScanNum;
                }
                if (maxLcScan < spec.ScanNum)
                {
                    maxLcScan = spec.ScanNum;
                }
                scanMetadata.Add(new ScanMetadata(spec.ScanNum, spec.MsLevel, spec.ElutionTime));
                PbfLcMsRun.WriteSpectrum(spec, writer);
            }

            // Precursor ion chromatograms
            ms1PeakList.Sort();
            var offsetBeginPrecursorChromatogram = writer.BaseStream.Position;

            //var minMzIndex = msdr.Ms1PeakList.Any() ? PbfLcMsRun.GetMzBinIndex(msdr.Ms1PeakList[0].Mz) : 0;
            //var maxMzIndex = msdr.Ms1PeakList.Any() ? PbfLcMsRun.GetMzBinIndex(msdr.Ms1PeakList[msdr.Ms1PeakList.Count - 1].Mz) : -1;
            var minMzIndex = ms1PeakList.Any() ? PbfLcMsRun.GetMzBinIndex(ms1PeakList[0].Mz) : 0;
            var maxMzIndex = ms1PeakList.Any() ? PbfLcMsRun.GetMzBinIndex(ms1PeakList[ms1PeakList.Count - 1].Mz) : -1;

            var chromMzIndexToOffset = new long[maxMzIndex - minMzIndex + 1];
            var prevMzIndex = -1;
            counter = 0;
            //countTotal = msdr.Ms1PeakList.Count;
            countTotal = ms1PeakList.Count;
            progressData.Status = "Writing precursor chromatograms";
            progressData.StepRange(42.9 + 15.7); // Approximately 16% of total file size
            //foreach (var peak in msdr.Ms1PeakList)
            foreach (var peak in ms1PeakList)
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
            //var ms2PeakList = new List<LcMsPeak>();
            //counter = 0;
            //countTotal = countMS2Spec;
            //progressData.Status = "Processing product chromatograms";
            //progressData.StepRange(42.9 + 15.7 + (41.2 / 2)); // Approximately 41% of total file size
            //foreach (var ms2ScanNum in msdr.GetScanNumbers(2))
            //{
            //    progress.Report(progressData.UpdatePercent((double)counter / countTotal * 100.0));
            //    counter++;
            //    var productSpec = msdr.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            //    if (productSpec == null) continue;
            //    foreach (var peak in productSpec.Peaks)
            //    {
            //        ms2PeakList.Add(new LcMsPeak(peak.Mz, peak.Intensity, ms2ScanNum));
            //    }
            //}
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

            //writer.Write(msdr.MinLcScan);
            //writer.Write(msdr.MaxLcScan);
            writer.Write(minLcScan);
            writer.Write(maxLcScan);
            scanMetadata.Sort();
            //for (var scanNum = msdr.MinLcScan; scanNum <= msdr.MaxLcScan; scanNum++)
            foreach (var scan in scanMetadata)
            {
                //var msLevel = msdr.GetMsLevel(scanNum);
                //writer.Write(msdr.GetMsLevel(scanNum));
                //writer.Write(msdr.GetElutionTime(scanNum));
                var msLevel = scan.MsLevel;
                writer.Write(scan.MsLevel);
                writer.Write(scan.ElutionTime);

                if (msLevel == 2)
                {
                    float minMz = 0;
                    float maxMz = 0;

                    //if (scanNum - msdr.MinLcScan < 0 || scanNum - msdr.MinLcScan >= scanNumToIsolationWindow.Length
                    if (scan.ScanNum - minLcScan < 0 || scan.ScanNum - minLcScan >= scanNumToIsolationWindow.Count)
                    {
                        if (!warnedInvalidScanNum)
                        {
                            //Console.WriteLine("\nWriteAsPbf encountered an invalid scan number: " + scanNum + "; " +
                            Console.WriteLine("\nWriteAsPbf encountered an invalid scan number: " + scan.ScanNum + "; " +
                                              "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown");
                            warnedInvalidScanNum = true;
                        }
                    }
                    else
                    {
                        //if (scanNumToIsolationWindow[scanNum - msdr.MinLcScan] == null)
                        if (scanNumToIsolationWindow[scan.ScanNum] == null)
                        {
                            if (!warnedNullScanToIsolationWindow)
                            {
                                //Console.WriteLine("\nWriteAsPbf encountered a Null entry in scanNumToIsolationWindow for scan " + scanNum + "; " +
                                Console.WriteLine("\nWriteAsPbf encountered a Null entry in scanNumToIsolationWindow for scan " + scan.ScanNum + "; " +
                                                  "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown");
                                warnedNullScanToIsolationWindow = true;
                            }
                        }
                        else
                        {
                            //minMz = (float)scanNumToIsolationWindow[scanNum - msdr.MinLcScan].MinMz;
                            //maxMz = (float)scanNumToIsolationWindow[scanNum - msdr.MinLcScan].MaxMz;
                            minMz = (float)scanNumToIsolationWindow[scan.ScanNum].MinMz;
                            maxMz = (float)scanNumToIsolationWindow[scan.ScanNum].MaxMz;
                        }
                    }

                    writer.Write(minMz);
                    writer.Write(maxMz);
                }
                //writer.Write(scanNumToSpecOffset[scanNum - msdr.MinLcScan]);
                writer.Write(scanNumToSpecOffset[scan.ScanNum]);
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

        private void WriteToPbf(IMassSpecDataReader msdr, BinaryWriter writer, IProgress<ProgressData> progress = null)
        {
            ScanNumToMsLevel = new Dictionary<int, int>(msdr.NumSpectra + 1);
            ScanNumElutionTimeMap = new Dictionary<int, double>(msdr.NumSpectra + 1);
            IsolationMzBinToScanNums = new Dictionary<int, int[]>();
            _scanNumToSpecOffset = new Dictionary<int, long>(msdr.NumSpectra + 1);

            MinLcScan = int.MaxValue;
            MaxLcScan = int.MinValue;
            MinMsLevel = int.MaxValue;
            MaxMsLevel = int.MinValue;

            long countTotal = 1;
            long counter = 0;
            if (progress == null)
            {
                progress = new Progress<ProgressData>();
            }
            var progressData = new ProgressData();
            progressData.IsPartialRange = true;
            progressData.Status = "Writing spectra data";

            var scanNumToIsolationWindow = new Dictionary<int, IsolationWindow>(msdr.NumSpectra + 1);
            var ms1Scans = new List<int>();
            var ms2Scans = new List<int>();

            // Spectra
            //var ms1PeakList = new List<LcMsPeak>();
            //var ms2PeakList = new List<LcMsPeak>();
            long ms1PeakCount = 0;
            long ms2PeakCount = 0;
            var scanMetadata = new List<ScanMetadata>(msdr.NumSpectra);
            countTotal = msdr.NumSpectra;
            counter = 0;
            progressData.MaxPercentage = 42.9; // SpecData: Approximately 43% of total file size
            foreach (var spec in msdr.ReadAllSpectra())
            {
                progress.Report(progressData.UpdatePercent((double)counter / countTotal * 100.0));
                counter++;
                _scanNumToSpecOffset.Add(spec.ScanNum, writer.BaseStream.Position);
                var productSpec = spec as ProductSpectrum;
                scanNumToIsolationWindow[spec.ScanNum] = null;
                if (productSpec != null)
                {
                    scanNumToIsolationWindow[spec.ScanNum] = productSpec.IsolationWindow;
                    ms2Scans.Add(productSpec.ScanNum);
                    //ms2PeakList.AddRange(productSpec.Peaks.Select(peak => new LcMsPeak(peak.Mz, peak.Intensity, productSpec.ScanNum)));
                    ms2PeakCount += productSpec.Peaks.Length;
                }
                else
                {
                    ms1Scans.Add(spec.ScanNum);
                    //ms1PeakList.AddRange(spec.Peaks.Select(peak => new LcMsPeak(peak.Mz, peak.Intensity, spec.ScanNum)));
                    ms1PeakCount += spec.Peaks.Length;
                }
                if (spec.ScanNum < MinLcScan)
                {
                    MinLcScan = spec.ScanNum;
                }
                if (MaxLcScan < spec.ScanNum)
                {
                    MaxLcScan = spec.ScanNum;
                }
                scanMetadata.Add(new ScanMetadata(spec.ScanNum, spec.MsLevel, spec.ElutionTime));
                PbfLcMsRun.WriteSpectrum(spec, writer);
            }

            ms1Scans.Sort();
            ms2Scans.Sort();

            // Precursor ion chromatograms
            _offsetPrecursorChromatogramStart = writer.BaseStream.Position;
            progressData.Status = "Writing precursor chromatograms";
            progressData.StepRange(42.9 + 15.7); // Approximately 16% of total file size

            var prog = new Progress<ProgressData>(p =>
            {
                progressData.StatusInternal = p.Status;
                progress.Report(progressData.UpdatePercent(p.Percent));
            });
            //CreateAndOutputMsXChromatogram(writer, ms1Scans, ms1PeakList, true, prog);
            //CreateAndOutputMsXChromatogram(writer, ms1Scans, null, true, prog);
            CreateAndOutputMsXChromatogram_Merge(writer, ms1Scans, ms1PeakCount, null, true, prog);

            // Product ion chromatograms
            _offsetProductChromatogramBegin = writer.BaseStream.Position;
            _offsetPrecursorChromatogramEnd = _offsetProductChromatogramBegin;
            progressData.Status = "Writing product chromatograms";
            progressData.StepRange(42.9 + 15.7 + 41.2); // Approximately 41% of total file size

            prog = new Progress<ProgressData>(p =>
            {
                progressData.StatusInternal = p.Status;
                progress.Report(progressData.UpdatePercent(p.Percent));
            });
            //CreateAndOutputMsXChromatogram(writer, ms2Scans, ms2PeakList, false, prog);
            //CreateAndOutputMsXChromatogram(writer, ms2Scans, null, false, prog);
            CreateAndOutputMsXChromatogram_Merge(writer, ms2Scans, ms2PeakCount, null, false, prog);

            // Meta information
            _offsetProductChromatogramEnd = writer.BaseStream.Position;
            var offsetBeginMetaInformation = _offsetProductChromatogramEnd;
            progressData.Status = "Writing metadata";
            progressData.IsPartialRange = false;
            progress.Report(progressData.UpdatePercent(99.8)); // Metadata: Approximately 0.2% of total file size

            var warnedInvalidScanNum = false;
            var warnedNullScanToIsolationWindow = false;

            writer.Write(MinLcScan);
            writer.Write(MaxLcScan);
            scanMetadata.Sort();

            var isoWindowSet = new HashSet<IsolationWindow>();
            var isDda = false;
            var isolationMzBinToScanNums = new Dictionary<int, List<int>>();

            MinMsLevel = scanMetadata.Select(scan => scan.MsLevel).Min();
            MaxMsLevel = scanMetadata.Select(scan => scan.MsLevel).Max();
            foreach (var scan in scanMetadata)
            {
                var msLevel = scan.MsLevel;
                writer.Write(scan.MsLevel);
                writer.Write(scan.ElutionTime);
                ScanNumToMsLevel[scan.ScanNum] = scan.MsLevel;
                ScanNumElutionTimeMap[scan.ScanNum] = scan.ElutionTime;

                if (msLevel == 2)
                {
                    float minMz = 0;
                    float maxMz = 0;

                    if (scan.ScanNum - MinLcScan < 0 || scan.ScanNum - MinLcScan >= scanNumToIsolationWindow.Count)
                    {
                        if (!warnedInvalidScanNum)
                        {
                            Console.WriteLine("\nWriteAsPbf encountered an invalid scan number: " + scan.ScanNum + "; " +
                                              "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown");
                            warnedInvalidScanNum = true;
                        }
                    }
                    else
                    {
                        if (scanNumToIsolationWindow[scan.ScanNum] == null)
                        {
                            if (!warnedNullScanToIsolationWindow)
                            {
                                Console.WriteLine("\nWriteAsPbf encountered a Null entry in scanNumToIsolationWindow for scan " + scan.ScanNum + "; " +
                                                  "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown");
                                warnedNullScanToIsolationWindow = true;
                            }
                        }
                        else
                        {
                            minMz = (float)scanNumToIsolationWindow[scan.ScanNum].MinMz;
                            maxMz = (float)scanNumToIsolationWindow[scan.ScanNum].MaxMz;
                        }
                    }

                    writer.Write(minMz);
                    writer.Write(maxMz);

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
                        scanNumList.Add(scan.ScanNum);
                    }
                }
                writer.Write(_scanNumToSpecOffset[scan.ScanNum]);
            }

            IsDiaOrNull = !isDda;

            foreach (var entry in isolationMzBinToScanNums)
            {
                var binNum = entry.Key;
                entry.Value.Sort();
                var scanNumList = entry.Value.ToArray();
                IsolationMzBinToScanNums[binNum] = scanNumList;
            }

            // Precursor chromatogram index
            writer.Write(_minMzIndex);   // min index
            writer.Write(_maxMzIndex);
            progress.Report(progressData.UpdatePercent(99.9)); // Metadata: Approximately 0.2% of total file size

            var prevOffset = offsetBeginMetaInformation;
            for (var i = _chromMzIndexToOffset.Length - 2; i >= 0; i--)
            {
                if (_chromMzIndexToOffset[i] < _offsetPrecursorChromatogramStart)
                    _chromMzIndexToOffset[i] = prevOffset;
                else
                    prevOffset = _chromMzIndexToOffset[i];
            }

            foreach (var offset in _chromMzIndexToOffset.Take(_chromMzIndexToOffset.Length - 1))
            {
                writer.Write(offset);
            }

            _chromMzIndexToOffset[_chromMzIndexToOffset.Length - 1] = _offsetPrecursorChromatogramEnd;

            writer.Write(_offsetPrecursorChromatogramStart); // 8
            writer.Write(_offsetProductChromatogramBegin); // 8
            writer.Write(offsetBeginMetaInformation); // 8
            progress.Report(progressData.UpdatePercent(100.0));
            writer.Write(FileFormatId); // 4
        }

        private void CreateAndOutputMsXChromatogram(BinaryWriter writer, List<int> scansForMsLevelX, List<LcMsPeak> peaksList, bool isMs1List, IProgress<ProgressData> progress)
        {
            // Other thought for a slower, but really low-memory chromatogram creator:
            //   Make sure peaks are sorted ascending when written, and then jump through all spectra, performing a massive merge sort on the peaks and outputting the lowest spectra
            var progData = new ProgressData();
            var count = 0;
            var countTotal = scansForMsLevelX.Count;
            if (peaksList == null || peaksList.Count == 0)
            {
                progData.MaxPercentage = 50;
                peaksList = new List<LcMsPeak>();
                if (_reader == null)
                {
                    _reader = new BinaryReader(new BufferedStream(File.Open(PbfFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)));
                }
                count = 0;
                progData.Status = "Reading MSn peaks for product chromatogram";
                if (isMs1List)
                {
                    progData.Status = "Reading MS1 peaks for precursor chromatogram";
                }
                foreach (var scan in scansForMsLevelX)
                {
                    progress.Report(progData.UpdatePercent((double)count / countTotal * 100));
                    count++;
                    peaksList.AddRange(GetPeaksForSpectrum(scan));
                }
                progData.StepRange(100);
            }
            peaksList.Sort();

            if (isMs1List)
            {
                _minMs1Mz = peaksList.First().Mz;
                _maxMs1Mz = peaksList.Last().Mz;

                _minMzIndex = peaksList.Any() ? PbfLcMsRun.GetMzBinIndex(peaksList[0].Mz) : 0;
                _maxMzIndex = peaksList.Any() ? PbfLcMsRun.GetMzBinIndex(peaksList[peaksList.Count - 1].Mz) : -1;

                _chromMzIndexToOffset = new long[_maxMzIndex - _minMzIndex + 2];
            }

            count = 0;
            countTotal = peaksList.Count;
            var prevMzIndex = -1;
            progData.Status = "Writing product chromatogram";
            if (isMs1List)
            {
                progData.Status = "Writing precursor chromatogram";
            }
            foreach (var peak in peaksList)
            {
                progress.Report(progData.UpdatePercent((double)count / countTotal * 100.0));
                count++;
                if (isMs1List)
                {
                    var mzIndex = GetMzBinIndex(peak.Mz);
                    if (mzIndex > prevMzIndex)
                    {
                        _chromMzIndexToOffset[mzIndex - _minMzIndex] = writer.BaseStream.Position;
                        prevMzIndex = mzIndex;
                    }
                }
                writer.Write(peak.Mz);
                writer.Write((float)peak.Intensity);
                writer.Write(peak.ScanNum);
            }
        }

        private void CreateAndOutputMsXChromatogram_Merge(BinaryWriter writer, List<int> scansForMsLevelX, long totalPeaksCount, List<LcMsPeak> peaksList, bool isMs1List, IProgress<ProgressData> progress)
        {
            // Other thought for a slower, but really low-memory chromatogram creator:
            //   Make sure peaks are sorted ascending when written, and then jump through all spectra, performing a massive merge sort on the peaks and outputting the lowest spectra
            var progData = new ProgressData();
            var count = 0;
            var countTotal = scansForMsLevelX.Count;
            var prevMzIndex = -1;

            //var avgPeaksPerSpec = (double)totalPeaksCount / scansForMsLevelX.Count;
            var totalKBytesInMem = (totalPeaksCount * (20 + 16 + 20)) / 1024;
            double memoryFreeKB = 0;
            // For pre-Vista: "SELECT * FROM Win32_LogicalMemoryConfiguration", and a different property.
            foreach (var item in new System.Management.ManagementObjectSearcher("SELECT * FROM CIM_OperatingSystem").Get())
            {
                memoryFreeKB += int.Parse(item["FreePhysicalMemory"].ToString());
                //foreach (var p in item.Properties)
                //{
                //    Console.WriteLine("{0}: {1}", p.Name, p.Value);
                //}
            }
            Console.WriteLine("memoryFreeKB: " + memoryFreeKB);
            Console.WriteLine("totalKBytesInMem: " + totalKBytesInMem);

            //if ((peaksList == null || peaksList.Count == 0) && (totalKBytesInMem > memoryFreeKB / 4))
            if ((peaksList == null || peaksList.Count == 0))
            {
                // Perform a massive merge-sort style read/write, to minimize memory usage - set a minimum of 5
                int maxInMemoryPerSpec = (int)(memoryFreeKB * 1024 / 2 / scansForMsLevelX.Count / (20 + 16 + 20));
                if (maxInMemoryPerSpec < 5)
                    maxInMemoryPerSpec = 5;

                Console.WriteLine("maxInMemoryPerSpec: " + maxInMemoryPerSpec);

                progData.Status = "Writing product chromatogram";
                if (isMs1List)
                {
                    progData.Status = "Writing precursor chromatogram";
                }
                var peaks = new SortedSet<LcMsPeak>();
                var peaksCount = 0;
                var metadata = new Dictionary<int, ScanPeakMetaData>(scansForMsLevelX.Count);
                lock (_filelock)
                {
                    if (_reader == null)
                    {
                        _reader =
                            new BinaryReader(
                                new BufferedStream(File.Open(PbfFilePath, FileMode.Open, FileAccess.Read,
                                    FileShare.ReadWrite)));
                    }

                    progData.MaxPercentage = isMs1List ? 10 : 20;
                    // Read in metadata, and the first 5 peaks of each scan
                    foreach (var scan in scansForMsLevelX)
                    {
                        progress.Report(progData.UpdatePercent((double) count / countTotal * 100));
                        count++;
                        var datum = GetPeakMetaDataForSpectrum(scan);
                        peaksCount += datum.NumPeaks;
                        metadata.Add(datum.ScanNum, datum);
                        foreach (var peak in datum.ReadPeaks(_reader, maxInMemoryPerSpec))
                        {
                            peaks.Add(peak);
                            datum.Count++;
                        }
                    }

                    if (peaksCount == 0)
                    {
                        if (isMs1List)
                        {
                            _minMs1Mz = int.MaxValue;
                            _maxMs1Mz = int.MinValue;

                            _minMzIndex = 0;
                            _maxMzIndex = -1;
                            _chromMzIndexToOffset = new long[0];
                        }
                        return;
                    }

                    progData.StepRange(100);
                    progress.Report(progData.UpdatePercent(0));
                    var peaksToRemove = new List<LcMsPeak>();
                    count = 0;
                    prevMzIndex = -1;
                    if (isMs1List)
                    {
                        _minMs1Mz = peaks.Min.Mz;
                        _minMzIndex = GetMzBinIndex(_minMs1Mz);
                    }
                    var chromMzIndexToOffset = new Dictionary<int, long>();

                    while (peaks.Count > 0)
                    {
                        foreach (var peak in peaks)
                        {
                            progress.Report(progData.UpdatePercent((double) count / peaksCount * 100.0));
                            count++;
                            if (isMs1List)
                            {
                                _maxMs1Mz = peak.Mz;
                                var mzIndex = GetMzBinIndex(peak.Mz);
                                if (mzIndex > prevMzIndex)
                                {
                                    chromMzIndexToOffset.Add(mzIndex - _minMzIndex, writer.BaseStream.Position);
                                    prevMzIndex = mzIndex;
                                }
                            }
                            writer.Write(peak.Mz);
                            writer.Write((float) peak.Intensity);
                            writer.Write(peak.ScanNum);

                            peaksToRemove.Add(peak);

                            var datum = metadata[peak.ScanNum];
                            datum.Count--;

                            if (datum.Count < 1 && datum.MorePeaksToRead)
                            {
                                break;
                            }
                        }

                        // Remove the output peaks from the sorted set
                        foreach (var peak in peaksToRemove)
                        {
                            peaks.Remove(peak);
                        }

                        peaksToRemove.Clear();

                        // add new entries back into the list from the spectra that the peaks came from
                        foreach (var datum in metadata.Values)
                        {
                            if (datum.Count < maxInMemoryPerSpec && datum.MorePeaksToRead)
                            {
                                foreach (var peak in datum.ReadPeaks(_reader, maxInMemoryPerSpec - datum.Count))
                                {
                                    peaks.Add(peak);
                                    datum.Count++;
                                }
                            }
                        }
                    }

                    if (isMs1List)
                    {
                        //_minMzIndex = GetMzBinIndex(_minMs1Mz);
                        _maxMzIndex = GetMzBinIndex(_maxMs1Mz);

                        _chromMzIndexToOffset = new long[_maxMzIndex - _minMzIndex + 2];
                        for (int i = 0; i < _chromMzIndexToOffset.Length; i++)
                        {
                            if (chromMzIndexToOffset.ContainsKey(i))
                            {
                                _chromMzIndexToOffset[i] = chromMzIndexToOffset[i];
                            }
                        }
                    }
                }
                return;
            }

            if (peaksList == null || peaksList.Count == 0)
            {
                progData.MaxPercentage = 50;
                peaksList = new List<LcMsPeak>();
                if (_reader == null)
                {
                    _reader = new BinaryReader(new BufferedStream(File.Open(PbfFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)));
                }
                count = 0;
                progData.Status = "Reading MSn peaks for product chromatogram";
                if (isMs1List)
                {
                    progData.Status = "Reading MS1 peaks for precursor chromatogram";
                }
                foreach (var scan in scansForMsLevelX)
                {
                    progress.Report(progData.UpdatePercent((double)count / countTotal * 100));
                    count++;
                    peaksList.AddRange(GetPeaksForSpectrum(scan));
                }
                progData.StepRange(100);
            }
            peaksList.Sort();

            if (isMs1List)
            {
                _minMs1Mz = peaksList.First().Mz;
                _maxMs1Mz = peaksList.Last().Mz;

                _minMzIndex = peaksList.Any() ? PbfLcMsRun.GetMzBinIndex(peaksList[0].Mz) : 0;
                _maxMzIndex = peaksList.Any() ? PbfLcMsRun.GetMzBinIndex(peaksList[peaksList.Count - 1].Mz) : -1;

                _chromMzIndexToOffset = new long[_maxMzIndex - _minMzIndex + 2];
            }

            count = 0;
            countTotal = peaksList.Count;
            prevMzIndex = -1;
            progData.Status = "Writing product chromatogram";
            if (isMs1List)
            {
                progData.Status = "Writing precursor chromatogram";
            }
            foreach (var peak in peaksList)
            {
                progress.Report(progData.UpdatePercent((double)count / countTotal * 100.0));
                count++;
                if (isMs1List)
                {
                    var mzIndex = GetMzBinIndex(peak.Mz);
                    if (mzIndex > prevMzIndex)
                    {
                        _chromMzIndexToOffset[mzIndex - _minMzIndex] = writer.BaseStream.Position;
                        prevMzIndex = mzIndex;
                    }
                }
                writer.Write(peak.Mz);
                writer.Write((float)peak.Intensity);
                writer.Write(peak.ScanNum);
            }
        }

        private List<LcMsPeak> GetPeaksForSpectrum(int scanNum)
        {
            long offset;
            if (!_scanNumToSpecOffset.TryGetValue(scanNum, out offset))
            {
                // Empty list won't cause exception, but null will.
                return new List<LcMsPeak>();
            }

            var peakList = new List<LcMsPeak>();
            lock (_filelock)
            {
                _reader.BaseStream.Seek(offset, SeekOrigin.Begin);

                var rScanNum = _reader.ReadInt32();
                var msLevel = _reader.ReadByte();
                var elutionTime = _reader.ReadDouble();

                if (msLevel > 1)
                {
                    _reader.BaseStream.Seek(8 + 4 + 1 + 8 + 8 + 8, SeekOrigin.Current);
                    //double? precursorMass = _reader.ReadDouble();
                    //int? precursorCharge = _reader.ReadInt32();
                    //var activationMethod = (ActivationMethod) _reader.ReadByte();
                    //var isolationWindowTargetMz = _reader.ReadDouble();
                    //var isolationWindowLowerOffset = _reader.ReadDouble();
                    //var isolationWindowUpperOffset = _reader.ReadDouble();
                }

                var numPeaks = _reader.ReadInt32();

                for (var i = 0; i < numPeaks; i++)
                {
                    var mz = _reader.ReadDouble();
                    var intensity = _reader.ReadSingle();
                    peakList.Add(new LcMsPeak(mz, intensity, rScanNum));
                }
            }
            return peakList;
        }

        private ScanPeakMetaData GetPeakMetaDataForSpectrum(int scanNum)
        {
            long offset;
            if (!_scanNumToSpecOffset.TryGetValue(scanNum, out offset))
            {
                return null;
            }

            ScanPeakMetaData data;
            lock (_filelock)
            {
                _reader.BaseStream.Seek(offset, SeekOrigin.Begin);

                var rScanNum = _reader.ReadInt32();
                var msLevel = _reader.ReadByte();
                var elutionTime = _reader.ReadDouble();

                if (msLevel > 1)
                {
                    _reader.BaseStream.Seek(8 + 4 + 1 + 8 + 8 + 8, SeekOrigin.Current);
                    //double? precursorMass = _reader.ReadDouble();
                    //int? precursorCharge = _reader.ReadInt32();
                    //var activationMethod = (ActivationMethod) _reader.ReadByte();
                    //var isolationWindowTargetMz = _reader.ReadDouble();
                    //var isolationWindowLowerOffset = _reader.ReadDouble();
                    //var isolationWindowUpperOffset = _reader.ReadDouble();
                }

                data = new ScanPeakMetaData(rScanNum, _reader);
            }
            return data;
        }

        private class ScanMetadata : IComparable<ScanMetadata>
        {
            public int ScanNum { get; private set; }
            public int MsLevel { get; private set; }
            public double ElutionTime { get; private set; }

            public ScanMetadata(int scanTime, int msLevel, double elutionTime)
            {
                ScanNum = scanTime;
                MsLevel = msLevel;
                ElutionTime = elutionTime;
            }

            public int CompareTo(ScanMetadata other)
            {
                if (ScanNum.CompareTo(other.ScanNum) == 0)
                {
                    if (MsLevel.CompareTo(other.MsLevel) == 0)
                    {
                        return ElutionTime.CompareTo(other.ElutionTime);
                    }
                    return MsLevel.CompareTo(other.MsLevel);
                }
                return ScanNum.CompareTo(other.ScanNum);
            }
        }

        private class ScanPeakMetaData
        {
            public int ScanNum { get; private set; }
            public int NumPeaks { get; private set; }
            public int NumPeaksRead { get; private set; }
            public long NextPeakOffset { get; private set; }
            public int Count { get; set; }

            public bool MorePeaksToRead
            {
                get { return NumPeaksRead < NumPeaks; }
            }

            public ScanPeakMetaData(int scanNum, BinaryReader reader)
            {
                Count = 0;
                ScanNum = scanNum;
                NumPeaks = reader.ReadInt32();
                NumPeaksRead = 0;
                NextPeakOffset = reader.BaseStream.Position;
            }

            public IEnumerable<LcMsPeak> ReadPeaks(BinaryReader reader, int numPeaksToRead)
            {
                var peaks = new List<LcMsPeak>();
                if (NumPeaksRead == NumPeaks)
                {
                    return peaks;
                }
                reader.BaseStream.Seek(NextPeakOffset, SeekOrigin.Begin);

                for (int i = 0; i < numPeaksToRead && NumPeaksRead < NumPeaks; i++, NumPeaksRead++)
                {
                    var mz = reader.ReadDouble();
                    var intensity = reader.ReadSingle();
                    peaks.Add(new LcMsPeak(mz, intensity, ScanNum));
                }
                NextPeakOffset = reader.BaseStream.Position;
                return peaks;
            }
        }

        private long GetOffset(double minMz, double maxMz, long beginOffset, long endOffset)
        {
            var minOffset = beginOffset;
            var maxOffset = endOffset;
            var curOffset = -1L;
            lock (_filelock)
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
            lock (_filelock)
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
