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
        public const string FileExtension = ".pbf";

        // This constant should be incremented by 1 if the binary file format is changed
        public const int FileFormatId = 150605;
        private const int EarliestSupportedFileFormatId = 150604;

        /** File format id history
         * 150604: Earliest supported, has all data needed by InformedProteomics projects
         * 150605: Added Total Ion Current and Native ID fields to spectrum output.
         * 
         */

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
            var specReader = MassSpecDataReaderFactory.GetMassSpecDataReader(specFilePath);
            return GetLcMsRun(specFilePath, specReader, precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold, progress);
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
            // ReSharper disable once CanBeReplacedWithTryCastAndCheckForNull (being silly)
            if (specReader is PbfLcMsRun)
                return (LcMsRun)specReader;

            return new PbfLcMsRun(specFilePath, specReader, "", precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold, progress);
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
        [Obsolete("Use GetLcMsRun() for an optimized pbf creation process", false)]
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
        [Obsolete("Use GetLcMsRun() for an optimized pbf creation process", false)]
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
            var progData = new ProgressData(progress);
            progData.StepRange(75.0);
            if (progress != null)
            {
                prog = new Progress<ProgressData>(p =>
                {
                    progData.Status = p.Status;
                    progData.Report(p.Percent);
                });
            }

            bool isCurrent;
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
                if (File.Exists(tempPath) && PbfLcMsRun.CheckFileFormatVersion(tempPath, out isCurrent) && isCurrent)
                {
                    return tempPath;
                }
            }

            if (!File.Exists(pbfPath) || !(PbfLcMsRun.CheckFileFormatVersion(pbfPath, out isCurrent) && isCurrent))
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
                    if (!File.Exists(tempPath) || !(PbfLcMsRun.CheckFileFormatVersion(tempPath, out isCurrent) && isCurrent))
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
            bool isCurrent;
            if (File.Exists(pbfPath) && PbfLcMsRun.CheckFileFormatVersion(pbfPath, out isCurrent))
            {
                return pbfPath;
            }

            // Return the temp path if the pbf file of proper format already exists in the temp directory
            if (File.Exists(tempPath) && PbfLcMsRun.CheckFileFormatVersion(tempPath, out isCurrent))
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

            _reader = new BinaryReader(File.Open(specFileName, FileMode.Open, FileAccess.Read, FileShare.Read));

            _precursorSignalToNoiseRatioThreshold = precursorSignalToNoiseRatioThreshold;
            _productSignalToNoiseRatioThreshold = productSignalToNoiseRatioThreshold;

            lock (_filelock)
            {
                if (ReadMetaInfo() == false)
                {
                    throw new FormatException("Illegal pbf file format!");
                }

                if (_offsetPrecursorChromatogramBegin == _offsetPrecursorChromatogramEnd)
                {
                    // No MS1 data
                    _minMs1Mz = 0;
                    _maxMs1Mz = 0;
                }
                else
                {
                    _reader.BaseStream.Seek(_offsetPrecursorChromatogramBegin, SeekOrigin.Begin);
                    _minMs1Mz = _reader.ReadDouble();

                    if (_offsetPrecursorChromatogramEnd - NumBytePeak < 0)
                    {
                        throw new FormatException("Corrupt pbf file (_offsetPrecursorChromatogramEnd is < 0)");
                    }

                    if (_offsetPrecursorChromatogramEnd - NumBytePeak >= specFile.Length)
                    {
                        throw new FormatException(
                            "Corrupt pbf file (_offsetPrecursorChromatogramEnd is past the end of the file)");
                    }

                    _reader.BaseStream.Seek(_offsetPrecursorChromatogramEnd - NumBytePeak, SeekOrigin.Begin);
                    _maxMs1Mz = _reader.ReadDouble();
                }
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
            bool isCurrent;
            if (specFileName.EndsWith(FileExtension) || File.Exists(pbfPath) && CheckFileFormatVersion(pbfPath, out isCurrent) && isCurrent)
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

                _reader = new BinaryReader(File.Open(specFileName, FileMode.Open, FileAccess.Read, FileShare.Read));

                _precursorSignalToNoiseRatioThreshold = precursorSignalToNoiseRatioThreshold;
                _productSignalToNoiseRatioThreshold = productSignalToNoiseRatioThreshold;

                lock (_filelock)
                {
                    if (ReadMetaInfo() == false)
                    {
                        throw new FormatException("Illegal pbf file format!");
                    }
                    _reader.BaseStream.Seek(_offsetPrecursorChromatogramBegin, SeekOrigin.Begin);
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
            _reader = new BinaryReader(File.Open(PbfFilePath, FileMode.Open, FileAccess.Read, FileShare.Read));

            CreatePrecursorNextScanMap();
        }

        #endregion

        #region Member Variables

        public string PbfFilePath { get; private set; }
        private string _rawFilePath;
        private int _fileFormatId = FileFormatId; // For internal checks and backwards compatibility usages.
        private const int NativeIdLength = 50;

        private readonly object _filelock = new object();
        private BinaryReader _reader;

        private readonly double _precursorSignalToNoiseRatioThreshold;
        private readonly double _productSignalToNoiseRatioThreshold;

        private double _minMs1Mz;
        private double _maxMs1Mz;

        private long _offsetPrecursorChromatogramBegin;
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

        // Each peak is a double, a float, and an int, representing mass, intensity, and scan number
        private const int NumBytePeak = 16;

        private List<XicPoint> _precursorChromatogramCache = new List<XicPoint>();

        #endregion

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

        public static bool CheckFileFormatVersion(string filePath, out bool isCurrent)
        {
            isCurrent = false;
            var pbfFile = new FileInfo(filePath);
            if (!pbfFile.Exists || pbfFile.Length < sizeof(int))
                return false;

            var fs = new FileStream(pbfFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            using (var reader = new BinaryReader(fs))
            {
                fs.Seek(-1 * sizeof(int), SeekOrigin.End);

                var fileFormatId = reader.ReadInt32();
                if (fileFormatId > FileFormatId && fileFormatId < EarliestSupportedFileFormatId)
                    return false;
                if (fileFormatId == FileFormatId)
                    isCurrent = true;
            }
            return true;
        }

        public void Close()
        {
            _reader.Close();
        }

        #endregion

        #region LcMsRun Public function overrides

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
            if (_scanNumToSpecOffset.TryGetValue(scanNum, out offset))
            {
                var spec = ReadSpectrum(offset, false) as ProductSpectrum;
                if (spec != null)
                {
                    return spec.IsolationWindow;
                }
            }
            return null;
        }

        /// <summary>
        /// Returns a xic for the chosen range that covers the entire run.
        /// </summary>
        /// <param name="minMz"></param>
        /// <param name="maxMz"></param>
        /// <param name="precursorMz"></param>
        /// <returns></returns>
        public override Xic GetFullProductExtractedIonChromatogram(double minMz, double maxMz, double precursorMz)
        {
            var targetOffset = GetOffset(minMz, maxMz, _offsetProductChromatogramBegin, _offsetProductChromatogramEnd);
            if (targetOffset < _offsetProductChromatogramBegin)
            {
                return new Xic();
            }
            var xic = GetXicPointsWithin(minMz, maxMz, _offsetProductChromatogramBegin, _offsetProductChromatogramEnd, targetOffset);
            if (!xic.Any())
            {
                return xic;
            }

            var scanToXicPoint = new XicPoint[MaxLcScan - MinLcScan + 1];
            foreach (var xicPoint in xic)
            {
                var prev = scanToXicPoint[xicPoint.ScanNum - MinLcScan];
                if (prev == null || xicPoint.Intensity > prev.Intensity)
                {
                    scanToXicPoint[xicPoint.ScanNum - MinLcScan] = xicPoint;
                }
            }

            var newXic = new Xic();

            newXic.AddRange(GetFragmentationSpectraScanNums(precursorMz).Select(scanNum => scanToXicPoint[scanNum - MinLcScan] ?? new XicPoint(scanNum, 0, 0)));
            return newXic;
        }

        /// <summary>
        /// Returns selected peaks between minMz and maxMz. The biggest peak per scan is selected.
        /// </summary>
        /// <param name="minMz"></param>
        /// <param name="maxMz"></param>
        /// <returns></returns>
        public override Xic GetPrecursorExtractedIonChromatogram(double minMz, double maxMz)
        {
            if (_precursorChromatogramCache.Count > 0 && _precursorChromatogramCache.First().Mz < minMz && _precursorChromatogramCache.Last().Mz > maxMz)
            {
                var xicl = new Xic();
                xicl.AddRange(_precursorChromatogramCache.Where(peak => minMz <= peak.Mz && peak.Mz <= maxMz));
                return Xic.GetSelectedXic(xicl);
            }

            var minBinIndex = GetMzBinIndex(minMz);
            var maxBinIndex = GetMzBinIndex(maxMz);

            long targetOffset;
            if (minBinIndex == maxBinIndex)
            {
                if (maxBinIndex < _minMzIndex || maxBinIndex > _maxMzIndex) return new Xic();
                var offset = _chromMzIndexToOffset[maxBinIndex - _minMzIndex];
                if (offset < _offsetPrecursorChromatogramBegin) return new Xic();

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

            if (targetOffset < _offsetPrecursorChromatogramBegin) return new Xic();
            var xic = GetXic(minMz, maxMz, _offsetPrecursorChromatogramBegin, _offsetPrecursorChromatogramEnd, targetOffset);
            //if (!xic.Any()) return xic;
            return xic;
        }

        /// <summary>
        /// Returns all peaks between minMz and maxMz, including multiple peaks per scan
        /// </summary>
        /// <param name="minMz"></param>
        /// <param name="maxMz"></param>
        /// <returns></returns>
        public override Xic GetPrecursorChromatogramRange(double minMz, double maxMz)
        {
            if (_precursorChromatogramCache.Count > 0 && _precursorChromatogramCache.First().Mz < minMz && _precursorChromatogramCache.Last().Mz > maxMz)
            {
                var xicl = new Xic();
                xicl.AddRange(_precursorChromatogramCache.Where(peak => minMz <= peak.Mz && peak.Mz <= maxMz));
                return xicl;
            }

            var minBinIndex = GetMzBinIndex(minMz);
            var maxBinIndex = GetMzBinIndex(maxMz);

            long targetOffset;
            if (minBinIndex == maxBinIndex)
            {
                if (maxBinIndex < _minMzIndex || maxBinIndex > _maxMzIndex) return new Xic();
                var offset = _chromMzIndexToOffset[maxBinIndex - _minMzIndex];
                if (offset < _offsetPrecursorChromatogramBegin) return new Xic();

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

            if (targetOffset < _offsetPrecursorChromatogramBegin) return new Xic();
            var xic = GetChromatogramRange(minMz, maxMz, _offsetPrecursorChromatogramBegin, _offsetPrecursorChromatogramEnd, targetOffset);
            if (!xic.Any()) return xic;
            xic.Sort();
            return xic;
        }

        public override Ms1Spectrum GetMs1Spectrum(int scanNum)
        {
            long offset;
            if (!_scanNumToSpecOffset.TryGetValue(scanNum, out offset)) return null;

            var ms1ScanNums = GetMs1ScanVector();
            var ms1ScanIndex = Array.BinarySearch(ms1ScanNums, scanNum);
            if (ms1ScanIndex < 0) return null;

            return ReadSpectrum(offset, true, ms1ScanIndex) as Ms1Spectrum;
        }

        #endregion

        #region Metadata reading

        private bool ReadMetaInfo()
        {
            ScanNumToMsLevel = new Dictionary<int, int>();
            ScanNumElutionTimeMap = new Dictionary<int, double>();
            IsolationMzBinToScanNums = new Dictionary<int, int[]>();

            // Read the file format integer from the end of the file
            _reader.BaseStream.Seek(-1 * sizeof(int), SeekOrigin.End);
            _fileFormatId = _reader.ReadInt32();
            if (_fileFormatId > FileFormatId || _fileFormatId < EarliestSupportedFileFormatId)
                return false;

            // Backup 10 bytes
            _reader.BaseStream.Seek(-3 * sizeof (long) - 1 * sizeof (int), SeekOrigin.End);

            // Read the byte offset of the start of the precursor chromatogram
            _offsetPrecursorChromatogramBegin = _reader.ReadInt64();

            // Read the byte offset of the start of the product chromatogram (which is the end ofthe precursor chromatogram)
            _offsetPrecursorChromatogramEnd = _offsetProductChromatogramBegin = _reader.ReadInt64();

            // Read the byte offset of the end of the product chromatogram
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

        #endregion

        #region Read/Write Spectrum

        private Spectrum ReadSpectrum(long offset, bool includePeaks = true, int specIndex = -1)
        {
            lock (_filelock)
            {
                _reader.BaseStream.Seek(offset, SeekOrigin.Begin);
                while (_reader.BaseStream.Position != (_reader.BaseStream.Length - sizeof(int)))
                {
                    var spec = ReadSpectrum(_reader, _fileFormatId, includePeaks, specIndex);
                    return spec;
                }
                return null;
            }
        }

        // Must reflect all changes to WriteSpectrum
        private static Spectrum ReadSpectrum(BinaryReader reader, int fileFormatId, bool includePeaks = true, int specIndex = -1)
        {
            var scanNum = reader.ReadInt32();
            string nativeId = String.Empty;
            if (fileFormatId > 150604)
            {
                var c = new char[NativeIdLength];
                reader.Read(c, 0, NativeIdLength);
                nativeId = (new string(c)).Trim();
            }
            var msLevel = reader.ReadByte();
            var elutionTime = reader.ReadDouble();
            double tic = -1;
            if (fileFormatId > 150604)
            {
                tic = reader.ReadSingle();
            }

            Spectrum spec;

            double calcTic;
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
                var peakList = ReadPeakList<Peak>(reader, fileFormatId, out calcTic, includePeaks);
                if (tic < 0)
                {
                    tic = calcTic;
                }
                spec =  new ProductSpectrum(peakList, scanNum)
                {
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
                // specIndex should only be set by GetMs1Spectrum, and should otherwise be negative
                if (specIndex >= 0)
                {
                    var peakList = ReadPeakList<Ms1Peak>(reader, fileFormatId, out calcTic, includePeaks, (ushort)specIndex);
                    if (tic < 0)
                    {
                        tic = calcTic;
                    }
                    spec = new Ms1Spectrum(scanNum, specIndex, peakList.ToArray());
                }
                else
                {
                    var peakList = ReadPeakList<Peak>(reader, fileFormatId, out calcTic, includePeaks);
                    if (tic < 0)
                    {
                        tic = calcTic;
                    }
                    spec = new Spectrum(peakList, scanNum);
                }
            }
            spec.MsLevel = msLevel;
            spec.ElutionTime = elutionTime;
            spec.NativeId = nativeId;
            spec.TotalIonCurrent = tic;
            return spec;
        }

        // All changes made here must be duplicated to ReadSpectrum() and GetPeakMetadataForSpectrum()
        public static void WriteSpectrum(Spectrum spec, BinaryWriter writer)
        {
            // scan number: 4
            writer.Write(spec.ScanNum);

            // NativeID: 50
            // pad or truncate to keep in limit (may have to change in future...)
            writer.Write(spec.NativeId.PadRight(NativeIdLength).ToCharArray(0, NativeIdLength), 0, NativeIdLength);

            // ms level: 1
            writer.Write(Convert.ToByte(spec.MsLevel));

            // elution time: 4
            writer.Write(spec.ElutionTime);

            // Total Ion Current: 4
            writer.Write(Convert.ToSingle(spec.TotalIonCurrent));

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
            //var peaks = spec.Peaks.ToList();
            Array.Sort(spec.Peaks);
            //peaks.Sort();
            // Number of peaks: 4
            writer.Write(spec.Peaks.Length);
            //writer.Write(peaks.Count);
            foreach (var peak in spec.Peaks)
            //foreach (var peak in peaks)
            {
                // m/z: 8
                writer.Write(peak.Mz);
                // intensity: 4
                writer.Write(Convert.ToSingle(peak.Intensity));
            }
        }

        // Must reflect all changes to WriteSpectrum
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
                // skip nativeId
                _reader.BaseStream.Seek(NativeIdLength, SeekOrigin.Current);
                var msLevel = _reader.ReadByte();
                var elutionTime = _reader.ReadDouble();
                var tic = _reader.ReadSingle();

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

        private static List<T> ReadPeakList<T>(BinaryReader reader, int fileFormatId, out double tic, bool includePeaks = true, ushort specIndex = 0)
            where T: Peak, new() // only allow Peak and derived classes, but requires public parameterless constructor
        {
            var peakList = new List<T>();
            var numPeaks = reader.ReadInt32();
            // Only used if fileFormatId < 150605
            tic = 0;

            // Skip the read if peaks aren't requested
            if (!includePeaks)
            {
                // first the number of peaks, then 12 bytes per peak (mz, double, 8 bytes, then intensity, single, 4 bytes)
                if (fileFormatId < 150605)
                {
                    // Read for calculating the tic
                    for (var i = 0; i < numPeaks; i++)
                    {
                        // skip 8 bytes from the mz
                        reader.ReadDouble();
                        // add up the intensities
                        tic += reader.ReadSingle();
                    }
                }
                else
                {
                    reader.BaseStream.Seek(numPeaks * 12, SeekOrigin.Current);
                }
                return peakList;
            }

            for (var i = 0; i < numPeaks; i++)
            {
                var mz = reader.ReadDouble();
                var intensity = reader.ReadSingle();
                // Only used if fileFormatId < 150605
                tic += intensity;
                var peak = new T();
                var ms1Peak = peak as Ms1Peak;
                if (ms1Peak != null)
                {
                    ms1Peak.SetMzIntensityIndices(mz, intensity, i, specIndex);
                }
                else
                {
                    peak.SetMzAndIntensity(mz, intensity);
                }
                peakList.Add(peak);
            }
            return peakList;
        }

        #endregion

        #region WriteAsPbf (obsolete)

        [Obsolete("Use PbfLcMsRun(string, IMassSpecDataReader, ...) for an optimized pbf creation process")]
        public static void WriteAsPbf(InMemoryLcMsRun imlr, string outputFilePath, IProgress<ProgressData> progress = null)
        {
            using (var writer = new BinaryWriter(File.Open(outputFilePath, FileMode.Create)))
            {
                WriteAsPbf(imlr, writer, progress);
            }
        }

        [Obsolete("Use PbfLcMsRun(string, IMassSpecDataReader, ...) for an optimized pbf creation process")]
        public static void WriteAsPbf(InMemoryLcMsRun imlr, BinaryWriter writer, IProgress<ProgressData> progress = null)
        {
            long countTotal = 1;
            long counter = 0;
            var progressData = new ProgressData(progress);
            progressData.Status = "Writing spectra data";

            var scanNumToSpecOffset = new long[imlr.NumSpectra + 1];
            var scanNumToIsolationWindow = new IsolationWindow[imlr.NumSpectra + 1];

            // Spectra
            countTotal = imlr.NumSpectra;
            counter = 0;
            progressData.StepRange(42.9); // SpecData: Approximately 43% of total file size
            long countMS2Spec = 0;
            for (var scanNum = imlr.MinLcScan; scanNum <= imlr.MaxLcScan; scanNum++)
            {
                progressData.Report(counter, countTotal);
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

            // Precursor ion chromatogram (MS1 spectra)
            var offsetBeginPrecursorChromatogram = writer.BaseStream.Position;

            var minMzIndex = imlr.Ms1PeakList.Any() ? PbfLcMsRun.GetMzBinIndex(imlr.Ms1PeakList[0].Mz) : 0;
            var maxMzIndex = imlr.Ms1PeakList.Any() ? PbfLcMsRun.GetMzBinIndex(imlr.Ms1PeakList[imlr.Ms1PeakList.Count - 1].Mz) : -1;

            var chromMzIndexToOffset = new long[maxMzIndex - minMzIndex + 1];
            var prevMzIndex = -1;
            counter = 0;
            countTotal = imlr.Ms1PeakList.Count;

            // All MS1 data sorted by mass, then scan number
            // In rare instances, imlr.Ms1PeakList will be blank (no MS1 spectra); that's OK
            progressData.Status = "Writing precursor chromatogram";
            progressData.StepRange(42.9 + 15.7); // Approximately 16% of total file size
            foreach (var peak in imlr.Ms1PeakList)
            {
                progressData.Report(counter, countTotal);
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

            // Product ion chromatogram (MSn spectra)
            var ms2PeakList = new List<LcMsPeak>();
            counter = 0;
            countTotal = countMS2Spec;
            progressData.Status = "Processing product ion chromatogram";
            progressData.StepRange(42.9 + 15.7 + (41.2 / 2)); // Approximately 41% of total file size
            foreach (var ms2ScanNum in imlr.GetScanNumbers(2))
            {
                progressData.Report(counter, countTotal);
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
            progressData.Status = "Writing product ion chromatogram";
            progressData.StepRange(42.9 + 15.7 + 41.2); // Approximately 41% of total file size
            foreach (var peak in ms2PeakList)
            {
                progressData.Report(counter, countTotal);
                counter++;
                writer.Write(peak.Mz);
                writer.Write((float)peak.Intensity);
                writer.Write(peak.ScanNum);
            }

            // Meta information
            var offsetBeginMetaInformation = writer.BaseStream.Position;
            progressData.Status = "Writing metadata";
            progressData.IsPartialRange = false;
            progressData.Report(99.8); // Metadata: Approximately 0.2% of total file size

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
            progressData.Report(99.9); // Metadata: Approximately 0.2% of total file size

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
            progressData.Report(100.0);
            writer.Write(FileFormatId); // 4
        }

        #endregion

        #region Pbf Creation

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
            var progressData = new ProgressData(progress);
            progressData.Status = "Writing spectra data";

            var scanNumToIsolationWindow = new Dictionary<int, IsolationWindow>(msdr.NumSpectra + 1);
            var ms1Scans = new List<int>();
            var ms2Scans = new List<int>();

            // Spectra
            long ms1PeakCount = 0;
            long ms2PeakCount = 0;
            var maxMs1Mz = double.MinValue;
            var minMs1Mz = double.MaxValue;
            var maxMs2Mz = double.MinValue;
            var minMs2Mz = double.MaxValue;
            var scanMetadata = new List<ScanMetadata>(msdr.NumSpectra);
            countTotal = msdr.NumSpectra;
            counter = 0;
            progressData.StepRange(42.9); // SpecData: Approximately 43% of total file size
            foreach (var spec in msdr.ReadAllSpectra())
            {
                progressData.Report(counter, countTotal);
                counter++;
                // Store offset, and write spectrum now
                _scanNumToSpecOffset.Add(spec.ScanNum, writer.BaseStream.Position);
                PbfLcMsRun.WriteSpectrum(spec, writer);

                // Handle other metadata stuff.
                var maxMz = double.MinValue;
                var minMz = double.MaxValue;
                if (spec.Peaks.Length > 0)
                {
                    minMz = spec.Peaks[0].Mz;
                    maxMz = spec.Peaks[spec.Peaks.Length - 1].Mz;
                }
                var productSpec = spec as ProductSpectrum;
                scanNumToIsolationWindow[spec.ScanNum] = null;
                if (productSpec != null)
                {
                    scanNumToIsolationWindow[spec.ScanNum] = productSpec.IsolationWindow;
                    ms2Scans.Add(productSpec.ScanNum);
                    //ms2PeakList.AddRange(productSpec.Peaks.Select(peak => new LcMsPeak(peak.Mz, peak.Intensity, productSpec.ScanNum)));
                    ms2PeakCount += productSpec.Peaks.Length;
                    if (minMz < minMs2Mz)
                    {
                        minMs2Mz = minMz;
                    }
                    if (maxMs2Mz < maxMz)
                    {
                        maxMs2Mz = maxMz;
                    }
                }
                else
                {
                    ms1Scans.Add(spec.ScanNum);
                    //ms1PeakList.AddRange(spec.Peaks.Select(peak => new LcMsPeak(peak.Mz, peak.Intensity, spec.ScanNum)));
                    ms1PeakCount += spec.Peaks.Length;
                    if (minMz < minMs1Mz)
                    {
                        minMs1Mz = minMz;
                    }
                    if (maxMs1Mz < maxMz)
                    {
                        maxMs1Mz = maxMz;
                    }
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
            }

            ms1Scans.Sort();
            ms2Scans.Sort();

            // Precursor ion chromatogram (MS1 spectra)
            _offsetPrecursorChromatogramBegin = writer.BaseStream.Position;

            if (ms1PeakCount > 0)
            {
                progressData.Status = "Writing precursor chromatogram";
                if (ms2PeakCount > 0)
                {
                    progressData.StepRange(42.9 + 15.7); // Approximately 16% of total file size, on standard LCMS file
                }
                else
                {
                    progressData.StepRange(42.9 + 15.7 + 41.2); // Use MS2 reserved chunk also, no MS2 spectra
                }
                var prog = new Progress<ProgressData>(p =>
                {
                    progressData.StatusInternal = p.Status;
                    progressData.Report(p.Percent);
                });
                CreateAndOutputMsXChromatogram_Merge(writer, ms1Scans, ms1PeakCount, true, minMs1Mz, maxMs1Mz, prog);
            }
            else
            {
                // Initialize this to an empty, zero-length array to prevent null reference exceptions
                _chromMzIndexToOffset = new long[0];
            }

            // Product ion chromatogram (MSn spectra)
            _offsetProductChromatogramBegin = writer.BaseStream.Position;
            _offsetPrecursorChromatogramEnd = _offsetProductChromatogramBegin;

            if (ms2PeakCount > 0)
            {
                progressData.Status = "Writing product chromatogram";
                progressData.StepRange(42.9 + 15.7 + 41.2); // Approximately 41% of total file size, on standard LCMS file
                var prog = new Progress<ProgressData>(p =>
                {
                    progressData.StatusInternal = p.Status;
                    progressData.Report(p.Percent);
                });
                CreateAndOutputMsXChromatogram_Merge(writer, ms2Scans, ms2PeakCount, false, minMs2Mz, maxMs2Mz, prog);
            }

            // Meta information
            _offsetProductChromatogramEnd = writer.BaseStream.Position;
            var offsetBeginMetaInformation = _offsetProductChromatogramEnd;
            progressData.Status = "Writing metadata";
            progressData.IsPartialRange = false;
            progressData.Report(99.8); // Metadata: Approximately 0.2% of total file size

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
            progressData.Report(99.9); // Metadata: Approximately 0.2% of total file size

            var prevOffset = offsetBeginMetaInformation;
            if (ms1PeakCount > 0)
            {
                for (var i = _chromMzIndexToOffset.Length - 2; i >= 0; i--)
                {
                    if (_chromMzIndexToOffset[i] < _offsetPrecursorChromatogramBegin)
                        _chromMzIndexToOffset[i] = prevOffset;
                    else
                        prevOffset = _chromMzIndexToOffset[i];
                }

                foreach (var offset in _chromMzIndexToOffset.Take(_chromMzIndexToOffset.Length - 1))
                {
                    writer.Write(offset);
                }

                _chromMzIndexToOffset[_chromMzIndexToOffset.Length - 1] = _offsetPrecursorChromatogramEnd;
            }

            writer.Write(_offsetPrecursorChromatogramBegin); // 8
            writer.Write(_offsetProductChromatogramBegin); // 8
            writer.Write(offsetBeginMetaInformation); // 8
            progressData.Report(100.0);
            writer.Write(FileFormatId); // 4
        }

        private int CalculateMaxSpecsInMem(double memFreeKB, int scansCount)
        {
            //Approximate Maximum amount of memory used for a set number of peaks:
            //     3.5 million: 500 MB
            //     10  million: 1.1 GB
            //     79  million: 6.6 GB
            const int max = 25000000;
            int maxInMemoryPerSpec = (int)(memFreeKB * 1024 / 2 / scansCount / (20 + 8));
            if (maxInMemoryPerSpec * scansCount > max) // Set a hard limit at 10 millions peaks in memory at once (only exception is the minimum)
                maxInMemoryPerSpec = max / scansCount;
            if (maxInMemoryPerSpec < 5)
                maxInMemoryPerSpec = 5;
            return maxInMemoryPerSpec;
        }

        private void CreateAndOutputMsXChromatogram_Merge(BinaryWriter writer, List<int> scansForMsLevelX, long totalPeaksCount, bool isMs1List, double minMz, double maxMz, IProgress<ProgressData> progress)
        {
            // Other thought for a slower, but really low-memory chromatogram creator:
            //   Make sure peaks are sorted ascending when written, and then jump through all spectra, performing a massive merge sort on the peaks and outputting the lowest spectra
            var progData = new ProgressData(progress);
            var count = 0;
            var countTotal = scansForMsLevelX.Count;
            var prevMzIndex = -1;

            // List overhead per item: is array-backed, so item (+ pointer to it)
            // item size: double, double, int, so 20 bytes
            //var avgPeaksPerSpec = (double)totalPeaksCount / scansForMsLevelX.Count;
            // Testing: Approx. 10 peaks per KB - it underestimates on the non-low-memory implementation, but only by 2-3 peaks/KB
            //var totalKBytesInMem = (totalPeaksCount * (20 + 8)) / 1024;
            var totalKBytesInMem = totalPeaksCount / 10;
            double memoryFreeKB = 0;
            double totalMemKB = 0;
            //var osVersion = Environment.OSVersion.Version;
            //if (osVersion < new Version(6, 0)) // Older than Vista
            //{
            //    // For pre-Vista: "SELECT * FROM Win32_LogicalMemoryConfiguration", and a different property.
            //    // Have no good systems to test on (Sequest head nodes??)
            //}
            foreach (var item in new System.Management.ManagementObjectSearcher("SELECT * FROM CIM_OperatingSystem").Get())
            {
                memoryFreeKB += double.Parse(item["FreePhysicalMemory"].ToString());
                //foreach (var p in item.Properties)
                //{
                //    Console.WriteLine("{0}: {1}", p.Name, p.Value);
                //}
            }
            // Get total physical memory
            foreach (var item in new System.Management.ManagementObjectSearcher("SELECT * FROM Win32_ComputerSystem").Get())
            {
                // TotalPhysicalMemory is in Bytes, so divide by 1024
                totalMemKB += double.Parse(item["TotalPhysicalMemory"].ToString()) / 1024.0;
                //foreach (var p in item.Properties)
                //{
                //    Console.WriteLine("{0}: {1}", p.Name, p.Value);
                //}
            }

            //Console.WriteLine("memoryFreeKB: " + memoryFreeKB);
            //Console.WriteLine("totalKBytesInMem: " + totalKBytesInMem);
            //Console.WriteLine("NumScans: " + scansForMsLevelX.Count);
            //Console.WriteLine("NumPeaks: " + totalPeaksCount);

            // Configure a reserve amount to avoid using all physical memory, which has the cost of excessive paging
            var quarterTotalPhysicalMemory = totalMemKB / 4;
            // Cut the reserve down on systems with large amounts of physical memory (16GB and greater)
            while (quarterTotalPhysicalMemory > 4194304)
            {
                quarterTotalPhysicalMemory = quarterTotalPhysicalMemory / 2;
            }
            var memFreeLessReserve = memoryFreeKB - quarterTotalPhysicalMemory;

            // Use a custom split-list implementation: Reduce the number of items to sort each time
            // The enumerator access is O(1)
            // 8 byte (pointer) overhead per item, plus ~30 bytes per list (number of lists is numSpectra / 1,000,000 + 2)
            // item size: double, double, int, so 20 bytes
            // Perform a massive merge-sort style read/write, to minimize memory usage - set a minimum of 5
            int maxInMemoryPerSpec = CalculateMaxSpecsInMem(memFreeLessReserve, scansForMsLevelX.Count);

            //Console.WriteLine("maxInMemoryPerSpec: " + maxInMemoryPerSpec);

            progData.Status = "Writing product chromatogram";
            if (isMs1List)
            {
                progData.Status = "Writing precursor chromatogram";
            }
            lock (_filelock)
            {
                var peaks = new SplitLcMsPeakLists(minMz, maxMz, totalPeaksCount);
                var peaksCount = 0;
                var metadata = new Dictionary<int, ScanPeakMetaData>(scansForMsLevelX.Count);

                if (_reader == null)
                {
                    _reader =
                        new BinaryReader(File.Open(PbfFilePath, FileMode.Open, FileAccess.Read,
                                FileShare.ReadWrite));
                }

                progData.StepRange(5);
                // Read in metadata, and the first "maxInMemoryPerSpec" peaks of each scan (min 5, max 25000000 / numSpectra)
                foreach (var scan in scansForMsLevelX)
                {
                    progData.Report(count, countTotal);
                    count++;
                    var datum = GetPeakMetaDataForSpectrum(scan);
                    peaksCount += datum.NumPeaks;
                    metadata.Add(datum.ScanNum, datum);
                    peaks.AddRange(datum.ReadPeaks(_reader, maxInMemoryPerSpec));
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
                var peaksToRemove = new List<LcMsPeak>();
                count = 0;
                prevMzIndex = -1;
                if (isMs1List)
                {
                    _minMs1Mz = minMz;
                    _minMzIndex = GetMzBinIndex(_minMs1Mz);
                    _maxMs1Mz = maxMz;
                    _maxMzIndex = GetMzBinIndex(_maxMs1Mz);
                }
                var chromMzIndexToOffset = new Dictionary<int, long>(_maxMzIndex - _minMzIndex);

                var peaksEnum = peaks.GetReusableEnumerator();
                while (peaks.Count > 0)
                {
                    while (peaksEnum.MoveNext())
                    {
                        var peak = peaksEnum.Current;
                        progData.Report(count, peaksCount);
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
                        datum.ReadOne();

                        if (datum.Count < 1)
                        {
                            if (datum.MorePeaksToRead)
                            {
                                break;
                            }
                            else
                            {
                                metadata.Remove(datum.ScanNum);
                                var numDead = datum.NumPeaks < maxInMemoryPerSpec ? datum.NumPeaks : maxInMemoryPerSpec;
                                if (numDead >= peaksToRemove.Count - 500)
                                {
                                    peaksToRemove.Clear();
                                }
                                else
                                {
                                    peaksToRemove.RemoveRange(peaksToRemove.Count - numDead - 1, numDead);
                                }
                            }
                        }
                    }
                    maxInMemoryPerSpec = CalculateMaxSpecsInMem(memFreeLessReserve, metadata.Count);

                    // add new entries back into the list from the spectra that the peaks came from
                    foreach (var datum in metadata.Values)
                    {
                        if (datum.Count < maxInMemoryPerSpec && datum.MorePeaksToRead)
                        {
                            peaks.AddRange(datum.ReadPeaksReuse(_reader, maxInMemoryPerSpec - datum.Count, peaksToRemove));
                        }
                    }

                    peaksEnum = peaks.GetReusableEnumerator();
                }

                if (isMs1List)
                {
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

            GC.Collect(GC.MaxGeneration, GCCollectionMode.Forced, false); // We just killed a large number of objects when the list of peaks went out of scope.
        }

        private class SplitLcMsPeakLists
        {
            private readonly List<List<LcMsPeak>> _lists;
            private readonly List<bool> _sorted;
            private int _mod = 0;
            private double _mzMod = 0;
            private const int Divisor = 1000000;
            private double _minMz;
            private double _maxMz;

            public SplitLcMsPeakLists(double minMz, double maxMz, long peakCount)
            {
                _minMz = minMz;
                _maxMz = maxMz;
                _mod = (int)(peakCount / Divisor + 1);
                _mzMod = (_maxMz - _minMz) / (_mod - 1);
                _lists = new List<List<LcMsPeak>>(_mod + 2);
                _sorted = new List<bool>(_mod + 2);
                for (int i = 0; i <= _mod + 1; i++)
                {
                    _lists.Add(new List<LcMsPeak>(1000));
                    _sorted.Add(false);
                }
            }

            public void Add(LcMsPeak peak)
            {
                if (peak == null)
                {
                    return;
                }
                RemoveEnumerated();
                var list = (int)((peak.Mz - _minMz) / _mzMod);
                if (list < 0)
                {
                    list = 0;
                }
                if (list > _mod + 1)
                {
                    list = _mod + 1;
                }
                _lists[list].Add(peak);
                _sorted[list] = false;
            }

            public void AddRange(IEnumerable<LcMsPeak> peaks)
            {
                foreach (var peak in peaks)
                {
                    Add(peak);
                }
            }

            private void RemoveEnumerated()
            {
                if (_enumerated > 0)
                {
                    // remove the already-enumerated peaks
                    for (int i = 0; i < _lists.Count && _enumerated > 0; i++)
                    {
                        if (_lists[i].Count <= 0)
                        {
                            continue;
                        }
                        if (_lists[i].Count <= _enumerated)
                        {
                            _enumerated -= _lists[i].Count;
                            _lists[i].Clear();
                        }
                        else
                        {
                            _lists[i].RemoveRange(0, (int)_enumerated);
                            _enumerated = 0;
                        }
                    }
                }
            }

            private void SortList(int index)
            {
                if (!_sorted[index])
                {
                    _lists[index].Sort();
                    _sorted[index] = true;
                }
            }

            public int Count
            {
                get
                {
                    return _lists.Sum(list => list.Count);
                }
            }

            private long _enumerated = 0;

            public IEnumerator<LcMsPeak> GetReusableEnumerator()
            {
                for (var i = 0; i < _lists.Count; i++)
                {
                    if (_lists[i].Count > 0)
                    {
                        SortList(i);
                        foreach (var peak in _lists[i])
                        {
                            _enumerated++;
                            yield return peak;
                        }
                        _enumerated -= _lists[i].Count; // prevent creep
                        _lists[i].Clear();
                        _lists[i] = new List<LcMsPeak>(10);
                    }
                }
            }
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
            public int PeaksRead { get { return _numPeaksRead; } }
            private int _numPeaksRead;
            private long _nextPeakOffset;
            public int Count { get; private set; }

            public bool MorePeaksToRead
            {
                get { return _numPeaksRead < NumPeaks; }
            }

            public ScanPeakMetaData(int scanNum, BinaryReader reader)
            {
                Count = 0;
                ScanNum = scanNum;
                NumPeaks = reader.ReadInt32();
                _numPeaksRead = 0;
                _nextPeakOffset = reader.BaseStream.Position;
            }

            public IEnumerable<LcMsPeak> ReadPeaks(BinaryReader reader, int numPeaksToRead)
            {
                if (_numPeaksRead >= NumPeaks)
                {
                    yield return null;
                }
                reader.BaseStream.Seek(_nextPeakOffset, SeekOrigin.Begin);

                for (int i = 0; i < numPeaksToRead && _numPeaksRead < NumPeaks; i++)
                {
                    var mz = reader.ReadDouble();
                    var intensity = reader.ReadSingle();
                    Count++;
                    _numPeaksRead++;
                    yield return new LcMsPeak(mz, intensity, ScanNum);
                }
                _nextPeakOffset = reader.BaseStream.Position;
            }

            public IEnumerable<LcMsPeak> ReadPeaksReuse(BinaryReader reader, int numPeaksToRead, List<LcMsPeak> peakStock)
            {
                if (_numPeaksRead >= NumPeaks)
                {
                    yield return null;
                }
                reader.BaseStream.Seek(_nextPeakOffset, SeekOrigin.Begin);

                for (int i = 0; i < numPeaksToRead && _numPeaksRead < NumPeaks; i++)
                {
                    var mz = reader.ReadDouble();
                    var intensity = reader.ReadSingle();
                    Count++;
                    _numPeaksRead++;
                    var lastIndex = peakStock.Count - 1;
                    if (lastIndex >= 0)
                    {
                        var peak = peakStock[lastIndex];
                        peakStock.RemoveAt(lastIndex);
                        yield return peak.ReplaceData(mz, intensity, ScanNum);
                    }
                    else
                    {
                        yield return new LcMsPeak(mz, intensity, ScanNum);
                    }
                }
                _nextPeakOffset = reader.BaseStream.Position;
            }

            public void ReadOne()
            {
                Count--;
            }
        }

        #endregion

        #region XIC reading functions

        public static int GetMzBinIndex(double mz)
        {
            return (int)(mz / MzBinSize);
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

            return Xic.GetSelectedXic(xic);
        }

        // beginOffset: inclusive, endOffset: exclusive
        private Xic GetChromatogramRange(double minMz, double maxMz, long beginOffset, long endOffset, long targetOffset)
        {
            var xic = GetXicPointsWithin(minMz, maxMz, beginOffset, endOffset, targetOffset);
            if (!xic.Any()) return xic;

            return xic;
        }

        // Must reflect any changes made in the chromatogram creation
        private Xic GetXicPointsWithin(double minMz, double maxMz, long beginOffset, long endOffset,
            long targetOffset)
        {
            var xic = new Xic();
            var curOffset = targetOffset - NumBytePeak;
            var cacheHigher = false;
            var cacheLower = false;
            if (endOffset <= _offsetPrecursorChromatogramEnd)
            {
                cacheHigher = HigherPrecursorChromatogramCacheSize >= 20;
                cacheLower = LowerPrecursorChromatogramCacheSize >= 20;
                _precursorChromatogramCache.Clear();
                if (cacheHigher)
                {
                    if (endOffset < targetOffset + HigherPrecursorChromatogramCacheSize * NumBytePeak)
                    {
                        endOffset = targetOffset + HigherPrecursorChromatogramCacheSize * NumBytePeak;
                    }
                    if (endOffset > _offsetPrecursorChromatogramEnd)
                    {
                        endOffset = _offsetPrecursorChromatogramEnd;
                    }
                }
                if (cacheLower)
                {
                    if (beginOffset > targetOffset - LowerPrecursorChromatogramCacheSize * NumBytePeak)
                    {
                        beginOffset = targetOffset - LowerPrecursorChromatogramCacheSize * NumBytePeak;
                    }
                    if (beginOffset < _offsetPrecursorChromatogramBegin)
                    {
                        beginOffset = _offsetPrecursorChromatogramBegin;
                    }
                }
            }
            var doCache = cacheLower || cacheHigher;

            lock (_filelock)
            {
                var cacheCount = 0;
                // go down
                while (curOffset >= beginOffset)
                {
                    _reader.BaseStream.Seek(curOffset, SeekOrigin.Begin);
                    var mz = _reader.ReadDouble();
                    var intensity = _reader.ReadSingle();
                    var scanNum = _reader.ReadInt32();
                    if (mz < minMz)
                    {
                        if (!cacheLower || cacheCount >= LowerPrecursorChromatogramCacheSize)
                        {
                            break;
                        }
                        _precursorChromatogramCache.Add(new XicPoint(scanNum, mz, intensity));
                        cacheCount++;
                    }
                    else
                    {
                        xic.Add(new XicPoint(scanNum, mz, intensity));
                        if (doCache)
                        {
                            _precursorChromatogramCache.Add(new XicPoint(scanNum, mz, intensity));
                        }
                    }
                    curOffset -= NumBytePeak;
                }

                if (doCache && curOffset < _offsetPrecursorChromatogramBegin)
                {
                    _precursorChromatogramCache.Add(new XicPoint(int.MinValue, double.NegativeInfinity, 0));
                }

                cacheCount = 0;
                // go up
                curOffset = targetOffset;
                while (curOffset < endOffset)
                {
                    _reader.BaseStream.Seek(curOffset, SeekOrigin.Begin);
                    var mz = _reader.ReadDouble();
                    var intensity = _reader.ReadSingle();
                    var scanNum = _reader.ReadInt32();
                    if (mz > maxMz)
                    {
                        if (!cacheHigher || cacheCount >= HigherPrecursorChromatogramCacheSize)
                        {
                            break;
                        }
                        _precursorChromatogramCache.Add(new XicPoint(scanNum, mz, intensity));
                        cacheCount++;
                    }
                    else
                    {
                        xic.Add(new XicPoint(scanNum, mz, intensity));
                        if (doCache)
                        {
                            _precursorChromatogramCache.Add(new XicPoint(scanNum, mz, intensity));
                        }
                    }
                    curOffset += NumBytePeak;
                }

                if (doCache && curOffset >= _offsetPrecursorChromatogramEnd)
                {
                    _precursorChromatogramCache.Add(new XicPoint(int.MaxValue, double.PositiveInfinity, 0));
                }
            }

            _precursorChromatogramCache.Sort((x,y) => x.Mz.CompareTo(y.Mz));

            return xic;
        }

        #endregion
    }
}
