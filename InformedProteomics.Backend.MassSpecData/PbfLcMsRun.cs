using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using InformedProteomics.Backend.Data.Spectrometry;
using PRISM;
using PSI_Interface.CV;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// A file-backed object for accessing mass spec data. Data is read from the vendor format
    /// to a binary PBF file, and that file is used for fast access to spectra and extracted ion chromatograms
    /// </summary>
    public class PbfLcMsRun : LcMsRun, IPbfLcMsRun
    {
        // Ignore Spelling: Baf, LcMs, Nums, overridable, Pbf, struct, structs

        /// <summary>
        /// File extension
        /// </summary>
        public const string FileExtensionConst = ".pbf";

        /// <summary>
        /// File extension used for this type
        /// </summary>
        protected virtual string FileExtensionVirtual => FileExtensionConst;

        /// <summary>
        /// File extension - overridable. Returns <see cref="FileExtensionConst"/> for current type. See <see cref="FileExtensionConst"/> for static access.
        /// </summary>
        public string FileExtension => FileExtensionVirtual;

        /// <summary>
        /// True if the file contains precursor or product chromatograms
        /// </summary>
        public virtual bool ContainsChromatograms => true;

        /// <summary>
        /// The current FileFormatId, which is written to the file, and checked before a file is read
        /// </summary>
        /// <remarks>This constant should be incremented by 1 if the binary file format is changed</remarks>
        public const int FileFormatId = 150608;
        private const int EarliestSupportedFileFormatId = 150604;

        /* File format id history
         * 150604: Earliest supported, has all data needed by InformedProteomics projects
         * 150605: Added Total Ion Current and Native ID fields to spectrum output.
         * 150606: Added ScanNum to the output in the metadata section to support skipped scans
         * 150607: Added the name/path of the source file, and the format of the native id
         * 150608: Added a file checksum (SHA-1)
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
            return MassSpecDataReaderFactory.ChangeExtension(specFileName, FileExtensionConst);
        }

        /// <summary>
        /// Convert a spec file to pbf, and return an LcMsRun that uses the pbf file
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
        /// Convert a spec file to pbf, and return an LcMsRun that uses the pbf file
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
        /// Convert a spec file to pbf, and return an LcMsRun that uses the pbf file
        /// </summary>
        /// <param name="specFilePath"></param>
        /// <param name="specReader">Data reader; if not a PbfLcMsRun, it will be closed when pbf file creation is finished</param>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        /// <param name="progress"></param>
        /// <returns></returns>
        /// <remarks>It is recommended that "MassSpecDataReaderFactory.NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path</remarks>
        public static LcMsRun GetLcMsRun(string specFilePath, IMassSpecDataReader specReader, double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold,
            IProgress<ProgressData> progress = null)
        {
            if (specReader is PbfLcMsRun lcMsRun)
            {
                return lcMsRun;
            }

            return new PbfLcMsRun(specFilePath, specReader, string.Empty, precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold, progress);
        }

        /// <summary>
        /// Convert a spec file to pbf
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
        public static string ConvertToPbf(string specFilePath, IMassSpecDataReader specReader, string pbfFilePath = null,
            IProgress<ProgressData> progress = null)
        {
            if (specFilePath.ToLower().EndsWith(FileExtensionConst, StringComparison.OrdinalIgnoreCase))
            {
                return specFilePath;
            }

            var pbfPath = pbfFilePath;
            var fileName = string.Empty;
            var tempPath = string.Empty;

            var prog = new Progress<ProgressData>();
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
                pbfPath = GetPbfFileName(MassSpecDataReaderFactory.NormalizeDatasetPath(specFilePath));
                fileName = Path.GetFileName(pbfPath);
                if (String.IsNullOrEmpty(fileName))
                {
                    throw new ArgumentException("Cannot create .pbf cache file", nameof(specFilePath));
                }

                tempPath = Path.Combine(Path.GetTempPath(), fileName);
                // Return the temp path if the pbf file of proper format already exists in the temp directory
                if (File.Exists(tempPath) && CheckFileFormatVersion(tempPath, out isCurrent) && isCurrent)
                {
                    return tempPath;
                }
            }

            if (!File.Exists(pbfPath) || !(CheckFileFormatVersion(pbfPath, out isCurrent) && isCurrent))
            {
                if (specReader == null)
                {
                    throw new Exception("Unsupported file format!");
                }
                var run = new InMemoryLcMsRun(specReader, 0, 0, prog);
                try
                {
                    progData.StepRange(100.0);
                    WriteAsPbf(run, pbfPath, prog);
                }
                catch (UnauthorizedAccessException) // Cannot write to same directory, attempt to write to temp directory
                {
                    // Fail out if the output path was specified, and we cannot write to it.
                    if (!string.IsNullOrWhiteSpace(pbfFilePath))
                    {
                        throw;
                    }
                    //var fileName = Path.GetFileName(pbfFilePath);
                    if (string.IsNullOrEmpty(fileName))
                    {
                        throw; // invalid path?
                    }
                    //var tempPath = Path.Combine(Path.GetTempPath(), fileName);
                    if (!File.Exists(tempPath) || !(CheckFileFormatVersion(tempPath, out isCurrent) && isCurrent))
                    {
                        WriteAsPbf(run, tempPath, prog);
                    }

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
        public virtual string GetCheckPbfFilePath(string specFilePath, out string pbfPath, out string fileName, out string tempPath)
        {
            return GetCheckPbfFilePath(specFilePath, out pbfPath, out fileName, out tempPath, FileExtensionConst);
        }

        /// <summary>
        /// Gets valid possible pbf file paths
        /// </summary>
        /// <param name="specFilePath">Path to the spectra file</param>
        /// <param name="pbfPath">Path to the default pbf file (in the same folder as the spectra file dataset)</param>
        /// <param name="fileName"></param>
        /// <param name="tempPath"></param>
        /// <param name="extension">The extension expected for the pbf file</param>
        /// <returns>The default path to the pbf file, unless a valid pbf file exists at the temp path</returns>
        protected internal string GetCheckPbfFilePath(string specFilePath, out string pbfPath, out string fileName, out string tempPath, string extension)
        {
            // Calls "NormalizeDatasetPath" to make sure we save the file to the containing directory
            pbfPath = MassSpecDataReaderFactory.ChangeExtension(MassSpecDataReaderFactory.NormalizeDatasetPath(specFilePath), extension);
            fileName = Path.GetFileName(pbfPath);
            if (String.IsNullOrEmpty(fileName))
            {
                throw new ArgumentException("Cannot create .pbf cache file", nameof(specFilePath));
            }

            tempPath = Path.Combine(Path.GetTempPath(), fileName);
            if (File.Exists(pbfPath) && CheckFileFormatVersion(pbfPath, out _))
            {
                return pbfPath;
            }

            // Return the temp path if the pbf file of proper format already exists in the temp directory
            if (File.Exists(tempPath) && CheckFileFormatVersion(tempPath, out _))
            {
                return tempPath;
            }
            return pbfPath;
        }

        #endregion

        #region Constructors

        /// <summary>
        /// Constructor for opening a PBF file
        /// </summary>
        /// <param name="specFileName"></param>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        public PbfLcMsRun(string specFileName, double precursorSignalToNoiseRatioThreshold = 0.0, double productSignalToNoiseRatioThreshold = 0.0)
        {
            _precursorSignalToNoiseRatioThreshold = precursorSignalToNoiseRatioThreshold;
            _productSignalToNoiseRatioThreshold = productSignalToNoiseRatioThreshold;

            OpenPbfFile(specFileName);
        }

        /// <summary>
        /// Constructor for creating and/or opening a PBF file
        /// </summary>
        /// <param name="specFileName"></param>
        /// <param name="msDataReader"></param>
        /// <param name="pbfFileName"></param>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        /// <param name="progress"></param>
        /// <param name="keepDataReaderOpen">use 'true' if the data reader should not be closed when finished creating the PBF file</param>
        /// <param name="scanStart">Minimum scan number to include in the .PBF file; 0 to disable this filter</param>
        /// <param name="scanEnd">Maximum scan number to include to the .PBF file; 0 to disable this filter</param>
        public PbfLcMsRun(string specFileName, IMassSpecDataReader msDataReader, string pbfFileName = null,
                          double precursorSignalToNoiseRatioThreshold = 0.0, double productSignalToNoiseRatioThreshold = 0.0,
                          IProgress<ProgressData> progress = null, bool keepDataReaderOpen = false, int scanStart = 0, int scanEnd = 0)
        {
            _precursorSignalToNoiseRatioThreshold = precursorSignalToNoiseRatioThreshold;
            _productSignalToNoiseRatioThreshold = productSignalToNoiseRatioThreshold;

            GetPbfFile(specFileName, msDataReader, pbfFileName, progress, keepDataReaderOpen, scanStart, scanEnd);
        }

        /// <summary>
        /// Constructor for creating and/or opening a PBF file
        /// </summary>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        protected internal PbfLcMsRun(double precursorSignalToNoiseRatioThreshold = 0.0, double productSignalToNoiseRatioThreshold = 0.0)
        {
            _precursorSignalToNoiseRatioThreshold = precursorSignalToNoiseRatioThreshold;
            _productSignalToNoiseRatioThreshold = productSignalToNoiseRatioThreshold;
        }

        [Obsolete("Only here for compatibility purposes with the obsolete WriteAsPbf function", true)]
        private PbfLcMsRun()
        {
        }

        #endregion

        #region Constructor support functions

        /// <summary>
        /// Given a spec file path and other information, either open an existing pbf corresponding to the spec file path, or create a new one
        /// </summary>
        /// <param name="specFileName"></param>
        /// <param name="msDataReader"></param>
        /// <param name="pbfFileName"></param>
        /// <param name="progress"></param>
        /// <param name="keepDataReaderOpen"></param>
        /// <param name="scanStart">Minimum scan number to include in the .PBF file; 0 to disable this filter</param>
        /// <param name="scanEnd">Maximum scan number to include to the .PBF file; 0 to disable this filter</param>
        protected internal void GetPbfFile(
            string specFileName,
            IMassSpecDataReader msDataReader,
            string pbfFileName,
            IProgress<ProgressData> progress,
            bool keepDataReaderOpen = false,
            int scanStart = 0,
            int scanEnd = 0)
        {
            var pbfPath = GetCheckPbfFilePath(specFileName, out _, out _, out var tempPath);
            if (!string.IsNullOrWhiteSpace(pbfFileName))
            {
                pbfPath = pbfFileName;
            }

            if (specFileName.EndsWith(FileExtension, StringComparison.OrdinalIgnoreCase) ||
                File.Exists(pbfPath) && CheckFileFormatVersion(pbfPath, out var isCurrent) && isCurrent)
            {
                // The existing pbf file is valid; use it

                if (!specFileName.EndsWith(FileExtension, StringComparison.OrdinalIgnoreCase))
                {
                    specFileName = pbfPath;
                }

                OpenPbfFile(specFileName);
                if (!keepDataReaderOpen)
                {
                    msDataReader?.Dispose();
                }
                return;
            }

            // Create a new pbf file
            BuildPbfFile(specFileName, msDataReader, pbfPath, tempPath, progress, keepDataReaderOpen, scanStart, scanEnd);
        }

        /// <summary>
        /// Code for writing a PBF file. Should only be called from a constructor.
        /// </summary>
        /// <param name="specFileName"></param>
        /// <param name="msDataReader"></param>
        /// <param name="pbfPath"></param>
        /// <param name="tempPath"></param>
        /// <param name="progress"></param>
        /// <param name="keepDataReaderOpen"></param>
        /// <param name="scanStart">Minimum scan number to include in the .PBF file; 0 to disable this filter</param>
        /// <param name="scanEnd">Maximum scan number to include to the .PBF file; 0 to disable this filter</param>
        protected internal void BuildPbfFile(
            string specFileName,
            IMassSpecDataReader msDataReader,
            string pbfPath,
            string tempPath,
            IProgress<ProgressData> progress,
            bool keepDataReaderOpen = false,
            int scanStart = 0,
            int scanEnd = 0)
        {
            if (msDataReader == null)
            {
                msDataReader = MassSpecDataReaderFactory.GetMassSpecDataReader(specFileName);
            }
            NumSpectra = msDataReader.NumSpectra;
            RawFilePath = specFileName;
            NativeIdFormat = msDataReader.NativeIdFormat;
            NativeFormat = msDataReader.NativeFormat;
            SrcFileChecksum = msDataReader.SrcFileChecksum.ToLower().Replace("-", "");

            try
            {
                var pbfFile = new FileInfo(pbfPath);
                if (pbfFile.Directory != null && !pbfFile.Directory.Exists)
                {
                    pbfFile.Directory.Create();
                }

                PbfFilePath = pbfFile.FullName;

                using (var writer =
                    new BinaryWriter(File.Open(pbfFile.FullName, FileMode.Create, FileAccess.ReadWrite, FileShare.Read)))
                {
                    WriteToPbf(msDataReader, writer, scanStart, scanEnd, progress);
                }
            }
            catch (UnauthorizedAccessException)
            {
                PbfFilePath = tempPath;
                using (var writer =
                    new BinaryWriter(File.Open(tempPath, FileMode.Create, FileAccess.ReadWrite, FileShare.Read)))
                {
                    WriteToPbf(msDataReader, writer, scanStart, scanEnd, progress);
                }
            }
            finally
            {
                if (!keepDataReaderOpen)
                {
                    msDataReader.Dispose();
                }
            }
            FileFormatVersion = FileFormatId.ToString();

            try
            {
                _reader = new BinaryReader(File.Open(PbfFilePath, FileMode.Open, FileAccess.Read, FileShare.Read));
            }
            catch
            {
                // For when something else, like a file scanner, has opened the file with read-write access
                _reader = new BinaryReader(File.Open(PbfFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));
            }

            CreatePrecursorNextScanMap();
        }

        /// <summary>
        /// Code for opening a PBF file. Should only be called from the constructors.
        /// </summary>
        /// <param name="specFileName"></param>
        protected internal void OpenPbfFile(string specFileName)
        {
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

            lock (_fileLock)
            {
                if (!ReadMetaInfo())
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

        /// <summary>
        /// Check the file format version of the specified file to see if it is readable with the current version.
        /// </summary>
        /// <param name="filePath">path to the file to check</param>
        /// <param name="isCurrent">true if the format is the same as the current format</param>
        /// <returns>True if the format is readable</returns>
        public static bool CheckFileFormatVersion(string filePath, out bool isCurrent)
        {
            isCurrent = false;
            var pbfFile = new FileInfo(filePath);
            if (!pbfFile.Exists || pbfFile.Length < sizeof(int))
            {
                return false;
            }

            var fs = new FileStream(pbfFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
            using (var reader = new BinaryReader(fs))
            {
                fs.Seek(-1 * sizeof(int), SeekOrigin.End);

                var fileFormatId = reader.ReadInt32();
                if (fileFormatId > FileFormatId && fileFormatId < EarliestSupportedFileFormatId)
                {
                    return false;
                }

                if (fileFormatId == FileFormatId)
                {
                    isCurrent = true;
                }
            }
            return true;
        }

        #endregion

        #region Member Variables

        /// <summary>
        /// The path to this pbf file
        /// </summary>
        public string PbfFilePath { get; private set; }

        /// <summary>
        /// Raw file path (as of when/where the file was created)
        /// </summary>
        public string RawFilePath { get; private set; }

        /// <summary>
        /// SHA-1 Checksum of the pbf file, calculated on first access to this property - lowercase, hex only
        /// </summary>
        public string PbfFileChecksum
        {
            get
            {
                if (string.IsNullOrWhiteSpace(_pbfFileChecksum))
                {
                    using (var fs = new FileStream(PbfFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
                    using (var sha1 = new SHA1Managed())
                    {
                        var hash = sha1.ComputeHash(fs);
                        _pbfFileChecksum = BitConverter.ToString(hash).ToLower().Replace("-", "");
                    }
                }
                return _pbfFileChecksum;
            }
        }

        private string _pbfFileChecksum = string.Empty;
        private const int FileChecksumLength = 40;
        private const int RawFilePathLength = 200;
        private int _fileFormatId = FileFormatId; // For internal checks and backwards compatibility usages.
        /// <summary>
        /// The length of the Native ID field in the binary file
        /// </summary>
        protected internal const int NativeIdLength = 50;

        private readonly object _fileLock = new object();
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

        private readonly List<XicPoint> _precursorChromatogramCache = new List<XicPoint>();

        #endregion

        #region IMassSpecDataReader implementation

        /// <summary>
        /// Try to make the reader random access capable
        /// </summary>
        /// <returns>true if is random access capable, false if not</returns>
        public override bool TryMakeRandomAccessCapable()
        {
            return true;
        }

        /// <summary>
        /// Close the reader
        /// </summary>
        public override void Close()
        {
            _reader.Close();
        }

        /// <summary>
        /// Properly dispose of all unmanaged resources (specifically, file handles)
        /// </summary>
        public override void Dispose()
        {
            _reader.Dispose();
        }

        /// <summary>
        /// Path to the file; is <see cref="string.Empty"/> if the reader is in-memory
        /// </summary>
        public override string FilePath
        {
            get => PbfFilePath;
            protected set
            {
                // DO NOT USE!!!
            }
        }

        /// <summary>
        /// SHA-1 Checksum of the original input file (raw, mzML, .d folder, etc.)
        /// </summary>
        public override string SrcFileChecksum { get; protected set; }

        /// <summary>
        /// Version of the immediate prior input file (raw, mzML, .d folder, etc.)
        /// </summary>
        public override string FileFormatVersion { get; protected set; }

        #endregion

        #region LcMsRun Public function overrides

        /// <summary>
        /// The smallest MS1 m/z
        /// </summary>
        public override double MinMs1Mz => _minMs1Mz;

        /// <summary>
        /// The largest MS1 m/z
        /// </summary>
        public override double MaxMs1Mz => _maxMs1Mz;

        /// <summary>
        /// Get the native ID of the specified scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns></returns>
        public override string GetNativeId(int scanNum)
        {
            if (!ScanNumNativeIdMap.TryGetValue(scanNum, out var nativeId))
            {
                nativeId = string.Empty;
                var spec = GetSpectrum(scanNum, false);
                if (spec != null)
                {
                    nativeId = spec.NativeId;
                }
                ScanNumNativeIdMap.Add(scanNum, nativeId);
            }

            return nativeId;
        }

        /// <summary>
        /// Read the specified spectrum from the file, optionally reading only the metadata
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="includePeaks"></param>
        /// <returns></returns>
        public override Spectrum GetSpectrum(int scanNum, bool includePeaks = true)
        {
            if (!_scanNumToSpecOffset.TryGetValue(scanNum, out var offset))
            {
                return null;
            }

            var spec = ReadSpectrum(offset, includePeaks);
            if (spec.MsLevel == 1 && _precursorSignalToNoiseRatioThreshold > 0.0)
            {
                spec.FilterNoise(_precursorSignalToNoiseRatioThreshold);
            }
            else if (_productSignalToNoiseRatioThreshold > 0.0)
            {
                spec.FilterNoise(_productSignalToNoiseRatioThreshold);
            }

            return spec;
        }

        /// <summary>
        /// Read and return the isolation window for the specified scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns></returns>
        public override IsolationWindow GetIsolationWindow(int scanNum)
        {
            if (_scanNumToSpecOffset.TryGetValue(scanNum, out var offset))
            {
                if (ReadSpectrum(offset, false) is ProductSpectrum spec)
                {
                    return spec.IsolationWindow;
                }
            }
            return null;
        }

        /// <summary>
        /// Returns a XIC for the chosen range that covers the entire run.
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

            var scanToXicPoint = new XicPoint[NumSpectra + 1];
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
                var localXic = new Xic();
                localXic.AddRange(_precursorChromatogramCache.Where(peak => minMz <= peak.Mz && peak.Mz <= maxMz));
                return Xic.GetSelectedXic(localXic);
            }

            var minBinIndex = GetMzBinIndex(minMz);
            var maxBinIndex = GetMzBinIndex(maxMz);

            long targetOffset;
            if (minBinIndex == maxBinIndex)
            {
                if (maxBinIndex < _minMzIndex || maxBinIndex > _maxMzIndex)
                {
                    return new Xic();
                }

                var offset = _chromMzIndexToOffset[maxBinIndex - _minMzIndex];
                if (offset < _offsetPrecursorChromatogramBegin)
                {
                    return new Xic();
                }

                // binary search
                var beginOffset = offset;
                var endOffset = _chromMzIndexToOffset[maxBinIndex - _minMzIndex + 1];
                targetOffset = GetOffset(minMz, maxMz, beginOffset, endOffset);
            }
            else
            {
                if (maxBinIndex < _minMzIndex || minBinIndex > _maxMzIndex)
                {
                    return new Xic();
                }

                targetOffset = maxBinIndex > _maxMzIndex ? _offsetPrecursorChromatogramEnd : _chromMzIndexToOffset[maxBinIndex - _minMzIndex];
            }

            if (targetOffset < _offsetPrecursorChromatogramBegin)
            {
                return new Xic();
            }

            var xic = GetXic(minMz, maxMz, _offsetPrecursorChromatogramBegin, _offsetPrecursorChromatogramEnd, targetOffset);
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
                var localXic = new Xic();
                localXic.AddRange(_precursorChromatogramCache.Where(peak => minMz <= peak.Mz && peak.Mz <= maxMz));
                return localXic;
            }

            var minBinIndex = GetMzBinIndex(minMz);
            var maxBinIndex = GetMzBinIndex(maxMz);

            long targetOffset;
            if (minBinIndex == maxBinIndex)
            {
                if (maxBinIndex < _minMzIndex || maxBinIndex > _maxMzIndex)
                {
                    return new Xic();
                }

                var offset = _chromMzIndexToOffset[maxBinIndex - _minMzIndex];
                if (offset < _offsetPrecursorChromatogramBegin)
                {
                    return new Xic();
                }

                // binary search
                var beginOffset = offset;
                var endOffset = _chromMzIndexToOffset[maxBinIndex - _minMzIndex + 1];
                targetOffset = GetOffset(minMz, maxMz, beginOffset, endOffset);
            }
            else
            {
                if (maxBinIndex < _minMzIndex || minBinIndex > _maxMzIndex)
                {
                    return new Xic();
                }

                targetOffset = maxBinIndex > _maxMzIndex ? _offsetPrecursorChromatogramEnd : _chromMzIndexToOffset[maxBinIndex - _minMzIndex];
            }

            if (targetOffset < _offsetPrecursorChromatogramBegin)
            {
                return new Xic();
            }

            var xic = GetChromatogramRange(minMz, maxMz, _offsetPrecursorChromatogramBegin, _offsetPrecursorChromatogramEnd, targetOffset);
            if (!xic.Any())
            {
                return xic;
            }

            xic.Sort();
            return xic;
        }

        /// <summary>
        /// If <paramref name="scanNum"/> is a MS1 scan, return it; otherwise, return null.
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="ms1ScanIndex"></param>
        /// <returns></returns>
        public override Spectrum GetMs1Spectrum(int scanNum, out int ms1ScanIndex)
        {
            ms1ScanIndex = -1;
            if (!_scanNumToSpecOffset.TryGetValue(scanNum, out var offset))
            {
                return null;
            }

            var ms1ScanNums = GetMs1ScanVector();
            ms1ScanIndex = Array.BinarySearch(ms1ScanNums, scanNum);
            if (ms1ScanIndex < 0)
            {
                return null;
            }

            return ReadSpectrum(offset, true);
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
            FileFormatVersion = _fileFormatId.ToString();
            if (_fileFormatId > FileFormatId || _fileFormatId < EarliestSupportedFileFormatId)
            {
                return false;
            }

            // Backup 10 bytes
            _reader.BaseStream.Seek(-3 * sizeof(long) - 1 * sizeof(int), SeekOrigin.End);

            // Temporarily store position (used to read raw file name, if available)
            var offsetDataPos = _reader.BaseStream.Position;

            // Read the byte offset of the start of the precursor chromatogram
            _offsetPrecursorChromatogramBegin = _reader.ReadInt64();

            // Read the byte offset of the start of the product chromatogram (which is the end of the precursor chromatogram)
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
                if (_fileFormatId > 150605)
                {
                    scanNum = _reader.ReadInt32();
                }
                var msLevel = _reader.ReadInt32();
                if (msLevel < minMsLevel)
                {
                    minMsLevel = msLevel;
                }

                if (msLevel > maxMsLevel)
                {
                    maxMsLevel = msLevel;
                }

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
                        if (isoWindowSet.Count >= NumUniqueIsolationWindowThresholdForDia)
                        {
                            isDda = true;
                        }
                    }
                    var minBinNum = (int)Math.Round(minMz * IsolationWindowBinningFactor);
                    var maxBinNum = (int)Math.Round(maxMz * IsolationWindowBinningFactor);
                    for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
                    {
                        if (!isolationMzBinToScanNums.TryGetValue(binNum, out var scanNumList))
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

            _minMzIndex = _reader.ReadInt32();
            _maxMzIndex = _reader.ReadInt32();
            // _maxMzIndex is less than _minMzIndex if there is no MS1 chromatogram data
            if (_maxMzIndex >= _minMzIndex)
            {
                _chromMzIndexToOffset = new long[_maxMzIndex - _minMzIndex + 2];

                for (var i = 0; i < _chromMzIndexToOffset.Length; i++)
                {
                    _chromMzIndexToOffset[i] = _reader.ReadInt64();
                }

                _chromMzIndexToOffset[_chromMzIndexToOffset.Length - 1] = _offsetPrecursorChromatogramEnd;
            }
            else
            {
                _chromMzIndexToOffset = new long[0];
            }

            if (_fileFormatId > 150606)
            {
                var rawPathLocation = offsetDataPos - RawFilePathLength - sizeof(int) * 2;
                _reader.BaseStream.Seek(rawPathLocation, SeekOrigin.Begin);
                var c = new char[RawFilePathLength];
                _reader.Read(c, 0, RawFilePathLength);
                RawFilePath = (new string(c)).Trim();
                NativeFormat = (CV.CVID)_reader.ReadInt32();
                NativeIdFormat = (CV.CVID)_reader.ReadInt32();

                var checksumLocation = rawPathLocation - 40; // stored as hex, so 40 bytes rather than 20.
                if (_fileFormatId > 150607)
                {
                    _reader.BaseStream.Seek(checksumLocation, SeekOrigin.Begin);
                    c = new char[FileChecksumLength];
                    _reader.Read(c, 0, FileChecksumLength);
                    SrcFileChecksum = (new string(c)).Trim();
                    // This code works if the checksum is written as bytes instead of as a string.
                    //var bytes = _reader.ReadBytes(20);
                    //FileChecksum = BitConverter.ToString(bytes);
                }
            }
            else
            {
                RawFilePath = string.Empty;
                NativeIdFormat = CV.CVID.CVID_Unknown;
                NativeFormat = CV.CVID.CVID_Unknown;

                // Take a guess at the Native ID Format, based on the native ID format of the first spectrum...
                var nativeId = GetSpectrum(MinLcScan, false).NativeId;
                if (string.IsNullOrWhiteSpace(nativeId))
                {
                    // We don't know, and can't guess with any accuracy.
                    NativeIdFormat = CV.CVID.CVID_Unknown;
                    NativeFormat = CV.CVID.CVID_Unknown;
                }
                else if (nativeId.StartsWith("controllerType="))
                {
                    NativeIdFormat = CV.CVID.MS_Thermo_nativeID_format;
                    NativeFormat = CV.CVID.MS_Thermo_RAW_format;
                }
                else if (nativeId.StartsWith("scanId="))
                {
                    NativeIdFormat = CV.CVID.MS_Agilent_MassHunter_nativeID_format;
                    NativeFormat = CV.CVID.MS_Agilent_MassHunter_format;
                }
                else if (nativeId.StartsWith("scan="))
                {
                    // Could be one of multiple - BrukerAgilentYep, BrukerBaf, ScanNumberOnly
                    // Just going to set BrukerBaf, since it seems more likely at current time than BrukerAgilentYep or ScanNumberOnly
                    NativeIdFormat = CV.CVID.MS_Bruker_BAF_nativeID_format;
                    NativeFormat = CV.CVID.MS_Bruker_BAF_format;
                }
                else if (nativeId.StartsWith("frame="))
                {
                    NativeIdFormat = CV.CVID.MS_UIMF_nativeID_format;
                    NativeFormat = CV.CVID.MS_UIMF_format;
                }
                else if (nativeId.StartsWith("function="))
                {
                    NativeIdFormat = CV.CVID.MS_Waters_nativeID_format;
                    NativeFormat = CV.CVID.MS_Waters_raw_format;
                }
                else if (nativeId.StartsWith("sample="))
                {
                    NativeIdFormat = CV.CVID.MS_WIFF_nativeID_format;
                    NativeFormat = CV.CVID.MS_ABI_WIFF_format;
                }
                else if (nativeId.StartsWith("index="))
                {
                    NativeIdFormat = CV.CVID.MS_multiple_peak_list_nativeID_format;
                    NativeFormat = CV.CVID.MS_Mascot_MGF_format;
                }
                // Include others if deemed necessary...
            }

            return true;
        }

        #endregion

        #region Read/Write Spectrum

        private Spectrum ReadSpectrum(long offset, bool includePeaks = true)
        {
            lock (_fileLock)
            {
                _reader.BaseStream.Seek(offset, SeekOrigin.Begin);
                while (_reader.BaseStream.Position != (_reader.BaseStream.Length - sizeof(int)))
                {
                    var spec = ReadSpectrum(_reader, includePeaks);
                    return spec;
                }
                return null;
            }
        }

        /// <summary>
        /// Read a spectrum from the current position in <paramref name="reader"/>, with the option to only read the metadata.
        /// </summary>
        /// <param name="reader"></param>
        /// <param name="includePeaks"></param>
        /// <returns></returns>
        protected internal virtual Spectrum ReadSpectrum(BinaryReader reader, bool includePeaks = true)
        {
            // Must reflect all changes to WriteSpectrum
            // reader is passed in because this function should be called within a lock to prevent reading/writing errors

            var scanNum = reader.ReadInt32();
            var nativeId = string.Empty;
            if (_fileFormatId > 150604)
            {
                var c = new char[NativeIdLength];
                reader.Read(c, 0, NativeIdLength);
                nativeId = (new string(c)).Trim();
            }
            var msLevel = reader.ReadByte();
            var elutionTime = reader.ReadDouble();
            double tic = -1;
            if (_fileFormatId > 150604)
            {
                tic = reader.ReadSingle();
            }

            Spectrum spec;

            double calcTic;
            if (msLevel > 1)
            {
                double? precursorMass = reader.ReadDouble();
                if (Math.Abs(precursorMass.Value) < float.Epsilon)
                {
                    precursorMass = null;
                }

                int? precursorCharge = reader.ReadInt32();
                if (precursorCharge == 0)
                {
                    precursorCharge = null;
                }

                var activationMethod = (ActivationMethod)reader.ReadByte();
                var isolationWindowTargetMz = reader.ReadDouble();
                var isolationWindowLowerOffset = reader.ReadDouble();
                var isolationWindowUpperOffset = reader.ReadDouble();
                var peakList = ReadPeakList(reader, out calcTic, includePeaks);
                if (tic < 0)
                {
                    tic = calcTic;
                }
                spec = new ProductSpectrum(peakList, scanNum)
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
                var peakList = ReadPeakList(reader, out calcTic, includePeaks);
                if (tic < 0)
                {
                    tic = calcTic;
                }
                spec = new Spectrum(peakList, scanNum);
            }
            spec.MsLevel = msLevel;
            spec.ElutionTime = elutionTime;
            spec.NativeId = nativeId;
            spec.TotalIonCurrent = tic;
            return spec;
        }

        /// <summary>
        /// Write the supplied spectrum to the current position in <paramref name="writer"/>
        /// </summary>
        /// <param name="spec"></param>
        /// <param name="writer"></param>
        protected internal virtual void WriteSpectrum(Spectrum spec, BinaryWriter writer)
        {
            // All changes made here must be duplicated to ReadSpectrum() and GetPeakMetadataForSpectrum()
            // writer is passed in because this function should be called within a lock to prevent reading/writing errors

            // scan number: 4
            writer.Write(spec.ScanNum);

            // NativeID: 50
            // pad or truncate to keep in limit (may have to change in future...)
            writer.Write(spec.NativeId.PadRight(NativeIdLength).ToCharArray(0, NativeIdLength), 0, NativeIdLength);

            // ms level: 1
            writer.Write(Convert.ToByte(spec.MsLevel));

            // elution time: 8
            writer.Write(spec.ElutionTime);

            // Total Ion Current: 4
            writer.Write(Convert.ToSingle(spec.TotalIonCurrent));

            if (spec is ProductSpectrum productSpec)    // product spectrum
            {
                var isolationWindow = productSpec.IsolationWindow;
                // precursor mass: 8
                writer.Write(isolationWindow.MonoisotopicMz ?? 0.0);
                // precursor charge: 4
                writer.Write(isolationWindow.Charge ?? 0);
                // Activation method: 1
                writer.Write((byte)productSpec.ActivationMethod);
                // Isolation window target m/z: 8
                writer.Write(isolationWindow.IsolationWindowTargetMz);
                // Isolation window lower offset: 8
                writer.Write(isolationWindow.IsolationWindowLowerOffset);
                // Isolation window upper offset: 8
                writer.Write(isolationWindow.IsolationWindowUpperOffset);
            }

            // Guarantee sorted peaks.
            Array.Sort(spec.Peaks);

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

        // Must reflect all changes to WriteSpectrum
        private ScanPeakMetaData GetPeakMetaDataForSpectrum(int scanNum)
        {
            if (!_scanNumToSpecOffset.TryGetValue(scanNum, out var offset))
            {
                return null;
            }

            ScanPeakMetaData data;
            lock (_fileLock)
            {
                _reader.BaseStream.Seek(offset, SeekOrigin.Begin);

                var rScanNum = _reader.ReadInt32();
                // skip nativeId
                _reader.BaseStream.Seek(NativeIdLength, SeekOrigin.Current);
                var msLevel = _reader.ReadByte();

                // ReSharper disable UnusedVariable
                var elutionTime = _reader.ReadDouble();
                var tic = _reader.ReadSingle();
                // ReSharper restore UnusedVariable

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

        private List<Peak> ReadPeakList(BinaryReader reader, out double tic, bool includePeaks = true)
        {
            var peakList = new List<Peak>();
            var numPeaks = reader.ReadInt32();
            // Only used if fileFormatId < 150605
            tic = 0;

            // Skip the read if peaks aren't requested
            if (!includePeaks)
            {
                // first the number of peaks, then 12 bytes per peak (mz, double, 8 bytes, then intensity, single, 4 bytes)
                if (_fileFormatId < 150605)
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
                peakList.Add(new Peak(mz, intensity));
            }
            return peakList;
        }

        #endregion

        #region WriteAsPbf (obsolete)

        /// <summary>
        /// Old PBF file creation workflow
        /// </summary>
        /// <param name="lcmsRun"></param>
        /// <param name="outputFilePath"></param>
        /// <param name="progress"></param>
        [Obsolete("Use PbfLcMsRun(string, IMassSpecDataReader, ...) for an optimized pbf creation process", true)]
        public static void WriteAsPbf(InMemoryLcMsRun lcmsRun, string outputFilePath, IProgress<ProgressData> progress = null)
        {
            using (var writer = new BinaryWriter(File.Open(outputFilePath, FileMode.Create)))
            {
                WriteAsPbf(lcmsRun, writer, progress);
            }
        }

        /// <summary>
        /// Old PBF file creation workflow
        /// </summary>
        /// <param name="lcmsRun"></param>
        /// <param name="writer"></param>
        /// <param name="progress"></param>
        [Obsolete("Use PbfLcMsRun(string, IMassSpecDataReader, ...) for an optimized pbf creation process", true)]
        public static void WriteAsPbf(InMemoryLcMsRun lcmsRun, BinaryWriter writer, IProgress<ProgressData> progress = null)
        {
            var pbfLcMsRun = new PbfLcMsRun();
            var progressData = new ProgressData(progress);

            var scanNumToSpecOffset = new long[lcmsRun.NumSpectra + 1];
            var scanNumToIsolationWindow = new IsolationWindow[lcmsRun.NumSpectra + 1];

            // Spectra
            long countTotal = lcmsRun.NumSpectra;
            long counter = 0;
            progressData.StepRange(42.9, "Writing spectra data"); // SpecData: Approximately 43% of total file size
            long countMS2Spec = 0;
            for (var scanNum = lcmsRun.MinLcScan; scanNum <= lcmsRun.MaxLcScan; scanNum++)
            {
                progressData.Report(counter, countTotal);
                counter++;
                scanNumToSpecOffset[scanNum - lcmsRun.MinLcScan] = writer.BaseStream.Position;
                var spec = lcmsRun.GetSpectrum(scanNum);
                if (spec == null)
                {
                    continue;
                }

                scanNumToIsolationWindow[scanNum - lcmsRun.MinLcScan] = null;
                if (spec is ProductSpectrum productSpec)
                {
                    scanNumToIsolationWindow[scanNum - lcmsRun.MinLcScan] = productSpec.IsolationWindow;
                    countMS2Spec++;
                }
                pbfLcMsRun.WriteSpectrum(spec, writer);
            }

            // Precursor ion chromatogram (MS1 spectra)
            var offsetBeginPrecursorChromatogram = writer.BaseStream.Position;

            var minMzIndex = lcmsRun.Ms1PeakList.Count > 0 ? GetMzBinIndex(lcmsRun.Ms1PeakList[0].Mz) : 0;
            var maxMzIndex = lcmsRun.Ms1PeakList.Count > 0 ? GetMzBinIndex(lcmsRun.Ms1PeakList[lcmsRun.Ms1PeakList.Count - 1].Mz) : -1;

            var chromMzIndexToOffset = new long[maxMzIndex - minMzIndex + 1];
            var prevMzIndex = -1;
            counter = 0;
            countTotal = lcmsRun.Ms1PeakList.Count;

            // All MS1 data sorted by mass, then scan number
            // In rare instances, lcmsRun.Ms1PeakList will be blank (no MS1 spectra); that's OK
            progressData.StepRange(42.9 + 15.7, "Writing precursor chromatogram"); // Approximately 16% of total file size
            foreach (var peak in lcmsRun.Ms1PeakList)
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
            progressData.StepRange(42.9 + 15.7 + (41.2 / 2), "Processing product ion chromatogram"); // Approximately 41% of total file size
            foreach (var ms2ScanNum in lcmsRun.GetScanNumbers(2))
            {
                progressData.Report(counter, countTotal);
                counter++;
                if (!(lcmsRun.GetSpectrum(ms2ScanNum) is ProductSpectrum productSpec))
                {
                    continue;
                }

                foreach (var peak in productSpec.Peaks)
                {
                    ms2PeakList.Add(new LcMsPeak(peak.Mz, peak.Intensity, ms2ScanNum));
                }
            }
            ms2PeakList.Sort();

            var offsetBeginProductChromatogram = writer.BaseStream.Position;
            counter = 0;
            countTotal = ms2PeakList.Count;
            progressData.StepRange(42.9 + 15.7 + 41.2, "Writing product ion chromatogram"); // Approximately 41% of total file size
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
            progressData.IsPartialRange = false;
            progressData.Report(99.8, "Writing metadata"); // Metadata: Approximately 0.2% of total file size

            var warnedInvalidScanNum = false;
            var warnedNullScanToIsolationWindow = false;

            writer.Write(lcmsRun.MinLcScan);
            writer.Write(lcmsRun.MaxLcScan);
            for (var scanNum = lcmsRun.MinLcScan; scanNum <= lcmsRun.MaxLcScan; scanNum++)
            {
                var msLevel = lcmsRun.GetMsLevel(scanNum);
                writer.Write(lcmsRun.GetMsLevel(scanNum));
                writer.Write(lcmsRun.GetElutionTime(scanNum));

                if (msLevel == 2)
                {
                    float minMz = 0;
                    float maxMz = 0;

                    if (scanNum - lcmsRun.MinLcScan < 0 || scanNum - lcmsRun.MinLcScan >= scanNumToIsolationWindow.Length)
                    {
                        if (!warnedInvalidScanNum)
                        {
                            Console.WriteLine();
                            ConsoleMsgUtils.ShowWarning(string.Format(
                                "WriteAsPbf encountered an invalid scan number: {0}; " +
                                "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown", scanNum));
                            warnedInvalidScanNum = true;
                        }
                    }
                    else
                    {
                        if (scanNumToIsolationWindow[scanNum - lcmsRun.MinLcScan] == null)
                        {
                            if (!warnedNullScanToIsolationWindow)
                            {
                                Console.WriteLine();
                                ConsoleMsgUtils.ShowWarning(string.Format(
                                    "WriteAsPbf encountered a Null entry in scanNumToIsolationWindow for scan {0}; " +
                                    "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown", scanNum));
                                warnedNullScanToIsolationWindow = true;
                            }
                        }
                        else
                        {
                            minMz = (float)scanNumToIsolationWindow[scanNum - lcmsRun.MinLcScan].MinMz;
                            maxMz = (float)scanNumToIsolationWindow[scanNum - lcmsRun.MinLcScan].MaxMz;
                        }
                    }

                    writer.Write(minMz);
                    writer.Write(maxMz);
                }
                writer.Write(scanNumToSpecOffset[scanNum - lcmsRun.MinLcScan]);
            }

            // Precursor chromatogram index
            writer.Write(minMzIndex);   // min index
            writer.Write(maxMzIndex);
            progressData.Report(99.9); // Metadata: Approximately 0.2% of total file size

            var prevOffset = offsetBeginMetaInformation;
            for (var i = chromMzIndexToOffset.Length - 1; i >= 0; i--)
            {
                if (chromMzIndexToOffset[i] < offsetBeginPrecursorChromatogram)
                {
                    chromMzIndexToOffset[i] = prevOffset;
                }
                else
                {
                    prevOffset = chromMzIndexToOffset[i];
                }
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

        /// <summary>
        /// Bulk of code to write a PBF file. Protected internal to support DPbfLcMsRun.
        /// </summary>
        /// <param name="msDataReader"></param>
        /// <param name="writer"></param>
        /// <param name="startScan"></param>
        /// <param name="endScan"></param>
        /// <param name="progress"></param>
        private void WriteToPbf(
            IMassSpecDataReader msDataReader,
            BinaryWriter writer,
            int startScan,
            int endScan,
            IProgress<ProgressData> progress = null)
        {
            ScanNumToMsLevel = new Dictionary<int, int>(msDataReader.NumSpectra + 1);
            ScanNumElutionTimeMap = new Dictionary<int, double>(msDataReader.NumSpectra + 1);
            IsolationMzBinToScanNums = new Dictionary<int, int[]>();
            _scanNumToSpecOffset = new Dictionary<int, long>(msDataReader.NumSpectra + 1);

            MinLcScan = int.MaxValue;
            MaxLcScan = int.MinValue;
            MinMsLevel = int.MaxValue;
            MaxMsLevel = int.MinValue;

            var progressData = new ProgressData(progress);

            var scanNumToIsolationWindow = new Dictionary<int, IsolationWindow>(msDataReader.NumSpectra + 1);
            var ms1Scans = new List<int>();
            var ms2Scans = new List<int>();

            // Spectra
            long ms1PeakCount = 0;
            long ms2PeakCount = 0;
            var maxMs1Mz = double.MinValue;
            var minMs1Mz = double.MaxValue;
            var maxMs2Mz = double.MinValue;
            var minMs2Mz = double.MaxValue;
            var scanMetadata = new List<ScanMetadata>(msDataReader.NumSpectra);
            int countTotal;

            if (endScan > 0 && endScan > msDataReader.NumSpectra)
            {
                endScan = 0;
            }

            if (startScan > 0)
            {
                if (endScan > 0)
                {
                    countTotal = endScan - startScan + 1;
                }
                else
                {
                    countTotal = msDataReader.NumSpectra - startScan + 1;
                }
            }
            else if (endScan > 0)
            {
                countTotal = endScan;
            }
            else
            {
                countTotal = msDataReader.NumSpectra;
            }

            long counter = 0;
            progressData.StepRange(42.9, "Writing spectra data"); // SpecData: Approximately 43% of total file size
            foreach (var spec in msDataReader.ReadAllSpectra())
            {
                if (startScan > 0 && spec.ScanNum < startScan)
                {
                    continue;
                }

                if (endScan > 0 && spec.ScanNum > endScan)
                {
                    break;
                }

                progressData.Report(counter, countTotal);
                counter++;

                // Store offset, and write spectrum now
                _scanNumToSpecOffset.Add(spec.ScanNum, writer.BaseStream.Position);
                WriteSpectrum(spec, writer);

                // Handle other metadata stuff.
                var maxMz = double.MinValue;
                var minMz = double.MaxValue;

                if (spec.Peaks.Length > 0)
                {
                    minMz = spec.Peaks[0].Mz;
                    maxMz = spec.Peaks[spec.Peaks.Length - 1].Mz;
                }
                scanNumToIsolationWindow[spec.ScanNum] = null;

                if (spec is ProductSpectrum productSpec)
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

            if (ms1PeakCount > 0 && ContainsChromatograms)
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

            if (ms2PeakCount > 0 && ContainsChromatograms)
            {
                progressData.StepRange(42.9 + 15.7 + 41.2, "Writing product chromatogram"); // Approximately 41% of total file size, on standard LCMS file
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
            progressData.IsPartialRange = false;
            progressData.Report(99.8, "Writing metadata"); // Metadata: Approximately 0.2% of total file size

            var warnedInvalidScanNum = false;
            var warnedNullScanToIsolationWindow = false;

            writer.Write(MinLcScan);
            writer.Write(MaxLcScan);
            scanMetadata.Sort();

            var isoWindowSet = new HashSet<IsolationWindow>();
            var isDda = false;
            var isolationMzBinToScanNums = new Dictionary<int, List<int>>();

            MinMsLevel = scanMetadata.Min(scan => scan.MsLevel);
            MaxMsLevel = scanMetadata.Max(scan => scan.MsLevel);
            foreach (var scan in scanMetadata)
            {
                var msLevel = scan.MsLevel;
                writer.Write(scan.ScanNum);
                writer.Write(scan.MsLevel);
                writer.Write(scan.ElutionTime);
                ScanNumToMsLevel[scan.ScanNum] = scan.MsLevel;
                ScanNumElutionTimeMap[scan.ScanNum] = scan.ElutionTime;

                if (msLevel == 2)
                {
                    float minMz = 0;
                    float maxMz = 0;

                    if (!scanNumToIsolationWindow.ContainsKey(scan.ScanNum))
                    {
                        if (!warnedInvalidScanNum)
                        {
                            Console.WriteLine();
                            ConsoleMsgUtils.ShowWarning(string.Format(
                                "WriteAsPbf encountered an invalid scan number: {0}; " +
                                "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown", scan.ScanNum));
                            warnedInvalidScanNum = true;
                        }
                    }
                    else
                    {
                        if (scanNumToIsolationWindow[scan.ScanNum] == null)
                        {
                            if (!warnedNullScanToIsolationWindow)
                            {
                                Console.WriteLine();
                                ConsoleMsgUtils.ShowWarning(string.Format(
                                    "WriteAsPbf encountered a Null entry in scanNumToIsolationWindow for scan {0}; " +
                                    "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown", scan.ScanNum));
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
                        if (isoWindowSet.Count >= NumUniqueIsolationWindowThresholdForDia)
                        {
                            isDda = true;
                        }
                    }
                    var minBinNum = (int)Math.Round(minMz * IsolationWindowBinningFactor);
                    var maxBinNum = (int)Math.Round(maxMz * IsolationWindowBinningFactor);
                    for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
                    {
                        if (!isolationMzBinToScanNums.TryGetValue(binNum, out var scanNumList))
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
            if (ms1PeakCount > 0 && ContainsChromatograms)
            {
                for (var i = _chromMzIndexToOffset.Length - 2; i >= 0; i--)
                {
                    if (_chromMzIndexToOffset[i] < _offsetPrecursorChromatogramBegin)
                    {
                        _chromMzIndexToOffset[i] = prevOffset;
                    }
                    else
                    {
                        prevOffset = _chromMzIndexToOffset[i];
                    }
                }

                foreach (var offset in _chromMzIndexToOffset.Take(_chromMzIndexToOffset.Length - 1))
                {
                    writer.Write(offset);
                }

                _chromMzIndexToOffset[_chromMzIndexToOffset.Length - 1] = _offsetPrecursorChromatogramEnd;
            }

            // Checksum: 40 bytes (could store in 20 bytes, but conversion from hex string to bytes isn't simple)
            SrcFileChecksum = msDataReader.SrcFileChecksum.ToLower().Replace("-", "");
            writer.Write(SrcFileChecksum.PadRight(FileChecksumLength).ToCharArray(0, FileChecksumLength), 0, FileChecksumLength);

            // RawFilePath: 200
            // pad or truncate to keep in limit (may have to change in future...)
            if (RawFilePath.Length <= RawFilePathLength)
            {
                // Length is within limit. Write the whole thing with padding.
                writer.Write(RawFilePath.PadRight(RawFilePathLength).ToCharArray(0, RawFilePathLength), 0, RawFilePathLength);
            }
            else
            {
                // Length is beyond limit. Write only the file name, with padding.
                writer.Write(Path.GetFileName(RawFilePath).PadRight(RawFilePathLength).ToCharArray(0, RawFilePathLength), 0, RawFilePathLength);
            }
            writer.Write((int)NativeFormat);
            writer.Write((int)NativeIdFormat);

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
            var maxInMemoryPerSpec = (int)(memFreeKB * 1024 / 2 / scansCount / (20 + 8));
            if (maxInMemoryPerSpec * scansCount > max) // Set a hard limit at 10 millions peaks in memory at once (only exception is the minimum)
            {
                maxInMemoryPerSpec = max / scansCount;
            }

            if (maxInMemoryPerSpec < 5)
            {
                maxInMemoryPerSpec = 5;
            }

            return maxInMemoryPerSpec;
        }

        private void CreateAndOutputMsXChromatogram_Merge(
            BinaryWriter writer,
            IReadOnlyCollection<int> scansForMsLevelX,
            long totalPeaksCount,
            bool isMs1List,
            double minMz,
            double maxMz,
            IProgress<ProgressData> progress)
        {
            // Other thought for a slower, but really low-memory chromatogram creator:
            //   Make sure peaks are sorted ascending when written, and then jump through all spectra, performing a massive merge sort on the peaks and outputting the lowest spectra

            // Limit progress reporting to 5 times per second - reporting progress for every output peak may overload consumers of the progress reporting
            // Overload was particularly noted in LcMsSpectator, causing UI lockups and extreme memory usage.
            const double minTimeForProgressSeconds = 0.2;

            var progData = new ProgressData(progress);
            var count = 0;
            var countTotal = scansForMsLevelX.Count;

            // List overhead per item: is array-backed, so item (+ pointer to it)
            // item size: double, double, int, so 20 bytes
            var memoryFreeKB = SystemInfo.GetFreeMemoryMB() / 1024.0;
            var totalMemKB = SystemInfo.GetTotalMemoryMB() / 1024.0;

            // Configure a reserve amount to avoid using all physical memory, which has the cost of excessive paging
            var quarterTotalPhysicalMemory = totalMemKB / 4;
            // Cut the reserve down on systems with large amounts of physical memory (16GB and greater)
            while (quarterTotalPhysicalMemory > 4194304)
            {
                quarterTotalPhysicalMemory /= 2;
            }
            var memFreeLessReserve = memoryFreeKB - quarterTotalPhysicalMemory;

            // Use a custom split-list implementation: Reduce the number of items to sort each time
            // The enumerator access is O(1)
            // 8 byte (pointer) overhead per item, plus ~30 bytes per list (number of lists is numSpectra / 1,000,000 + 2)
            // item size: double, double, int, so 20 bytes
            // Perform a massive merge-sort style read/write, to minimize memory usage - set a minimum of 5
            var maxInMemoryPerSpec = CalculateMaxSpecsInMem(memFreeLessReserve, scansForMsLevelX.Count);

            progData.Status = "Writing product chromatogram";
            if (isMs1List)
            {
                progData.Status = "Writing precursor chromatogram";
            }
            lock (_fileLock)
            {
                var peaks = new SplitLcMsPeakLists(minMz, maxMz, totalPeaksCount);
                var peaksCount = 0;
                var metadata = new Dictionary<int, ScanPeakMetaData>(scansForMsLevelX.Count);

                if (_reader == null)
                {
                    _reader = new BinaryReader(File.Open(PbfFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));
                }

                progData.StepRange(5);
                var lastProgressTime = DateTime.MinValue;
                // Read in metadata, and the first "maxInMemoryPerSpec" peaks of each scan (min 5, max 25000000 / numSpectra)
                foreach (var scan in scansForMsLevelX)
                {
                    if (lastProgressTime.AddSeconds(minTimeForProgressSeconds) < DateTime.Now)
                    {
                        progData.Report(count, countTotal);
                        lastProgressTime = DateTime.Now;
                    }
                    count++;
                    var datum = GetPeakMetaDataForSpectrum(scan);
                    peaksCount += datum.NumPeaks;
                    metadata.Add(datum.ScanNum, datum);
                    peaks.AddRange(datum.ReadPeaks(_reader, maxInMemoryPerSpec));
                }

                progData.Report(100);

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
                lastProgressTime = DateTime.MinValue;
                count = 0;
                var prevMzIndex = -1;
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
                        if (lastProgressTime.AddSeconds(minTimeForProgressSeconds) < DateTime.Now)
                        {
                            progData.Report(count, peaksCount);
                            lastProgressTime = DateTime.Now;
                        }
                        count++;

                        if (peak.Mz.Equals(0))
                        {
                            continue;
                        }

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
                        writer.Write((float)peak.Intensity);
                        writer.Write(peak.ScanNum);

                        var datum = metadata[peak.ScanNum];
                        datum.ReadOne();

                        if (datum.Count < 1)
                        {
                            if (datum.MorePeaksToRead)
                            {
                                break;
                            }

                            metadata.Remove(datum.ScanNum);
                        }
                    }
                    maxInMemoryPerSpec = CalculateMaxSpecsInMem(memFreeLessReserve, metadata.Count);

                    // add new entries back into the list from the spectra that the peaks came from
                    foreach (var datum in metadata.Values)
                    {
                        if (datum.Count < maxInMemoryPerSpec && datum.MorePeaksToRead)
                        {
                            peaks.AddRange(datum.ReadPeaks(_reader, maxInMemoryPerSpec - datum.Count));
                        }
                    }

                    peaksEnum = peaks.GetReusableEnumerator();
                }

                progData.Report(100);

                if (isMs1List)
                {
                    _chromMzIndexToOffset = new long[_maxMzIndex - _minMzIndex + 2];
                    for (var i = 0; i < _chromMzIndexToOffset.Length; i++)
                    {
                        if (chromMzIndexToOffset.ContainsKey(i))
                        {
                            _chromMzIndexToOffset[i] = chromMzIndexToOffset[i];
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Basically a struct version of <see cref="LcMsPeak"/> for avoiding some class memory allocation/collection issues, and pointer overhead
        /// </summary>
        /// <remarks>When using a class, there is a very large number of small objects placed in the generation 2 heap,
        /// which requires a blocking garbage collection to clear. The first workaround was to re-use peak objects to
        /// minimize the number of objects allocated (and necessary garbage collections), but that still tended to have
        /// some post-PBF creation memory usage issues. Changing to a struct significantly reduces those issues, even
        /// though these structs are slightly larger than Microsoft's recommendation at
        /// https://docs.microsoft.com/en-us/dotnet/standard/design-guidelines/choosing-between-class-and-struct</remarks>
        private readonly struct ChromatogramPeak : IComparable<ChromatogramPeak>
        {
            public ChromatogramPeak(double mz, double intensity, int scanNum)
            {
                Mz = mz;
                Intensity = intensity;
                ScanNum = scanNum;
            }

            public double Mz { get; }
            public double Intensity { get; }
            public int ScanNum { get; }

            public int CompareTo(ChromatogramPeak other)
            {
                // Sort by mz, then by scan number
                var mzCompare = Mz.CompareTo(other.Mz);
                if (mzCompare == 0)
                {
                    // Force a stable sort
                    return ScanNum.CompareTo(other.ScanNum);
                }
                return mzCompare;
            }
        }

        private class SplitLcMsPeakLists
        {
            private readonly List<List<ChromatogramPeak>> _lists;
            private readonly List<bool> _sorted;
            private readonly int _mod;
            private readonly double _mzMod;
            private const int Divisor = 1000000;
            private readonly double _minMz;

            public SplitLcMsPeakLists(double minMz, double maxMz, long peakCount)
            {
                _minMz = minMz;
                _mzMod = 0;
                _mod = (int)(peakCount / Divisor + 1);
                _mzMod = (maxMz - _minMz) / (_mod - 1);
                _lists = new List<List<ChromatogramPeak>>(_mod + 2);
                _sorted = new List<bool>(_mod + 2);
                for (var i = 0; i <= _mod + 1; i++)
                {
                    _lists.Add(new List<ChromatogramPeak>(1000));
                    _sorted.Add(false);
                }
            }

            private void Add(ChromatogramPeak peak)
            {
                if (peak.Mz.Equals(0))
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

            public void AddRange(IEnumerable<ChromatogramPeak> peaks)
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
                    for (var i = 0; i < _lists.Count && _enumerated > 0; i++)
                    {
                        if (_lists[i].Count == 0)
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

            private long _enumerated;

            public IEnumerator<ChromatogramPeak> GetReusableEnumerator()
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
                        _lists[i] = new List<ChromatogramPeak>(10);
                    }
                }
            }
        }

        private class ScanMetadata : IComparable<ScanMetadata>
        {
            public int ScanNum { get; }
            public int MsLevel { get; }
            public double ElutionTime { get; }

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
            public int ScanNum { get; }
            public int NumPeaks { get; }
            private int PeaksRead { get; set; }

            private long _nextPeakOffset;
            public int Count { get; private set; }

            public bool MorePeaksToRead => PeaksRead < NumPeaks;

            public ScanPeakMetaData(int scanNum, BinaryReader reader)
            {
                Count = 0;
                ScanNum = scanNum;
                NumPeaks = reader.ReadInt32();
                PeaksRead = 0;
                _nextPeakOffset = reader.BaseStream.Position;
            }

            public IEnumerable<ChromatogramPeak> ReadPeaks(BinaryReader reader, int numPeaksToRead)
            {
                if (PeaksRead >= NumPeaks)
                {
                    yield return new ChromatogramPeak();
                }
                reader.BaseStream.Seek(_nextPeakOffset, SeekOrigin.Begin);

                for (var i = 0; i < numPeaksToRead && PeaksRead < NumPeaks; i++)
                {
                    var mz = reader.ReadDouble();
                    var intensity = reader.ReadSingle();
                    Count++;
                    PeaksRead++;
                    yield return new ChromatogramPeak(mz, intensity, ScanNum);
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

        /// <summary>
        /// Get the MzBin index for the supplied m/z
        /// </summary>
        /// <param name="mz"></param>
        /// <returns></returns>
        public static int GetMzBinIndex(double mz)
        {
            return (int)(mz / MzBinSize);
        }

        private long GetOffset(double minMz, double maxMz, long beginOffset, long endOffset)
        {
            var minOffset = beginOffset;
            var maxOffset = endOffset;
            var curOffset = -1L;
            lock (_fileLock)
            {
                // binary search
                while (minOffset <= maxOffset)
                {
                    curOffset = minOffset + (maxOffset - minOffset) / NumBytePeak / 2 * NumBytePeak;
                    _reader.BaseStream.Seek(curOffset, SeekOrigin.Begin);
                    var curMz = _reader.ReadDouble();
                    if (curMz < minMz)
                    {
                        minOffset = curOffset + NumBytePeak;
                    }
                    else if (curMz > maxMz)
                    {
                        maxOffset = curOffset - NumBytePeak;
                    }
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
            if (!xic.Any())
            {
                return xic;
            }

            return Xic.GetSelectedXic(xic);
        }

        // beginOffset: inclusive, endOffset: exclusive
        private Xic GetChromatogramRange(double minMz, double maxMz, long beginOffset, long endOffset, long targetOffset)
        {
            var xic = GetXicPointsWithin(minMz, maxMz, beginOffset, endOffset, targetOffset);
            if (!xic.Any())
            {
                return xic;
            }

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

            lock (_fileLock)
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

            _precursorChromatogramCache.Sort((x, y) => x.Mz.CompareTo(y.Mz));

            return xic;
        }

        #endregion
    }
}
