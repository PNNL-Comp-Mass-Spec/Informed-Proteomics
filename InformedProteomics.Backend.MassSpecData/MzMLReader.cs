using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text;
using System.Xml;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;
using PRISM;
using PSI_Interface.CV;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Reader for mzML files
    /// </summary>
    /// <remarks>
    /// Can handle gzipped mzML files, and read in a forward-only fashion or in a random-access fashion
    /// </remarks>
    public sealed class MzMLReader : IMassSpecDataReader
    {
        // ReSharper disable CommentTypo

        // Ignore Spelling: bal, bool, centroided, dta, endian, foreach, gzipped, interchannel, multiphoton
        // Ignore Spelling: photodissociation, pkl, readonly, referenceable, struct, wiff, xcalibur, xsi, zlib

        // Ignore Spelling: Biosystems, Biotech, Bioworks, Bruker, Micromass, Phenyx, Proteinscape, Sciex, Shimadzu

        // ReSharper restore CommentTypo

        #region Private Members

        private readonly string _filePath;
        private string _srcFileChecksum = string.Empty;
        private string _fileFormatVersion = string.Empty;
        private CV.CVID _nativeIdFormat = CV.CVID.CVID_Unknown;
        private CV.CVID _nativeFormat = CV.CVID.CVID_Unknown;
        private Stream _file;
        private StreamReader _fileReader;
        private XmlReader _xmlReaderForYield;
        private bool _reduceMemoryUsage;
        //private long _artificialScanNum = 1;
        private long _numSpectra = -1;
        private readonly IndexList _spectrumOffsets = new IndexList() { IndexType = IndexList.IndexListType.Spectrum };
        private readonly IndexList _chromatogramOffsets = new IndexList() { IndexType = IndexList.IndexListType.Chromatogram };
        private long _indexListOffset;
        private bool _haveIndex;
        private bool _haveMetaData;
        private bool _isGzipped;
        private string _unzippedFilePath;
        private bool _randomAccess;
        private bool _allRead;
        private readonly XmlReaderSettings _xSettings = new XmlReaderSettings { IgnoreWhitespace = true };
        private Encoding _encoding;
        private readonly List<Spectrum> _spectra = new List<Spectrum>();

        #endregion

        #region Internal Objects

        /// <summary>
        /// Enumeration of common mzML versions
        /// </summary>
        private enum MzML_Version
        {
            mzML1_0_0,
            mzML1_1_0
        }

        /// <summary>
        /// Store the mzML version, so that we can use it to adjust how some things are processed
        /// </summary>
        private MzML_Version _version;

        private class ScanData
        {
            public double MonoisotopicMz;
            public double StartTime;

            public ScanData()
            {
                MonoisotopicMz = 0.0;
                StartTime = 0.0;
            }
        }

        private enum ParamType
        {
            cvParam,
            userParam
        }

        /// <summary>
        /// Parameter base class
        /// </summary>
        private abstract class Param
        {
            protected ParamType ParamType { get; set; }

            /// <summary>
            /// Parameter name
            /// </summary>
            /// <remarks>Required</remarks>
            public string Name;

            /// <summary>
            /// Parameter value
            /// </summary>
            /// <remarks>Optional</remarks>
            public string Value;

            /// <summary>
            /// Parameter unit name
            /// </summary>
            /// <remarks>Optional</remarks>
            public string UnitName;

            /// <summary>
            /// Parameter unit CVRef
            /// </summary>
            /// <remarks>Optional</remarks>
            public string UnitCVRef;

            /// <summary>
            /// Parameter accession
            /// </summary>
            /// <remarks>Optional</remarks>
            public string UnitAccession;

            protected string _cvRef;
            protected string _accession;
            protected string _type;

            public virtual string CVRef
            {
                get => string.Empty;
                set => _cvRef = string.Empty;
            }

            public virtual string Accession
            {
                get => string.Empty;
                set => _accession = string.Empty;
            }

            public virtual string Type
            {
                get => string.Empty;
                set => _type = string.Empty;
            }

            protected Param()
            {
                Name = string.Empty;
                Value = string.Empty;
                UnitAccession = string.Empty;
                UnitCVRef = string.Empty;
                UnitName = string.Empty;
                _cvRef = string.Empty;
                _accession = string.Empty;
                _type = string.Empty;
            }
        }

        /// <summary>
        /// Controlled vocabulary parameter
        /// </summary>
        private class CVParam : Param
        {
            /// <summary>
            /// Controlled vocabulary reference
            /// </summary>
            /// <remarks>Required</remarks>
            public override string CVRef
            {
                get => _cvRef;
                set => _cvRef = value;
            }

            /// <summary>
            /// Accession
            /// </summary>
            /// <remarks>Required</remarks>
            public override string Accession
            {
                get => _accession;
                set => _accession = value;
            }

            /// <summary>
            /// Constructor
            /// </summary>
            public CVParam()
            {
                ParamType = ParamType.cvParam;
            }
        }

        /// <summary>
        /// User parameter
        /// </summary>
        private class UserParam : Param
        {
            /// <summary>
            /// Parameter type
            /// </summary>
            /// <remarks>Optional</remarks>
            public override string Type
            {
                get => _type;
                set => _type = value;
            }

            /// <summary>
            /// Constructor
            /// </summary>
            public UserParam()
            {
                ParamType = ParamType.userParam;
            }
        }

        /// <summary>
        /// This class tracks byte offsets for spectra and chromatograms
        /// </summary>
        private class IndexList
        {
            private long _artificialScanNum = 1;
            public class IndexItem // A struct would be faster, but it can also be a pain since it is a value type
            {
                /// <summary>
                /// NativeID
                /// </summary>
                public readonly string Ref;

                /// <summary>
                /// Byte offset
                /// </summary>
                public readonly long Offset;

                /// <summary>
                /// Scan number (1-based) or Chromatogram number (1-based)
                /// </summary>
                /// <remarks>
                /// Even if the first scan in a file has scan number = 200, the value stored here will be 1
                /// </remarks>
                public readonly long IdNum;

                /// <summary>
                /// Constructor
                /// </summary>
                /// <param name="idRef">NativeID, e.g. controllerType=0 controllerNumber=1 scan=1</param>
                /// <param name="offset">Byte offset</param>
                /// <param name="idNum">Scan number (1-based) or Chromatogram number (1-based)</param>
                public IndexItem(string idRef, long offset, long idNum)
                {
                    Ref = idRef;
                    Offset = offset;
                    IdNum = idNum;
                }
            }

            public void Clear()
            {
                _artificialScanNum = 1;
                Offsets.Clear();
                OffsetsMapInt.Clear();
                // OffsetsMapNative.Clear();
                // IdToNativeMap.Clear();
                NativeToIdMap.Clear();
            }

            public IndexListType IndexType = IndexListType.Unknown;
            public enum IndexListType
            {
                Spectrum,
                Chromatogram,
                Unknown,
            }

            /// <summary>
            /// Spectrum offset info, in the order that it was read from the .mzML file
            /// </summary>
            public List<IndexItem> Offsets { get; } = new List<IndexItem>();

            // Unused
            // private readonly Dictionary<string, long> OffsetsMapNative = new Dictionary<string, long>();

            /// <summary>
            /// Keys in this dictionary are artificial scan number
            /// Values are file offsets, in bytes
            /// </summary>
            public readonly Dictionary<long, long> OffsetsMapInt = new Dictionary<long, long>();

            // Unused
            // private readonly Dictionary<long, string> IdToNativeMap = new Dictionary<long, string>();

            /// <summary>
            /// Keys are NativeID
            /// Values are artificial scan number
            /// </summary>
            private readonly Dictionary<string, long> NativeToIdMap = new Dictionary<string, long>();

            // Unused
            // private readonly Dictionary<long, long> ActualScanToIdMap = new Dictionary<long, long>();

            /// <summary>
            /// Store the offset info in the dictionaries
            /// </summary>
            /// <param name="idRef">NativeID</param>
            /// <param name="offset">Byte offset, as text</param>
            public void AddOffset(string idRef, string offset)
            {
                AddOffset(idRef, long.Parse(offset));
            }

            /// <summary>
            /// Store the offset info in the dictionaries
            /// </summary>
            /// <param name="idRef">NativeID</param>
            /// <param name="offset">Byte offset, as a long</param>
            public void AddOffset(string idRef, long offset)
            {
                var scanNum = _artificialScanNum++;

                // Try to convert the NativeId to a scan number
                // This is straightforward for Thermo datasets, but can be problematic for others
                // * Waters .raw files have multiple functions (and spectra) per scan number, unless a filter is used to only select one function
                // * Agilent datasets often have non-contiguous scan numbers (in NativeID)
                // * UIMF datasets have multiple frames, and each frame has its own list of scans
                //
                // if (NativeIdConversion.TryGetScanNumberLong(idRef, out var actualScanNumber))
                // {
                //     if (!ActualScanToIdMap.ContainsKey(actualScanNumber))
                //     {
                //         ActualScanToIdMap.Add(actualScanNumber, scanNum);
                //     }
                // }

                var item = new IndexItem(idRef, offset, scanNum);
                AddMapForOffset(item);
                Offsets.Add(item);
            }

            private void AddMapForOffset(IndexItem item)
            {
                if (IndexType == IndexListType.Chromatogram)
                {
                    // Update the dictionaries, but only if the key is not yet defined
                    AddMapForOffset(item, true);
                }
                else
                {
                    // Update the dictionaries
                    // If a duplicate key exists an exception will be thrown
                    AddMapForOffset(item, false);
                }
            }

            private void AddMapForOffset(IndexItem item, bool skipDuplicateKeys)
            {
                // if (!skipDuplicateKeys || !OffsetsMapNative.ContainsKey(item.Ref))
                // {
                //     OffsetsMapNative.Add(item.Ref, item.Offset);
                // }

                // This won't be sufficient until there is a valid parser for all forms of NativeID.
                // Using artificial scan number for now.

                OffsetsMapInt.Add(item.IdNum, item.Offset);
                // IdToNativeMap.Add(item.IdNum, item.Ref);

                if (!skipDuplicateKeys || !NativeToIdMap.ContainsKey(item.Ref))
                {
                    NativeToIdMap.Add(item.Ref, item.IdNum);
                }

                /*if (IndexType == IndexListType.Spectrum)
                {
                    long id = Int64.Parse(item.Ref.Substring(item.Ref.LastIndexOfAny(new char[] {'=', 'F'}) + 1));
                    OffsetsMapInt.Add(id, item.Offset);
                    IdNativeMap.Add(id, item.Ref);
                    item.IdNum = id;
                }*/
            }

            public void RegenerateMaps()
            {
                // OffsetsMapNative.Clear();
                OffsetsMapInt.Clear();
                // IdToNativeMap.Clear();
                NativeToIdMap.Clear();
                foreach (var offsetItem in Offsets)
                {
                    AddMapForOffset(offsetItem);
                }
            }
        }

        /// <summary>
        /// Helper class for converting between native IDs and scan numbers
        /// </summary>
        public static class NativeIdConversion
        {
            private static Dictionary<string, string> ParseNativeId(string nativeId)
            {
                var tokens = nativeId.Split(new[] { '\t', ' ', '\n', '\r' }, StringSplitOptions.RemoveEmptyEntries);
                var map = new Dictionary<string, string>();
                foreach (var token in tokens)
                {
                    var equals = token.IndexOf('=');
                    var name = token.Substring(0, equals);
                    var value = token.Substring(equals + 1);
                    map.Add(name, value);
                }
                return map;
            }

            /// <summary>
            /// Try to get the scan number for <paramref name="nativeId"/>, and parse it as a long int
            /// </summary>
            /// <param name="nativeId"></param>
            /// <param name="num"></param>
            /// <returns>True if a numeric scan number was found, otherwise false</returns>
            public static bool TryGetScanNumberLong(string nativeId, out long num)
            {
                return long.TryParse(GetScanNumber(nativeId), out num);
            }

            /// <summary>
            /// Try to get the scan number for <paramref name="nativeId"/>, and parse it as an int
            /// </summary>
            /// <param name="nativeId"></param>
            /// <param name="num"></param>
            /// <returns>True if a numeric scan number was found, otherwise false</returns>
            public static bool TryGetScanNumberInt(string nativeId, out int num)
            {
                return int.TryParse(GetScanNumber(nativeId), out num);
            }

            /// <summary>
            /// For the supplied <paramref name="nativeId"/>, get the corresponding scan number
            /// </summary>
            /// <param name="nativeId"></param>
            /// <returns>Scan number (as a string)</returns>
            /// <remarks>
            /// <para>
            /// If the scan number cannot be extracted, returns the original nativeId
            /// </para>
            /// <para>
            /// Code is ported from MSData.cpp in ProteoWizard
            /// </para>
            /// </remarks>
            public static string GetScanNumber(string nativeId)
            {
                // TODO: Add interpreter for Waters' S0F1, S1F1, S0F2,... format
                //switch (nativeIdFormat)
                //{
                //    case MS_spectrum_identifier_nativeID_format: // mzData
                //        return value(id, "spectrum");
                //
                //    case MS_multiple_peak_list_nativeID_format: // MGF
                //        return value(id, "index");
                //
                //    case MS_Agilent_MassHunter_nativeID_format:
                //        return value(id, "scanId");
                //
                //    case MS_Thermo_nativeID_format:
                //        // conversion from Thermo nativeIDs assumes default controller information
                //        if (id.find("controllerType=0 controllerNumber=1") != 0)
                //            return string.Empty;
                //
                //        // fall through to get scan
                //
                //    case MS_Bruker_Agilent_YEP_nativeID_format:
                //    case MS_Bruker_BAF_nativeID_format:
                //    case MS_scan_number_only_nativeID_format:
                //        return value(id, "scan");
                //
                //    default:
                //        if (bal::starts_with(id, "scan=")) return value(id, "scan");
                //        else if (bal::starts_with(id, "index=")) return value(id, "index");
                //        return string.Empty;
                //}
                if (nativeId.Contains("="))
                {
                    var map = ParseNativeId(nativeId);
                    if (map.ContainsKey("spectrum"))
                    {
                        return map["spectrum"];
                    }
                    if (map.ContainsKey("index"))
                    {
                        return map["index"];
                    }
                    if (map.ContainsKey("scanId"))
                    {
                        return map["scanId"];
                    }
                    if (map.ContainsKey("scan"))
                    {
                        return map["scan"];
                    }
                }

                // No equals sign, don't have parser breakdown
                // Or key data not found in breakdown of nativeId
                return nativeId;
            }

            //public static string GetNativeId(string scanNumber)
            //{
            //    switch (nativeIdFormat)
            //    {
            //        case MS_Thermo_nativeID_format:
            //            return "controllerType=0 controllerNumber=1 scan=" + scanNumber;
            //
            //        case MS_spectrum_identifier_nativeID_format:
            //            return "spectrum=" + scanNumber;
            //
            //        case MS_multiple_peak_list_nativeID_format:
            //            return "index=" + scanNumber;
            //
            //        case MS_Agilent_MassHunter_nativeID_format:
            //            return "scanId=" + scanNumber;
            //
            //        case MS_Bruker_Agilent_YEP_nativeID_format:
            //        case MS_Bruker_BAF_nativeID_format:
            //        case MS_scan_number_only_nativeID_format:
            //            return "scan=" + scanNumber;
            //
            //        default:
            //            return string.Empty;
            //    }
            //}
        }

        private readonly Dictionary<string, List<Param>> _referenceableParamGroups = new Dictionary<string, List<Param>>();

        /// <summary>
        /// Selected Ion m/z
        /// </summary>
        private class SelectedIon
        {
            /// <summary>
            /// Selected Ion m/z
            /// </summary>
            public double SelectedIonMz;

            /// <summary>
            /// Selected Ion charge
            /// </summary>
            public int Charge;

            // ReSharper disable once NotAccessedField.Local
            public int OldCharge;

            /// <summary>
            /// Constructor
            /// </summary>
            public SelectedIon()
            {
                SelectedIonMz = 0.0;
                Charge = 0;
                OldCharge = 0;
            }
        }

        private class Precursor
        {
            public List<SelectedIon> Ions;
            public double IsolationWindowTargetMz;
            public double IsolationWindowLowerOffset;
            public double IsolationWindowUpperOffset;
            public ActivationMethod Activation;

            public Precursor()
            {
                Ions = new List<SelectedIon>();
                IsolationWindowTargetMz = 0.0;
                IsolationWindowLowerOffset = 0.0;
                IsolationWindowUpperOffset = 0.0;
                Activation = ActivationMethod.Unknown;
            }
        }

        private enum Precision
        {
            Precision32,
            Precision64
        };

        private enum ArrayType
        {
            m_z_array,
            intensity_array,
            charge_array,
            signal_to_noise_array,
            time_array,
            wavelength_array,
            non_standard_data_array,
            flow_rate_array,
            pressure_array,
            temperature_array
        };

        private enum Instrument
        {
            ABI_WIFF_format, //MS_ABI_WIFF_format, "MS:1000562", "ABI WIFF format", "Applied Biosystems WIFF file format."
            Thermo_RAW_format, //MS_Thermo_RAW_format, "MS:1000563", "Thermo RAW format", "Thermo Scientific RAW file format."
            Waters_raw_format, //MS_Waters_raw_format, "MS:1000526", "Waters raw format", "Waters data file format found in a Waters RAW directory, generated from an MS acquisition."
            Unknown
        };

        //private Instrument _instrument;

        private class BinaryDataArray
        {
            public int ArrayLength;
            public Precision Precision;
            public ArrayType ArrayType;
            public double[] Data;

            public BinaryDataArray()
            {
                Data = null;
                Precision = Precision.Precision32;
                ArrayType = ArrayType.m_z_array;
                ArrayLength = 0;
            }
        }

        #endregion

        #region Constructor

        /// <summary>
        /// Initialize a MzMLReader object
        /// </summary>
        /// <param name="filePath">Path to mzML file</param>
        /// <param name="randomAccess">If mzML reader should be configured for random access</param>
        /// <param name="tryReducingMemoryUsage">
        /// Set to true if the mzML reader should try to avoid reading all spectra into memory.
        /// This will reduce memory usage for a non-random access MzMLReader, as long as ReadMassSpectrum(int) isn't used.
        /// </param>
        public MzMLReader(string filePath, bool randomAccess = false, bool tryReducingMemoryUsage = true)
        {
            _filePath = filePath;
            //_instrument = Instrument.Unknown;
            _version = MzML_Version.mzML1_1_0;
            _randomAccess = randomAccess;
            _reduceMemoryUsage = tryReducingMemoryUsage;
            _unzippedFilePath = _filePath;

            ConfigureFileHandles();
        }

        private void ConfigureFileHandles()
        {
            _file?.Dispose();

            var sourceFile = new FileInfo(_filePath);
            if (!sourceFile.Exists)
                throw new FileNotFoundException(".mzML file not found", _filePath);

            // Set a very large read buffer, it does decrease the read times for uncompressed files.
            _file = new FileStream(sourceFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite, 65536);

            /*****************************************************************************************************************************************************
             * TODO: Change how the file handles are used for safety purposes - open up each time, or what?
             *****************************************************************************************************************************************************/

            if (sourceFile.Name.Trim().EndsWith(".mzML.gz", StringComparison.OrdinalIgnoreCase))
            {
                _isGzipped = true;
                var zipStreamFile = new GZipStream(_file, CompressionMode.Decompress);
                if (!_randomAccess)
                {
                    _file = zipStreamFile;
                }
                else
                {
                    // Unzip the file to the temp path
                    _unzippedFilePath = Path.Combine(Path.GetTempPath(), Path.GetFileNameWithoutExtension(sourceFile.Name));
                    using (_file)
                    using (zipStreamFile)
                    using (
                        var tempFile = new FileStream(_unzippedFilePath, FileMode.Create, FileAccess.ReadWrite,
                            FileShare.None, 65536))
                    {
                        zipStreamFile.CopyTo(tempFile/*, 65536*/);
                    }
                    _file = new FileStream(_unzippedFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite, 65536);
                }
            }
            _fileReader = new StreamReader(_file, Encoding.UTF8, true, 65536);

            if (!_isGzipped || _randomAccess) // can't reset the position on a gzipped file...
            {
                // perform a read to perform encoding auto detection
                _fileReader.ReadLine();
                _encoding = _fileReader.CurrentEncoding;
                // Reset to beginning of file.
                _fileReader.DiscardBufferedData();
                _fileReader.BaseStream.Position = 0;
            }
        }

        #endregion

        #region Public interface functions

        /// <summary>
        /// The number of spectra in the file
        /// </summary>
        public int NumSpectra
        {
            get
            {
                RequireMetadata();
                return (int)_numSpectra;
            }
        }

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes
        /// </summary>
        /// <remarks>
        /// Child term of PSI-MS term MS:1000767, native spectrum identifier format
        /// </remarks>
        public CV.CVID NativeIdFormat
        {
            get
            {
                RequireMetadata();
                return _nativeIdFormat;
            }
        }

        /// <summary>
        /// The Native Format of the source file - needed for tracking purposes
        /// </summary>
        /// <remarks>
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </remarks>
        public CV.CVID NativeFormat
        {
            get
            {
                RequireMetadata();
                return _nativeFormat;
            }
        }

        /// <summary>
        /// Re-opens the file using random access
        /// </summary>
        /// <returns>true if successful, exception if an error</returns>
        public bool TryMakeRandomAccessCapable()
        {
            _randomAccess = true;
            ConfigureFileHandles(); // Reopen the files
            return true;
        }

        /// <summary>
        /// Path to the file; is <see cref="string.Empty"/> if the reader is in-memory
        /// </summary>
        // ReSharper disable once ConvertToAutoPropertyWithPrivateSetter
        public string FilePath => _unzippedFilePath;

        /// <summary>
        /// SHA-1 Checksum of the original input file (raw, mzML, .d folder, etc.)
        /// </summary>
        public string SrcFileChecksum
        {
            get
            {
                RequireMetadata();
                return _srcFileChecksum;
            }
        }

        /// <summary>
        /// Version of the immediate prior input file (raw, mzML, .d folder, etc.)
        /// </summary>
        public string FileFormatVersion
        {
            get
            {
                RequireMetadata();
                return _fileFormatVersion;
            }
        }

        /// <summary>
        /// Read the file-level metadata from the mzML file, without reading any spectra
        /// </summary>
        private void RequireMetadata()
        {
            if (!_haveMetaData)
            {
                var tempBool = _reduceMemoryUsage; // Set a flag to avoid reading the entire file before returning.
                ReadMzMl(); // Read the index and metadata so that the offsets get populated
                            // The number of spectra is an attribute in the spectrumList tag
                _reduceMemoryUsage = tempBool;
            }
        }

        /// <summary>
        /// Returns all mass spectra.
        /// ReadAllSpectraNonRandom and ReadAllSpectraRandom use "yield return" to allow processing one spectra at a time if called from a for each loop statement.
        /// </summary>
        /// <param name="includePeaks">true to include peak data; only used for random-access capable reading</param>
        /// <returns>Enumerable list of spectrum objects</returns>
        public IEnumerable<Spectrum> ReadAllSpectra(bool includePeaks = true)
        {
            if (!_randomAccess)
            {
                return ReadAllSpectraNonRandom();
            }
            return ReadAllSpectraRandom(includePeaks);
        }

        /// <summary>
        /// <para>
        /// Returns a single spectrum from the file
        /// </para>
        /// <para>
        /// If Random access is disabled, this method calls method ReadMassSpectrumNonRandom
        /// That will cause all spectra in the file to be loaded into memory
        /// </para>
        /// </summary>
        /// <param name="index">
        /// 1-based index of the spectrum in the file
        /// Using index = 1 will return the first spectrum in the file, regardless of its actual scan number
        /// </param>
        /// <param name="includePeaks">true to include peak data (ignored if _randomAccess is false)</param>
        /// <returns>Mass spectrum</returns>
        /// <remarks>
        /// <para>
        /// To enable random access, either set randomAccess to true when instantiating the MzMLReader class, or call TryMakeRandomAccessCapable
        /// </para>
        /// <para>
        /// If random access mode is turned on, this will respond quickly and use only as much memory as is needed to store the spectrum.
        /// If random access mode is off, this will cause the memory usage reducing mode to shut off, and all spectra will be read into memory.
        /// </para>
        /// </remarks>
        public Spectrum ReadMassSpectrum(int index, bool includePeaks = true)
        {
            if (!_randomAccess)
            {
                // Proper functionality when not random access
                return ReadMassSpectrumNonRandom(index);
            }
            return ReadMassSpectrumRandom(index, includePeaks);
        }

        /// <summary>
        /// Read the specified spectrum from the file, optionally reading only the metadata
        /// </summary>
        /// <param name="scanNum">
        /// 1-based index of the spectrum in the file
        /// Using index = 1 will return the first spectrum in the file, regardless of its actual scan number
        /// </param>
        /// <param name="includePeaks"></param>
        /// <returns>Spectrum object</returns>
        public Spectrum GetSpectrum(int scanNum, bool includePeaks = true)
        {
            return ReadMassSpectrum(scanNum, includePeaks);
        }

        #endregion

        #region Interface Helper Functions: Non-Random Access

        /// <summary>
        /// Read all mass spectra in the file, not using random access
        /// Uses "yield return" to use less memory when called from a "foreach" statement
        /// </summary>
        /// <returns>Enumerable list of spectrum objects</returns>
        private IEnumerable<Spectrum> ReadAllSpectraNonRandom()
        {
            if (_reduceMemoryUsage)
            {
                //_artificialScanNum = 1;
                if (!_haveMetaData)
                {
                    ReadMzMl();
                }

                while (_xmlReaderForYield.ReadState == ReadState.Interactive)
                {
                    // Handle exiting out properly at EndElement tags
                    if (_xmlReaderForYield.NodeType != XmlNodeType.Element)
                    {
                        _xmlReaderForYield.Read();
                        continue;
                    }
                    if (_xmlReaderForYield.Name == "spectrum")
                    {
                        // Schema requirements: zero to many instances of this element
                        // Use reader.ReadSubtree() to provide an XmlReader that is only valid for the element and child nodes
                        yield return ReadSpectrum(_xmlReaderForYield.ReadSubtree());
                        // "spectrum" might not have any child nodes
                        // We will either consume the EndElement, or the same element that was passed to ReadSpectrum (in case of no child nodes)
                        _xmlReaderForYield.Read();
                    }
                    else
                    {
                        _xmlReaderForYield.Skip();
                    }
                }
            }
            else
            {
                if (!_allRead)
                {
                    ReadMzMl();
                }
                foreach (var spec in _spectra)
                {
                    yield return spec;
                }
            }
        }

        /// <summary>
        /// Read a single mass spectrum and return it
        /// </summary>
        /// <param name="index">
        /// 1-based index of the spectrum in the file
        /// Using index = 1 will return the first spectrum in the file, regardless of its actual scan number
        /// </param>
        /// <remarks>Causes all spectra in the file to be loaded into memory</remarks>
        /// <returns>Spectrum object</returns>
        private Spectrum ReadMassSpectrumNonRandom(int index)
        {
            if (!_allRead)
            {
                //_artificialScanNum = 1;
                _reduceMemoryUsage = false; // They called this on a non-random access reader, now they suffer the consequences.
                ReadMzMl();
            }

            if (index == 0)
            {
                Console.WriteLine("Warning: the index argument for method ReadMassSpectrumNonRandom is 1-based");
                return null;
            }

            return _spectra[index - 1];
        }

        #endregion

        #region Interface Helper Functions: Random Access

        /// <summary>
        /// Read all mass spectra in the file, using random access
        /// Uses "yield return" to use less memory when called from a "foreach" statement
        /// </summary>
        /// <param name="includePeaks">true to include peak data</param>
        /// <returns>Enumerable list of spectrum objects</returns>
        private IEnumerable<Spectrum> ReadAllSpectraRandom(bool includePeaks = true)
        {
            if (!_haveIndex || !_haveMetaData)
            {
                ReadMzMl(); // Read the index and metadata so that the offsets get populated.
            }
            foreach (var specIndex in _spectrumOffsets.Offsets)
            {
                yield return ReadMassSpectrumRandom(specIndex.IdNum, includePeaks);
            }
        }

        /// <summary>
        /// Read a single mass spectrum and return it
        /// </summary>
        /// <param name="index">
        /// 1-based index of the spectrum in the file
        /// Using index = 1 will return the first spectrum in the file, regardless of its actual scan number
        /// </param>
        /// <param name="includePeaks"></param>
        private Spectrum ReadMassSpectrumRandom(long index, bool includePeaks = true)
        {
            if (!_haveIndex || !_haveMetaData)
            {
                ReadMzMl();
            }

            if (!_spectrumOffsets.OffsetsMapInt.ContainsKey(index))
            {
                return null;
            }

            _fileReader.DiscardBufferedData();
            _fileReader.BaseStream.Position = _spectrumOffsets.OffsetsMapInt[index];
            // Not allowed for a GZipStream.....
            using (var reader = XmlReader.Create(_fileReader, _xSettings))
            {
                reader.MoveToContent();
                return ReadSpectrum(reader.ReadSubtree(), includePeaks);
            }
        }

        #endregion

        #region Cleanup functions

        /// <summary>
        /// Close out the file handle and delete any temp files
        /// </summary>
        public void Close()
        {
            _xmlReaderForYield?.Dispose();
            _fileReader?.Dispose();
            _file?.Dispose();
        }

        /// <summary>
        /// Delete unzipped file, if we had to unzip the file to read it
        /// </summary>
        public void Cleanup()
        {
            if (_randomAccess && _isGzipped && File.Exists(_unzippedFilePath))
            {
                try
                {
                    File.Delete(_unzippedFilePath);
                }
                // ReSharper disable once EmptyGeneralCatchClause
                catch
                { } // We failed, it's not optimal, but it's okay.
            }
        }

        /// <summary>
        /// Close the file handles and cleanup any temp files
        /// </summary>
        public void Dispose()
        {
            Dispose(true);
        }

        /// <summary>
        /// Close the file handles and cleanup any temp files
        /// </summary>
        ~MzMLReader()
        {
            Dispose(false);
        }

        /// <summary>
        /// Dispose, with a decreased chance of issues/failure
        /// </summary>
        /// <param name="disposing"></param>
        private void Dispose(bool disposing)
        {
            if (disposing)
            {
                GC.SuppressFinalize(this);
            }
            Close();
            Cleanup();
        }

        /// <summary>
        /// Clear out cached data - keep the index information, if it is a random access reader
        /// </summary>
        public void ClearDataCache()
        {
            _spectra.Clear();
            _allRead = false;
        }

        #endregion

        #region Index reading functions

        /// <summary>
        /// Find and read the index information, starting at the end of the file
        /// </summary>
        private void ReadIndexFromEnd()
        {
            var stream = new FileStream(_unzippedFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite, 1);
            var testPos = stream.Length;

            var streamReader = new StreamReader(stream, Encoding.UTF8, true, 65536);
            streamReader.DiscardBufferedData();
            var haveOffset = false;

            while (!haveOffset)
            {
                // "<indexListOffset"
                const int bufSize = 512; //65536 (17 bits), 131072 (18 bits), 262144 (19 bits), 524288 (20 bits)
                testPos -= bufSize;
                stream.Position = testPos;
                var byteBuffer = new byte[bufSize];
                while (stream.Position < stream.Length && !haveOffset)
                {
                    var bufStart = stream.Position;
                    var bytesRead = stream.Read(byteBuffer, 0, bufSize);
                    var stringBuffer = _encoding.GetString(byteBuffer, 0, bytesRead);
                    // set up the rewind to ensure full tags
                    var lastTagEnd = stringBuffer.LastIndexOf('>');
                    var lastTagStart = stringBuffer.LastIndexOf('<');
                    if (lastTagStart != -1 && lastTagEnd != -1 && lastTagStart > lastTagEnd)
                    {
                        var endOfString = lastTagStart;
                        var rewindBy = _encoding.GetByteCount(stringBuffer.Substring(endOfString));
                        stringBuffer = stringBuffer.Substring(0, endOfString);
                        stream.Seek(-rewindBy, SeekOrigin.Current);
                    }

                    var found = stringBuffer.IndexOf("<indexListOffset", StringComparison.OrdinalIgnoreCase);
                    if (found >= 0)
                    {
                        var pos = bufStart + _encoding.GetByteCount(stringBuffer.Substring(0, found));
                        streamReader.DiscardBufferedData();
                        streamReader.BaseStream.Position = pos;
                        using (var reader = XmlReader.Create(streamReader, _xSettings))
                        {
                            reader.MoveToContent();
                            var reader2 = reader.ReadSubtree(); // Get past root element problems
                            reader2.MoveToContent();
                            _indexListOffset = reader2.ReadElementContentAsLong();
                            reader2.Dispose();
                        }
                        haveOffset = true;
                    }
                }
            }
            if (_indexListOffset < stream.Length / 2) // Probably invalid, now we must search...
            {
                // "<indexList"
                haveOffset = false;
                streamReader.DiscardBufferedData();
                const int bufSize = 524588; //65536 (17 bits), 131072 (18 bits), 262144 (19 bits), 524288 (20 bits)
                testPos = stream.Length;
                while (!haveOffset)
                {
                    testPos -= bufSize;
                    stream.Position = testPos;
                    var byteBuffer = new byte[bufSize];
                    while (stream.Position < stream.Length && !haveOffset)
                    {
                        var bufStart = stream.Position;
                        var bytesRead = stream.Read(byteBuffer, 0, bufSize);
                        var stringBuffer = _encoding.GetString(byteBuffer, 0, bytesRead);
                        // set up the rewind to ensure full tags
                        var lastTagEnd = stringBuffer.LastIndexOf('>');
                        var lastTagStart = stringBuffer.LastIndexOf('<');
                        if (lastTagStart != -1 && lastTagEnd != -1 && lastTagStart > lastTagEnd)
                        {
                            var endOfString = lastTagStart;
                            var rewindBy = _encoding.GetByteCount(stringBuffer.Substring(endOfString));
                            stringBuffer = stringBuffer.Substring(0, endOfString);
                            stream.Seek(-rewindBy, SeekOrigin.Current);
                        }

                        var found = stringBuffer.IndexOf("<indexList ", StringComparison.OrdinalIgnoreCase);
                        if (found >= 0)
                        {
                            var pos = bufStart + _encoding.GetByteCount(stringBuffer.Substring(0, found));
                            _indexListOffset = pos;
                            haveOffset = true;
                        }
                    }
                }
            }

            // Now we definitely have the offset of the indexOffsetList... (unless the file is invalid)
            // Create the XmlReader at the right position, and read.
            streamReader.DiscardBufferedData();
            streamReader.BaseStream.Position = _indexListOffset;
            using (var reader = XmlReader.Create(streamReader, _xSettings))
            {
                reader.MoveToContent();
                ReadIndexList(reader.ReadSubtree());
            }
            var isValid = true;
            // Validate the index - if there are duplicate offsets, it is probably invalid
            var collisions = new Dictionary<long, int>();
            foreach (var index in _spectrumOffsets.Offsets)
            {
                if (!collisions.ContainsKey(index.Offset))
                {
                    collisions.Add(index.Offset, 0);
                }
                else
                {
                    isValid = false;
                    collisions[index.Offset]++;
                }
            }
            _haveIndex = isValid;
        }

        /// <summary>
        /// Read the Checksum from the indexed mzML data
        /// </summary>
        private void ReadChecksum()
        {
            using (var stream = new FileStream(_unzippedFilePath, FileMode.Open, FileAccess.Read, FileShare.Read, 1))
            using (var streamReader = new StreamReader(stream, Encoding.UTF8, true, 65536))
            {
                try
                {
                    stream.Position = stream.Length - 500;
                    streamReader.DiscardBufferedData();
                    var data = streamReader.ReadToEnd();
                    var pos = data.IndexOf("<fileChecksum", StringComparison.OrdinalIgnoreCase);
                    if (pos >= 0)
                    {
                        data = data.Substring(pos);
                        pos = data.IndexOf('>') + 1;
                        data = data.Substring(pos);
                        pos = data.IndexOf('<');
                        data = data.Substring(0, pos);
                        if (data.Length == 40)
                        {
                            _srcFileChecksum = data;
                        }
                    }
                }
                catch
                {
                    // Dropping errors - if this doesn't work, we'll just checksum the whole file.
                }
            }
        }

        /// <summary>
        /// Handle the child nodes of the run element
        /// Called by IndexMzMl (XML hierarchy)
        /// </summary>
        private void ReadRunForOffsets()
        {
            // Set the buffer to 1 byte (minimum allowed value), since we need accurate positions for the indices
            Stream file = new FileStream(_unzippedFilePath, FileMode.Open, FileAccess.Read, FileShare.Read, 1);

            _chromatogramOffsets.Clear();
            _spectrumOffsets.Clear();

            using (file)
            {
                const string specTag = "spectrum";
                const string chromTag = "chromatogram";
                const int maxRead = 524288; //65536 (17 bits), 131072 (18 bits), 262144 (19 bits), 524288 (20 bits)
                var byteBuffer = new byte[maxRead];
                while (file.Position < file.Length)
                {
                    var bufStart = file.Position;
                    var bytesRead = file.Read(byteBuffer, 0, maxRead);
                    var stringBuffer = _encoding.GetString(byteBuffer, 0, bytesRead);
                    // set up the rewind to ensure full tags
                    var lastTagEnd = stringBuffer.LastIndexOf('>');
                    var lastTagStart = stringBuffer.LastIndexOf('<');
                    if (lastTagStart != -1 && lastTagEnd != -1 && lastTagStart > lastTagEnd)
                    {
                        var endOfString = lastTagStart;
                        var rewindBy = _encoding.GetByteCount(stringBuffer.Substring(endOfString));
                        stringBuffer = stringBuffer.Substring(0, endOfString);
                        file.Seek(-rewindBy, SeekOrigin.Current);
                    }

                    var searchPoint = 0;
                    while (searchPoint < stringBuffer.Length)
                    {
                        var foundSpec = stringBuffer.IndexOf("<" + specTag + " ", searchPoint, StringComparison.Ordinal);
                        var foundChrom = stringBuffer.IndexOf("<" + chromTag + " ", searchPoint, StringComparison.Ordinal);

                        if (foundSpec >= 0)
                        {
                            searchPoint = foundSpec;
                        }
                        else if (foundChrom >= 0)
                        {
                            searchPoint = foundChrom;
                        }
                        else
                        {
                            break;
                        }

                        var pos = bufStart + _encoding.GetByteCount(stringBuffer.Substring(0, searchPoint));
                        var end = stringBuffer.IndexOf('>', searchPoint + 1);
                        // Grab everything between '<' and the next '>'
                        var builder = stringBuffer.Substring(searchPoint + 1, end - 1 - (searchPoint + 1));
                        // Get the ID of the tag
                        var attributeName = "id";
                        if (_version == MzML_Version.mzML1_0_0)
                        {
                            attributeName = "nativeID";
                        }
                        var idIndex = builder.IndexOf(attributeName + "=\"", StringComparison.OrdinalIgnoreCase);
                        var idOpenQuote = builder.IndexOf("\"", idIndex, StringComparison.OrdinalIgnoreCase);
                        var idCloseQuote = builder.IndexOf("\"", idOpenQuote + 1, StringComparison.OrdinalIgnoreCase);
                        var length = idCloseQuote - idOpenQuote - 1;
                        var id = builder.Substring(idOpenQuote + 1, length);
                        // Add offset to the correct list
                        if (builder.StartsWith(specTag))
                        {
                            _spectrumOffsets.AddOffset(id, pos);
                        }
                        else if (builder.StartsWith(chromTag))
                        {
                            _chromatogramOffsets.AddOffset(id, pos);
                        }

                        // Force find of the next tag
                        searchPoint = end;
                    }
                }
                // Read through the entire file, searching for tags that start with "<s" or with "<c"
                /*while (file.Position < file.Length)
                {
                    while (file.Position < file.Length)
                    {
                        if (file.ReadByte() == '<')
                        {
                            position = file.Position - 1; // position of caret
                            int spaceCaretSlash = 0;
                            for (int i = 0; i < 13; i++) // 13: "chromatogram "
                            {
                                char c = (char)(file.ReadByte());
                                builder += c;
                                if (" >/".IndexOf(c) >= 0)
                                {
                                    spaceCaretSlash = i;
                                    break;
                                }
                            }
                            if (builder[0] == 's' || builder[0] == 'c')
                            {
                                if (spaceCaretSlash == specTag.Length || spaceCaretSlash == chromTag.Length)
                                {
                                    string tagName = builder.Substring(0, spaceCaretSlash);
                                    if (string.Equals(tagName.ToLower(), specTag) ||
                                        string.Equals(tagName.ToLower(), chromTag))
                                    {
                                        break;
                                    }
                                }
                            }
                            // reset to empty if we didn't break out
                            builder = string.Empty;
                        }
                    }
                    // We have a '<', followed by 's' or 'c'. Store position, and check the tag name.
                    // searching for "spectrum" and "chromatogram"
                    // We have "spectrum" or "chromatogram"
                    // Assemble the tag and attributes
                    while (file.Position < file.Length)
                    {
                        char c = (char)(file.ReadByte());
                        if (c == '>')
                        {
                            break;
                        }
                        builder += c;
                    }
                    // Get the ID of the tag
                    var idIndex = builder.IndexOf("id=\"");
                    var idOpenQuote = builder.IndexOf("\"", idIndex);
                    var idCloseQuote = builder.IndexOf("\"", idOpenQuote + 1);
                    var length = idCloseQuote - idOpenQuote - 1;
                    var id = builder.Substring(idOpenQuote + 1, length);
                    // Add offset to the correct list
                    if (builder.StartsWith(specTag))
                    {
                        _spectrumOffsets.AddOffset(id, position);
                    }
                    else if (builder.StartsWith(chromTag))
                    {
                        _chromatogramOffsets.AddOffset(id, position);
                    }
                }*/
                _haveIndex = true;
            }
        }

        /// <summary>
        /// Handle the child nodes of the indexed mzML element
        /// Called by IndexMzMl (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "indexList" element</param>
        private void ReadIndexList(XmlReader reader)
        {
            reader.MoveToContent();
            reader.ReadStartElement("indexList"); // Throws exception if we are not at the "run" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                switch (reader.Name)
                {
                    case "index":
                        // Schema requirements: zero to one instances of this element
                        // Use reader.ReadSubtree() to provide an XmlReader that is only valid for the element and child nodes
                        ReadIndex(reader.ReadSubtree());
                        // "spectrumList" might not have any child nodes
                        // We will either consume the EndElement, or the same element that was passed to ReadSpectrumList (in case of no child nodes)
                        reader.Read();
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            _haveIndex = true;
            reader.Dispose();
        }

        /// <summary>
        /// Handle the child nodes of the indexList element
        /// Called by ReadIndexList (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "index" element</param>
        private void ReadIndex(XmlReader reader)
        {
            reader.MoveToContent();
            var iType = reader.GetAttribute("name");
            var eType = IndexList.IndexListType.Unknown;
            if (iType != null && string.Equals(iType, "spectrum", StringComparison.OrdinalIgnoreCase))
            {
                eType = IndexList.IndexListType.Spectrum;
            }
            else if (iType != null && string.Equals(iType, "chromatogram", StringComparison.OrdinalIgnoreCase))
            {
                eType = IndexList.IndexListType.Chromatogram;
            }
            reader.ReadStartElement("index"); // Throws exception if we are not at the "run" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                switch (reader.Name)
                {
                    case "offset":
                        // Schema requirements: zero to one instances of this element
                        // Use reader.ReadSubtree() to provide an XmlReader that is only valid for the element and child nodes
                        var idRef = reader.GetAttribute("idRef");
                        var offset = reader.ReadElementContentAsString(); // Reads the start element, content, and end element
                        switch (eType)
                        {
                            case IndexList.IndexListType.Spectrum:
                                _spectrumOffsets.AddOffset(idRef, offset);
                                break;
                            case IndexList.IndexListType.Chromatogram:
                                _chromatogramOffsets.AddOffset(idRef, offset);
                                break;
                        }
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();
        }

        #endregion

        #region Root tag reader

        /// <summary>
        /// Read and parse a .mzML file
        /// Files are commonly larger than 100 MB, so use a streaming reader instead of a DOM reader
        /// Very conditional, depending on configuration
        /// </summary>
        private void ReadMzMl()
        {
            if (_randomAccess && _haveIndex && _haveMetaData)
            {
                return;
            }
            // Handle disposal of allocated object correctly

            if (!_isGzipped || _randomAccess) // can't reset the position on a gzipped file...
            {
                if (_fileReader.BaseStream.Position > 0)
                {
                    // Reset to beginning of file.
                    _fileReader.DiscardBufferedData();
                    _fileReader.BaseStream.Position = 0;
                }
            }

            var reader = XmlReader.Create(_fileReader, _xSettings);
            // Guarantee a move to the root node
            reader.MoveToContent();
            if (_encoding == null)
            {
                _encoding = _fileReader.CurrentEncoding;
            }
            XmlReader indexReader = null;
            if (reader.Name == "indexedmzML")
            {
                indexReader = reader;
                // Read to the mzML root tag, and ignore the extra indexed mzML data
                reader.ReadToDescendant("mzML");
                if (_randomAccess && !_haveIndex)
                {
                    // run to the end of the file (using stream.position = stream.length) and jump backwards to read the index first, and then read the file for needed data
                    ReadIndexFromEnd();
                }
                ReadChecksum();
                reader = reader.ReadSubtree();
                reader.MoveToContent();
            }
            if (string.IsNullOrWhiteSpace(_srcFileChecksum))
            {
                using (var fs = new FileStream(_unzippedFilePath, FileMode.Open, FileAccess.Read, FileShare.Read))
                using (var sha1 = new SHA1Managed())
                {
                    var hash = sha1.ComputeHash(fs);
                    _srcFileChecksum = BitConverter.ToString(hash).ToLower().Replace("-", string.Empty);
                }
            }
            var schemaName = reader.GetAttribute("xsi:schemaLocation");
            // We automatically assume it uses the mzML_1.1.0 schema. Check for the old version.
            //if (!schemaName.Contains("mzML1.1.0.xsd"))
            if (schemaName?.Contains("mzML1.0.0.xsd") == true)
            {
                _version = MzML_Version.mzML1_0_0;
            }
            _fileFormatVersion = reader.GetAttribute("version");
            // Consume the mzML root tag
            // Throws exception if we are not at the "mzML" tag.
            // This is a critical error; we want to stop processing for this file if we encounter this error
            reader.ReadStartElement("mzML");
            var continueReading = true;
            // Read the next node - should be the first child node
            while (reader.ReadState == ReadState.Interactive && continueReading)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                // Handle each 1st level as a chunk
                switch (reader.Name)
                {
                    case "cvList":
                        // Schema requirements: one instance of this element
                        reader.Skip();
                        break;
                    case "fileDescription":
                        // Schema requirements: one instance of this element
                        if (!_randomAccess || (_randomAccess && !_haveMetaData))
                        {
                            ReadFileDescription(reader.ReadSubtree());
                            reader.ReadEndElement(); // "fileDescription" must have child nodes
                        }
                        else
                        {
                            reader.Skip();
                        }
                        break;
                    case "referenceableParamGroupList":
                        // Schema requirements: zero to one instances of this element
                        if (!_randomAccess || (_randomAccess && !_haveMetaData))
                        {
                            ReadReferenceableParamGroupList(reader.ReadSubtree());
                            reader.ReadEndElement(); // "referenceableParamGroupList" must have child nodes
                        }
                        else
                        {
                            reader.Skip();
                        }
                        break;
                    case "sampleList":
                        // Schema requirements: zero to one instances of this element
                        reader.Skip();
                        break;
                    case "softwareList":
                        // Schema requirements: one instance of this element
                        reader.Skip();
                        break;
                    case "scanSettingsList":
                        // Schema requirements: zero to one instances of this element
                        reader.Skip();
                        break;
                    case "instrumentConfigurationList":
                        // Schema requirements: one instance of this element
                        reader.Skip();
                        break;
                    case "dataProcessingList":
                        // Schema requirements: one instance of this element
                        reader.Skip();
                        break;
                    case "acquisitionSettingsList": // mzML 1.0.0 compatibility
                        // Schema requirements: zero to one instances of this element
                        reader.Skip();
                        break;
                    case "run":
                        // Schema requirements: one instance of this element
                        // Use reader.ReadSubtree() to provide an XmlReader that is only valid for the element and child nodes
                        ReadRunData(reader.ReadSubtree());
                        if (_randomAccess || _reduceMemoryUsage)
                        {
                            // Kill the read, since we already have a valid index
                            continueReading = false;
                            // don't worry about the skip, since it can take some time.
                            //reader.Skip();
                        }
                        else
                        {
                            // "run" might not have any child nodes
                            // We will either consume the EndElement, or the same element that was passed to ReadRunData (in case of no child nodes)
                            reader.Read();
                        }
                        break;
                    default:
                        // We are not reading anything out of the tag, so bypass it
                        reader.Skip();
                        break;
                }
            }
            _haveMetaData = true;
            if (!_randomAccess && !_reduceMemoryUsage)
            {
                _allRead = true;
            }
            //_numSpectra = _spectrumOffsets.Offsets.Count;
            /* // Now read before any of the metadata.
            if (indexReader != null)
            {
                reader = indexReader;
                //_reader.ReadStartElement("mzML");
                // Read the next node - should be the first child node
                while (reader.ReadState == ReadState.Interactive)
                {
                    // Handle exiting out properly at EndElement tags
                    if (reader.NodeType != XmlNodeType.Element)
                    {
                        reader.Read();
                        continue;
                    }
                    // Handle each 1st level as a chunk
                    switch (reader.Name)
                    {
                        case "indexList":
                            // Schema requirements: one instance of this element
                            ReadIndexList(reader.ReadSubtree());
                            reader.ReadEndElement(); // "fileDescription" must have child nodes
                            break;
                        case "indexListOffset":
                            // Schema requirements: zero to one instances of this element
                            _indexListOffset = Int64.Parse(reader.ReadElementContentAsString());
                            break;
                        case "fileChecksum":
                            // Schema requirements: zero to one instances of this element
                            reader.Skip();
                            break;
                        default:
                            // We are not reading anything out of the tag, so bypass it
                            reader.Skip();
                            break;
                    }
                }
                reader.Dispose();
            } */
            if (!_reduceMemoryUsage)
            {
                // Don't worry about closing the subtree readers, just close the root reader.
                // reader is the root if it is not an indexed mzML file.
                if (indexReader == null)
                {
                    reader.Dispose();
                }
                else
                {
                    indexReader.Dispose();
                }
            }
        }

        #endregion

        #region Metadata tag readers

        /// <summary>
        /// Handle the child nodes of the fileDescription element
        /// Called by ReadMzML (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "fileDescription" element</param>
        private void ReadFileDescription(XmlReader reader)
        {
            reader.MoveToContent();
            reader.ReadStartElement("fileDescription"); // Throws exception if we are not at the "fileDescription" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                switch (reader.Name)
                {
                    case "fileContent":
                        // Schema requirements: one instance of this element
                        reader.Skip();
                        break;
                    case "sourceFileList":
                        // Schema requirements: zero to one instances of this element
                        ReadSourceFileList(reader.ReadSubtree());
                        reader.ReadEndElement(); // "sourceFileList" must have child nodes
                        break;
                    case "contact":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();
        }

        /// <summary>
        /// Handle a single sourceFileList element and child nodes
        /// Called by ReadMzML (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single sourceFileList element</param>
        private void ReadSourceFileList(XmlReader reader)
        {
            reader.MoveToContent();
            //int count = Convert.ToInt32(reader.GetAttribute("count"));
            reader.ReadStartElement("sourceFileList"); // Throws exception if we are not at the "sourceFileList" tag.

            while (reader.ReadState == ReadState.Interactive && (_nativeIdFormat == CV.CVID.CVID_Unknown || _nativeFormat == CV.CVID.CVID_Unknown))
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                if (reader.Name == "sourceFile")
                {
                    // Schema requirements: one to many instances of this element
                    var innerReader = reader.ReadSubtree();
                    innerReader.MoveToContent();
                    innerReader.ReadStartElement("sourceFile"); // Throws exception if we are not at the "sourceFile" tag.
                    while (innerReader.ReadState == ReadState.Interactive)
                    {
                        // Handle exiting out properly at EndElement tags
                        if (innerReader.NodeType != XmlNodeType.Element)
                        {
                            innerReader.Read();
                            continue;
                        }
                        switch (innerReader.Name)
                        {
                            case "referenceableParamGroupRef":
                                // Schema requirements: zero to many instances of this element
                                innerReader.Skip();
                                break;
                            case "cvParam":
                                // Schema requirements: zero to many instances of this element
                                // ReSharper disable CommentTypo
                                /* MUST supply a *child* term of MS:1000767 (native spectrum identifier format) only once
                                 *   e.g.: MS:1000768 (Thermo nativeID format)
                                 *   e.g.: MS:1000769 (Waters nativeID format)
                                 *   e.g.: MS:1000770 (WIFF nativeID format)
                                 *   e.g.: MS:1000771 (Bruker/Agilent YEP nativeID format)
                                 *   e.g.: MS:1000772 (Bruker BAF nativeID format)
                                 *   e.g.: MS:1000773 (Bruker FID nativeID format)
                                 *   e.g.: MS:1000774 (multiple peak list nativeID format)
                                 *   e.g.: MS:1000775 (single peak list nativeID format)
                                 *   e.g.: MS:1000776 (scan number only nativeID format)
                                 *   e.g.: MS:1000777 (spectrum identifier nativeID format)
                                 *   e.g.: MS:1000823 "Bruker U2 nativeID format"
                                 *   e.g.: MS:1000824 "no nativeID format"
                                 *   e.g.: MS:1000929 "Shimadzu Biotech nativeID format"
                                 *   e.g.: MS:1001480 "AB SCIEX TOF/TOF nativeID format"
                                 *   e.g.: MS:1001508 "Agilent MassHunter nativeID format"
                                 *   e.g.: MS:1001526 "spectrum from database integer nativeID format"
                                 *   e.g.: MS:1001528 "Mascot query number"
                                 *   e.g.: MS:1001531 "spectrum from ProteinScape database nativeID format"
                                 *   e.g.: MS:1001532 "spectrum from database string nativeID format"
                                 *   e.g.: MS:1001559 "AB SCIEX TOF/TOF T2D nativeID format"
                                 *   e.g.: MS:1001562 "Scaffold nativeID format"
                                 *   e.g.: MS:1002303 "Bruker Container nativeID format"
                                 *   etc.
                                 * MUST supply a *child* term of MS:1000561 (data file checksum type) one or more times
                                 *   e.g.: MS:1000568 (MD5)
                                 *   e.g.: MS:1000569 (SHA-1)
                                 * MUST supply a *child* term of MS:1000560 (mass spectrometer file format) only once
                                 *   e.g.: MS:1000526 (Waters raw file)             1.0.0: (MassLynx raw format)
                                 *   e.g.: MS:1000562 (ABI WIFF file)               1.0.0: (wiff file)
                                 *   e.g.: MS:1000563 (Thermo RAW file)             1.0.0: (Xcalibur RAW file)
                                 *   e.g.: MS:1000564 (PSI mzData file)             1.0.0: (mzData file)
                                 *   e.g.: MS:1000565 (Micromass PKL file)          1.0.0: (pkl file)
                                 *   e.g.: MS:1000566 (ISB mzXML file)              1.0.0: (mzXML file)
                                 *   e.g.: MS:1000567 (Bruker/Agilent YEP file)     1.0.0: (yep file)
                                 *   e.g.: MS:1000584 (mzML file)
                                 *   e.g.: MS:1000613 (DTA file)                    1.0.0: (dta file)
                                 *   e.g.: MS:1000614 (ProteinLynx Global Server mass spectrum XML file)
                                 *   e.g.: MS:1000740 "parameter file"
                                 *   e.g.: MS:1000742 "Bioworks SRF format"
                                 *   e.g.: MS:1000815 "Bruker BAF format"
                                 *   e.g.: MS:1000816 "Bruker U2 format"
                                 *   e.g.: MS:1000825 "Bruker FID format"
                                 *   e.g.: MS:1000930 "Shimadzu Biotech database entity"
                                 *   e.g.: MS:1001062 "Mascot MGF format"
                                 *   e.g.: MS:1001245 "PerSeptive PKS format"
                                 *   e.g.: MS:1001246 "Sciex API III format"
                                 *   e.g.: MS:1001247 "Bruker XML format"
                                 *   e.g.: MS:1001369 "text format"
                                 *   e.g.: MS:1001463 "Phenyx XML format"
                                 *   e.g.: MS:1001481 "AB SCIEX TOF/TOF database"
                                 *   e.g.: MS:1001509 "Agilent MassHunter format"
                                 *   e.g.: MS:1001527 "Proteinscape spectra"
                                 *   e.g.: MS:1001560 "AB SCIEX TOF/TOF T2D format"
                                 *   e.g.: MS:1001881 "mz5 format"
                                 *   e.g.: MS:1002302 "Bruker Container format"
                                 *   e.g.: MS:1002385 "SCiLS Lab format"
                                 *   e.g.: MS:1002441 "Andi-MS format"
                                 *   etc.
                                 */
                                // ReSharper restore CommentTypo
                                var cv = innerReader.GetAttribute("cvRef");
                                var accession = innerReader.GetAttribute("accession");
                                if (string.IsNullOrWhiteSpace(cv))
                                {
                                    cv = "MS";
                                }
                                if (cv.IndexOf("PSI", StringComparison.OrdinalIgnoreCase) >= 0)
                                {
                                    cv = "MS";
                                }
                                if (!string.IsNullOrWhiteSpace(accession) && CV.TermAccessionLookup[cv].ContainsKey(accession))
                                {
                                    var cvID = CV.TermAccessionLookup[cv][accession];
                                    if (CV.CvidIsA(cvID, CV.CVID.MS_native_spectrum_identifier_format))
                                    {
                                        _nativeIdFormat = cvID;
                                    }
                                    else if (CV.CvidIsA(cvID, CV.CVID.MS_mass_spectrometer_file_format))
                                    {
                                        _nativeFormat = cvID;
                                    }
                                }
                                /*switch (innerReader.GetAttribute("accession"))
                                {
                                    case "MS:1000768":
                                        // name="Thermo nativeID format"
                                        _instrument = Instrument.Thermo_RAW_format;
                                        break;
                                    case "MS:1000769":
                                        // name="Waters nativeID format"
                                        _instrument = Instrument.Waters_raw_format;
                                        break;
                                    case "MS:1000770":
                                        // name="WIFF nativeID format"
                                        _instrument = Instrument.ABI_WIFF_format;
                                        break;
                                    case "MS:1000771":
                                        // name="Bruker/Agilent YEP nativeID format"
                                        break;
                                    case "MS:1000772":
                                        // name="Bruker BAF nativeID format"
                                        break;
                                    case "MS:1000773":
                                        // name="Bruker FID nativeID format"
                                        break;
                                    case "MS:1000774":
                                        // name="multiple peak list nativeID format"
                                        break;
                                    case "MS:1000775":
                                        // name="single peak list nativeID format"
                                        break;
                                    case "MS:1000776":
                                        // name="scan number only nativeID format"
                                        break;
                                    case "MS:1000777":
                                        // name="spectrum identifier nativeID format"
                                        break;
                                    case "MS:1000823":
                                        // name="Bruker U2 nativeID format"
                                        break;
                                    case "MS:1000824":
                                        // name="no nativeID format"
                                        break;
                                    case "MS:1000929":
                                        // name="Shimadzu Biotech nativeID format"
                                        break;
                                    case "MS:1001480":
                                        // name="AB SCIEX TOF/TOF nativeID format"
                                        break;
                                    case "MS:1001508":
                                        // name="Agilent MassHunter nativeID format"
                                        break;
                                    case "MS:1001526":
                                        // name="spectrum from database integer nativeID format"
                                        break;
                                    case "MS:1001528":
                                        // name="Mascot query number"
                                        break;
                                    case "MS:1001531":
                                        // name="spectrum from ProteinScape database nativeID format"
                                        break;
                                    case "MS:1001532":
                                        // name="spectrum from database string nativeID format"
                                        break;
                                    case "MS:1001559":
                                        // name="AB SCIEX TOF/TOF T2D nativeID format"
                                        break;
                                    case "MS:1001562":
                                        // name="Scaffold nativeID format"
                                        break;
                                    case "MS:1002303":
                                        // name="Bruker Container nativeID format"
                                        break;
                                }*/
                                innerReader.Read(); // Consume the cvParam element (no child nodes)
                                break;
                            case "userParam":
                                // Schema requirements: zero to many instances of this element
                                innerReader.Skip();
                                break;
                            default:
                                innerReader.Skip();
                                break;
                        }
                    }
                    innerReader.Dispose();
                    reader.Read();
                }
                else
                {
                    reader.Read();
                }
            }
            reader.Dispose();
        }

        /// <summary>
        /// Handle the child nodes of the referenceableParamGroupList element
        /// Called by ReadMzML (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "referenceableParamGroupList" element</param>
        private void ReadReferenceableParamGroupList(XmlReader reader)
        {
            _referenceableParamGroups.Clear(); // In case of second read of file, clear out existing.
            reader.MoveToContent();
            // var count = Convert.ToInt32(reader.GetAttribute("count"));
            reader.ReadStartElement("referenceableParamGroupList"); // Throws exception if we are not at the "referenceableParamGroupList" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                if (reader.Name == "referenceableParamGroup")
                {
                    // Schema requirements: one to many instances of this element
                    var id = reader.GetAttribute("id");
                    var paramList = new List<Param>();
                    var innerReader = reader.ReadSubtree();
                    innerReader.MoveToContent();
                    innerReader.ReadStartElement("referenceableParamGroup"); // Throws exception if we are not at the "sourceFile" tag.
                    while (innerReader.ReadState == ReadState.Interactive)
                    {
                        // Handle exiting out properly at EndElement tags
                        if (innerReader.NodeType != XmlNodeType.Element)
                        {
                            innerReader.Read();
                            continue;
                        }
                        switch (innerReader.Name)
                        {
                            case "cvParam":
                                // Schema requirements: zero to many instances of this element
                                paramList.Add(ReadCvParam(innerReader.ReadSubtree()));
                                innerReader.Read(); // Consume the cvParam element (no child nodes)
                                break;
                            case "userParam":
                                // Schema requirements: zero to many instances of this element
                                paramList.Add(ReadUserParam(innerReader.ReadSubtree()));
                                innerReader.Read(); // Consume the userParam element (no child nodes)
                                break;
                            default:
                                innerReader.Skip();
                                break;
                        }
                    }
                    innerReader.Dispose();
                    reader.Read();
                    if (string.IsNullOrWhiteSpace(id))
                    {
                        ConsoleMsgUtils.ShowWarning("Encountered referenceable param group with null or empty id");
                    }
                    else
                    {
                        _referenceableParamGroups.Add(id, paramList);
                    }
                }
                else
                {
                    reader.Read();
                }
            }
            reader.Dispose();
        }

        /// <summary>
        /// Handle the cvParam element
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "cvParam" element</param>
        private CVParam ReadCvParam(XmlReader reader)
        {
            reader.MoveToContent();
            var cvParam = new CVParam
            {
                Accession = reader.GetAttribute("accession"),
                CVRef = reader.GetAttribute("cvRef"),
                Name = reader.GetAttribute("name"),
                Value = reader.GetAttribute("value"),
                UnitAccession = reader.GetAttribute("unitAccession"),
                UnitCVRef = reader.GetAttribute("unitCVRef"),
                UnitName = reader.GetAttribute("unitName")
            };

            reader.Dispose();
            return cvParam;
        }

        /// <summary>
        /// Handle the userParam element
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "userParam" element</param>
        private UserParam ReadUserParam(XmlReader reader)
        {
            reader.MoveToContent();
            var userParam = new UserParam
            {
                Name = reader.GetAttribute("name"),
                Type = reader.GetAttribute("type"),
                Value = reader.GetAttribute("value"),
                UnitAccession = reader.GetAttribute("unitAccession"),
                UnitCVRef = reader.GetAttribute("unitCVRef"),
                UnitName = reader.GetAttribute("unitName")
            };

            reader.Dispose();
            return userParam;
        }

        #endregion

        #region Run and SpectrumList Tags

        /// <summary>
        /// Handle the child nodes of the run element
        /// Called by ReadMzML (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "run" element</param>
        private void ReadRunData(XmlReader reader)
        {
            reader.MoveToContent();
            reader.ReadStartElement("run"); // Throws exception if we are not at the "run" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "sourceFileRefList": // mzML_1.0.0 compatibility
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "spectrumList":
                        // Schema requirements: zero to one instances of this element
                        // Use reader.ReadSubtree() to provide an XmlReader that is only valid for the element and child nodes
                        ReadSpectrumList(reader.ReadSubtree());
                        if (_randomAccess || _reduceMemoryUsage)
                        {
                            // Don't worry about reading anything more, and closing the XmlReader will take more time than it is worth.
                            return;
                        }
                        // "spectrumList" might not have any child nodes
                        // We will either consume the EndElement, or the same element that was passed to ReadSpectrumList (in case of no child nodes)
                        reader.Read();
                        break;
                    case "chromatogramList":
                        // Schema requirements: zero to one instances of this element
                        reader.Skip();
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();
        }

        /// <summary>
        /// Handle the child nodes of a spectrumList element
        /// Called by ReadRunData (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single spectrumList element</param>
        private void ReadSpectrumList(XmlReader reader)
        {
            reader.MoveToContent();
            _numSpectra = Convert.ToInt64(reader.GetAttribute("count"));
            if (_randomAccess)
            {
                // randomAccess: We only read to this point for the count of spectra.
                // We only want to read for offsets if we weren't able to get valid offsets from an index
                //reader.Dispose(); // Closing can be slow for a subtree...
                if (!_haveIndex)
                {
                    ReadRunForOffsets();
                }
                return;
            }
            reader.ReadStartElement("spectrumList"); // Throws exception if we are not at the "SpectrumIdentificationList" tag.
            if (_reduceMemoryUsage)
            {
                // Kill the read, we are at the first spectrum
                _xmlReaderForYield = reader;
                // If in the "ReadAllSpectra" call stack, we don't want the reader closed - we still need the subtree
                return;
            }
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                if (reader.Name == "spectrum")
                {
                    // Schema requirements: zero to many instances of this element
                    // Use reader.ReadSubtree() to provide an XmlReader that is only valid for the element and child nodes
                    _spectra.Add(ReadSpectrum(reader.ReadSubtree()));
                    // "spectrum" might not have any child nodes
                    // We will either consume the EndElement, or the same element that was passed to ReadSpectrum (in case of no child nodes)
                    reader.Read();
                }
                else
                {
                    reader.Skip();
                }
            }
            reader.Dispose();
        }

        #endregion

        #region Spectrum Tag

        /// <summary>
        /// Handle a single spectrum element and child nodes
        /// Called by ReadSpectrumList (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single spectrum element</param>
        /// <param name="includePeaks">Whether to read binary data arrays</param>
        /// <remarks>If the data is not centroided, this method uses InformedProteomics.Backend.Utils.Centroider to find peaks</remarks>
        private Spectrum ReadSpectrum(XmlReader reader, bool includePeaks = true)
        {
            reader.MoveToContent();
            var index = reader.GetAttribute("index");
            //Console.WriteLine("Reading spectrum indexed by " + index);
            // This is correct for Thermo files converted by msConvert, but need to implement for others as well
            var spectrumId = reader.GetAttribute("id"); // Native ID in mzML_1.1.0; unique identifier in mzML_1.0.0, often same as nativeID
            var nativeId = spectrumId;
            if (_version == MzML_Version.mzML1_0_0)
            {
                nativeId = reader.GetAttribute("nativeID"); // Native ID in mzML_1.0.0
            }

            if (string.IsNullOrWhiteSpace(nativeId))
            {
                throw new Exception("nativeID is empty for spectrum index " + index);
            }

            // Use the index (0-based) as a "scan number" (the first spectrum in the file typically has scan number = 1)
            // It's the only reliable way to have the scan number be contiguous and have no duplicates
            var scanNum = Convert.ToInt32(index) + 1;

            // Converting the NativeId to a scan number is nice and all, but also problematic:
            // * Waters .raw files have multiple functions (and spectra) per scan number, unless a filter is used to only select one function
            // * Agilent datasets often have non-contiguous scan numbers (in NativeID)
            // * UIMF datasets have multiple frames, and each frame has its own list of scans
            // * Other examples undoubtedly exist
            // Better to just rely on the "index + 1" as a 1-based index
            // If a random access reader, there is already a scan number stored, based on the order of the index. Use it instead.
            //if (_randomAccess)
            //{
            //    scanNum = (int)(_spectrumOffsets.NativeToIdMap[nativeId]);
            //}
            //else
            //{
            //    //scanNum = (int)(_artificialScanNum++);
            //    scanNum = Convert.ToInt32(index) + 1;
            //    // Interpret the NativeID (if the format has an interpreter) and use it instead of the artificial number.
            //    // TODO: Better handling than the artificial ID for other nativeIDs (ones currently not supported)
            //    if (NativeIdConversion.TryGetScanNumberInt(nativeId, out var value))
            //    {
            //        scanNum = value;
            //    }
            //}

            var defaultArraySize = Convert.ToInt32(reader.GetAttribute("defaultArrayLength"));
            reader.ReadStartElement("spectrum"); // Throws exception if we are not at the "spectrum" tag.
            var is_ms_ms = false;
            var msLevel = 0;
            var centroided = false;
            double tic = 0;
            var precursors = new List<Precursor>();
            var scans = new List<ScanData>();
            var binaryData = new List<BinaryDataArray>();
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                //////////////////////////////////////////////////////////////////////////////////////
                //
                // MS1 Spectra: only need Spectrum data: scanNum, MSLevel, ElutionTime, mzArray, IntensityArray
                //
                // MS2 Spectra: use ProductSpectrum; adds ActivationMethod and IsolationWindow
                //
                //////////////////////////////////////////////////////////////////////////////////////
                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        /* MAY supply a *child* term of MS:1000465 (scan polarity) only once
                         *   e.g.: MS:1000129 (negative scan)
                         *   e.g.: MS:1000130 (positive scan)
                         * MUST supply a *child* term of MS:1000559 (spectrum type) only once
                         *   e.g.: MS:1000322 (charge inversion mass spectrum)
                         *   e.g.: MS:1000325 (constant neutral gain spectrum)
                         *   e.g.: MS:1000326 (constant neutral loss spectrum)
                         *   e.g.: MS:1000328 (e/2 mass spectrum)
                         *   e.g.: MS:1000341 (precursor ion spectrum)
                         *   e.g.: MS:1000579 (MS1 spectrum)
                         *   e.g.: MS:1000580 (MSn spectrum)
                         *   e.g.: MS:1000581 (CRM spectrum)
                         *   e.g.: MS:1000582 (SIM spectrum)
                         *   e.g.: MS:1000583 (SRM spectrum)
                         *   e.g.: MS:1000620 (PDA spectrum)
                         *   e.g.: MS:1000627 (selected ion current chromatogram)
                         *   e.g.: MS:1000789 (enhanced multiply charged spectrum)
                         *   e.g.: MS:1000790 (time-delayed fragmentation spectrum)
                         *   etc.
                         * MUST supply term MS:1000525 (spectrum representation) or any of its children only once
                         *   e.g.: MS:1000127 (centroid spectrum)
                         *   e.g.: MS:1000128 (profile spectrum)
                         * MAY supply a *child* term of MS:1000499 (spectrum attribute) one or more times
                         *   e.g.: MS:1000285 (total ion current)
                         *   e.g.: MS:1000497 (zoom scan)
                         *   e.g.: MS:1000504 (base peak m/z)
                         *   e.g.: MS:1000505 (base peak intensity)
                         *   e.g.: MS:1000511 (ms level)
                         *   e.g.: MS:1000527 (highest observed m/z)
                         *   e.g.: MS:1000528 (lowest observed m/z)
                         *   e.g.: MS:1000618 (highest observed wavelength)
                         *   e.g.: MS:1000619 (lowest observed wavelength)
                         *   e.g.: MS:1000796 (spectrum title)
                         *   etc.
                         */
                        switch (reader.GetAttribute("accession"))
                        {
                            case "MS:1000127":
                                // name="centroid spectrum"
                                centroided = true;
                                break;
                            case "MS:1000128":
                                // name="profile spectrum"
                                centroided = false;
                                break;
                            case "MS:1000511":
                                // name="ms level"
                                msLevel = Convert.ToInt32(reader.GetAttribute("value"));
                                break;
                            case "MS:1000579":
                                // name="MS1 spectrum"
                                is_ms_ms = false;
                                break;
                            case "MS:1000580":
                                // name="MSn spectrum"
                                is_ms_ms = true;
                                break;
                            case "MS:1000285":
                                // name="total ion current"
                                tic = Convert.ToDouble(reader.GetAttribute("value"));
                                break;
                        }
                        reader.Read(); // Consume the cvParam element (no child nodes)
                        break;
                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "spectrumDescription": // mzML_1.0.0 compatibility
                        // Schema requirements: one instance of this element
                        ReadSpectrumDescription(reader.ReadSubtree(), ref scans, ref precursors, out centroided);
                        reader.ReadEndElement(); // "spectrumDescription" must have child nodes
                        break;
                    case "scanList":
                        // Schema requirements: zero to one instances of this element
                        scans.AddRange(ReadScanList(reader.ReadSubtree()));
                        reader.ReadEndElement(); // "scanList" must have child nodes
                        break;
                    case "precursorList":
                        // Schema requirements: zero to one instances of this element
                        precursors.AddRange(ReadPrecursorList(reader.ReadSubtree()));
                        reader.ReadEndElement(); // "precursorList" must have child nodes
                        break;
                    case "productList":
                        // Schema requirements: zero to one instances of this element
                        reader.Skip();
                        break;
                    case "binaryDataArrayList":
                        // Schema requirements: zero to one instances of this element
                        if (includePeaks)
                        {
                            binaryData.AddRange(ReadBinaryDataArrayList(reader.ReadSubtree(), defaultArraySize));
                            reader.ReadEndElement(); // "binaryDataArrayList" must have child nodes
                        }
                        else
                        {
                            reader.Skip();
                        }
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();

            // Process the spectrum data
            Spectrum spectrum;
            var mzs = new BinaryDataArray();
            var intensities = new BinaryDataArray();
            foreach (var dataArray in binaryData)
            {
                if (dataArray.ArrayType == ArrayType.m_z_array)
                {
                    mzs = dataArray;
                }
                else if (dataArray.ArrayType == ArrayType.intensity_array)
                {
                    intensities = dataArray;
                }
            }

            if (!centroided && includePeaks)
            {
                // Centroid spectrum
                // ProteoWizard
                var centroider = new Centroider(mzs.Data, intensities.Data);
                centroider.GetCentroidedData(out var centroidedMzs, out var centroidedIntensities);
                mzs.Data = centroidedMzs;
                intensities.Data = centroidedIntensities;
            }

            ScanData scan;
            if (scans.Count == 1)
            {
                scan = scans[0];
            }
            else if (scans.Count > 1)
            {
                // TODO: Should do something else to appropriately handle combinations...
                scan = scans[0];
            }
            else
            {
                scan = new ScanData();
            }

            if (is_ms_ms)
            {
                var precursor = new Precursor();
                if (precursors.Count == 1)
                {
                    precursor = precursors[0];
                }
                else if (precursors.Count > 1)
                {
                    // TODO: Should do something else to appropriately handle multiple precursors...
                    precursor = precursors[0];
                }
                var ion = new SelectedIon();
                if (precursor.Ions.Count == 1)
                {
                    ion = precursor.Ions[0];
                }
                else if (precursor.Ions.Count > 1)
                {
                    // TODO: Should do something else to appropriately handle multiple selected ions...
                    ion = precursor.Ions[0];
                }

                var productSpectrum = new ProductSpectrum(mzs.Data, intensities.Data, scanNum)
                {
                    ActivationMethod = precursor.Activation
                };

                // Select mz value to use based on presence of a Thermo-specific user param.
                // The user param has a slightly higher precision, if that matters.
                var mz = Math.Abs(scan.MonoisotopicMz) < float.Epsilon ? ion.SelectedIonMz : scan.MonoisotopicMz;
                var isoMz = precursor.IsolationWindowTargetMz > float.Epsilon ? precursor.IsolationWindowTargetMz : mz;

                productSpectrum.IsolationWindow = new IsolationWindow(isoMz, precursor.IsolationWindowLowerOffset, precursor.IsolationWindowUpperOffset, mz, ion.Charge);
                // productSpectrum.IsolationWindow.OldCharge = ion.OldCharge;
                // productSpectrum.IsolationWindow.SelectedIonMz = ion.SelectedIonMz;
                spectrum = productSpectrum;
            }
            else
            {
                spectrum = new Spectrum(mzs.Data, intensities.Data, scanNum);
            }
            spectrum.MsLevel = msLevel;
            spectrum.ElutionTime = scan.StartTime;
            spectrum.NativeId = nativeId;
            spectrum.TotalIonCurrent = tic;

            return spectrum;
        }

        #endregion

        #region Spectrum internal Tags

        /// <summary>
        /// mzML_1.0.0 compatibility
        /// Handle a single spectrumDescription element and child nodes
        /// Called by ReadSpectrumList (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single spectrum element</param>
        /// <param name="scans"></param>
        /// <param name="precursors"></param>
        /// <param name="centroided"></param>
        private void ReadSpectrumDescription(XmlReader reader, ref List<ScanData> scans, ref List<Precursor> precursors, out bool centroided)
        {
            reader.MoveToContent();
            // This is correct for Thermo files converted by msConvert, but need to implement for others as well
            reader.ReadStartElement("spectrumDescription"); // Throws exception if we are not at the "spectrumDescription" tag.
            centroided = false;
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }
                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        /* MUST supply a *child* term of MS:1000525 (spectrum representation) only once
                         *   e.g.: MS:1000127 (centroid mass spectrum)
                         *   e.g.: MS:1000128 (profile mass spectrum)
                         * MUST supply a *child* term of MS:1000499 (spectrum attribute) one or more times
                         *   e.g.: MS:1000285 (total ion current)
                         *   e.g.: MS:1000504 (base peak m/z)
                         *   e.g.: MS:1000505 (base peak intensity)
                         *   e.g.: MS:1000527 (highest m/z value)
                         *   e.g.: MS:1000528 (lowest m/z value)
                         *   e.g.: MS:1000618 (highest wavelength value)
                         *   e.g.: MS:1000619 (lowest wavelength value)
                         */
                        switch (reader.GetAttribute("accession"))
                        {
                            case "MS:1000127":
                                // name="centroid spectrum"
                                centroided = true;
                                break;
                        }
                        reader.Read(); // Consume the cvParam element (no child nodes)
                        break;
                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "acquisitionList":
                        // Schema requirements: zero to one instances of this element
                        // Very comparable to mzML_1.1.0's scanList. Use it.
                        scans.AddRange(ReadScanList(reader.ReadSubtree()));
                        reader.ReadEndElement(); // "acquisitionList" must have child nodes
                        break;
                    case "precursorList":
                        // Schema requirements: zero to one instances of this element
                        precursors.AddRange(ReadPrecursorList(reader.ReadSubtree()));
                        reader.ReadEndElement(); // "precursorList" must have child nodes
                        break;
                    case "scan":
                        // Schema requirements: zero to one instances of this element
                        scans.Add(ReadScan(reader.ReadSubtree()));
                        reader.Read(); // "scan" might not have child nodes
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();
        }

        /// <summary>
        /// Handle a single scanList element and child nodes
        /// Called by ReadSpectrum (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single scanList element</param>
        /// <returns>Enumerable list of scan data objects</returns>
        private IEnumerable<ScanData> ReadScanList(XmlReader reader)
        {
            reader.MoveToContent();
            // var count = Convert.ToInt32(reader.GetAttribute("count"));
            var scans = new List<ScanData>();
            if (_version == MzML_Version.mzML1_0_0)
            {
                reader.ReadStartElement("acquisitionList"); // Throws exception if we are not at the "scanList" tag.
            }
            else
            {
                reader.ReadStartElement("scanList"); // Throws exception if we are not at the "scanList" tag.
            }
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        /* MUST supply a *child* term of MS:1000570 (spectra combination) only once
                         *   e.g.: MS:1000571 (sum of spectra)
                         *   e.g.: MS:1000573 (median of spectra)
                         *   e.g.: MS:1000575 (mean of spectra)
                         *   e.g.: MS:1000795 (no combination)
                         */
                        switch (reader.GetAttribute("accession"))
                        {
                            case "MS:1000571":
                                // name="sum of spectra"
                                break;
                            case "MS:1000573":
                                // name="median of spectra"
                                break;
                            case "MS:1000575":
                                // name="mean of spectra"
                                break;
                            case "MS:1000795":
                                // name="no combination"
                                break;
                        }
                        reader.Read(); // Consume the cvParam element (no child nodes)
                        break;
                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "scan":
                    case "acquisition": // mzML_1.0.0 compatibility
                        // Schema requirements: one to many instances of this element
                        scans.Add(ReadScan(reader.ReadSubtree()));
                        reader.Read(); // "scan" might not have child nodes
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();
            return scans;
        }

        /// <summary>
        /// Handle a single scan element and child nodes
        /// Called by ReadSpectrum (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single scan element</param>
        /// <returns>Scan data object</returns>
        private ScanData ReadScan(XmlReader reader)
        {
            reader.MoveToContent();
            if (_version == MzML_Version.mzML1_0_0)
            {
                var name = reader.Name;
                if (!name.Equals("scan") && !name.Equals("acquisition"))
                {
                    throw new XmlException("Invalid schema");
                }
                reader.ReadStartElement(name);
            }
            else
            {
                reader.ReadStartElement("scan"); // Throws exception if we are not at the "scan" tag.
            }
            var scan = new ScanData();
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        reader.Skip();
                        break;
                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        /* MAY supply a *child* term of MS:1000503 (scan attribute) one or more times
                         *   e.g.: MS:1000011 (mass resolution)
                         *   e.g.: MS:1000015 (scan rate)
                         *   e.g.: MS:1000016 (scan start time)
                         *   e.g.: MS:1000502 (dwell time)
                         *   e.g.: MS:1000512 (filter string)
                         *   e.g.: MS:1000616 (preset scan configuration)
                         *   e.g.: MS:1000800 (mass resolving power)
                         *   e.g.: MS:1000803 (analyzer scan offset)
                         *   e.g.: MS:1000826 (elution time)
                         *   e.g.: MS:1000880 (interchannel delay)
                         * MAY supply a *child* term of MS:1000018 (scan direction) only once
                         *   e.g.: MS:1000092 (decreasing m/z scan)
                         *   e.g.: MS:1000093 (increasing m/z scan)
                         * MAY supply a *child* term of MS:1000019 (scan law) only once
                         *   e.g.: MS:1000094 (exponential)
                         *   e.g.: MS:1000095 (linear)
                         *   e.g.: MS:1000096 (quadratic)
                         */
                        switch (reader.GetAttribute("accession"))
                        {
                            case "MS:1000016":
                                // name="scan start time"
                                var time = Convert.ToDouble(reader.GetAttribute("value"));
                                var isSeconds = reader.GetAttribute("unitName") == "second";
                                // Should only see "second" and "minute"
                                scan.StartTime = isSeconds ? time / 60.0 : time;
                                //scan.StartTime = Convert.ToDouble(reader.GetAttribute("value"));
                                break;
                            case "MS:1000512":
                                // name="filter string"
                                break;
                            case "MS:1000616":
                                // name="preset scan configuration"
                                break;
                            case "MS:1000927":
                                // name="ion injection time"
                                break;
                            case "MS:1000826":
                                // name="elution time"
                                //startTime = Convert.ToDouble(reader.GetAttribute("value"));
                                break;
                        }
                        reader.Read(); // Consume the cvParam element (no child nodes)
                        break;
                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        if (reader.GetAttribute("name") == "[Thermo Trailer Extra]Monoisotopic M/Z:")
                        {
                            scan.MonoisotopicMz = Convert.ToDouble(reader.GetAttribute("value"));
                        }
                        reader.Read(); // Consume the userParam element (no child nodes)
                        break;
                    case "scanWindowList":
                        // Schema requirements: zero to one instances of this element
                        //ReadScanList(reader.ReadSubtree());
                        //reader.ReadEndElement(); // "scanWindowList" must have child nodes
                        reader.Skip();
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();
            return scan;
        }

        /// <summary>
        /// Handle a single precursorList element and child nodes
        /// Called by ReadSpectrum (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single precursorList element</param>
        /// <returns>Enumerable list of precursor objects</returns>
        private IEnumerable<Precursor> ReadPrecursorList(XmlReader reader)
        {
            reader.MoveToContent();
            // var count = Convert.ToInt32(reader.GetAttribute("count"));
            var precursors = new List<Precursor>();
            reader.ReadStartElement("precursorList"); // Throws exception if we are not at the "precursorList" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "precursor":
                        // Schema requirements: one to many instances of this element
                        precursors.Add(ReadPrecursor(reader.ReadSubtree()));
                        reader.ReadEndElement(); // "SpectrumIdentificationItem" must have child nodes
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();
            return precursors;
        }

        /// <summary>
        /// Handle a single precursor element and child nodes
        /// Called by ReadPrecursorList (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single precursor element</param>
        /// <returns>Precursor object</returns>
        private Precursor ReadPrecursor(XmlReader reader)
        {
            reader.MoveToContent();
            reader.ReadStartElement("precursor"); // Throws exception if we are not at the "precursor" tag.
            var precursor = new Precursor();
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                XmlReader innerReader;
                switch (reader.Name)
                {
                    case "isolationWindow":
                        // Schema requirements: zero to one instances of this element
                        innerReader = reader.ReadSubtree();
                        innerReader.MoveToContent();
                        innerReader.ReadStartElement("isolationWindow"); // Throws exception if we are not at the "selectedIon" tag.
                        while (innerReader.ReadState == ReadState.Interactive)
                        {
                            // Handle exiting out properly at EndElement tags
                            if (innerReader.NodeType != XmlNodeType.Element)
                            {
                                innerReader.Read();
                                continue;
                            }
                            switch (innerReader.Name)
                            {
                                case "referenceableParamGroupRef":
                                    // Schema requirements: zero to many instances of this element
                                    innerReader.Skip();
                                    break;
                                case "cvParam":
                                    // Schema requirements: zero to many instances of this element
                                    /* MUST supply a *child* term of MS:1000792 (isolation window attribute) one or more times
                                     *   e.g.: MS:1000827 (isolation window target m/z)
                                     *   e.g.: MS:1000828 (isolation window lower offset)
                                     *   e.g.: MS:1000829 (isolation window upper offset)
                                     */
                                    switch (innerReader.GetAttribute("accession"))
                                    {
                                        case "MS:1000827":
                                            // name="isolation window target m/z"
                                            precursor.IsolationWindowTargetMz = Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                        case "MS:1000828":
                                            // name="isolation window lower offset"
                                            precursor.IsolationWindowLowerOffset = Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                        case "MS:1000829":
                                            // name="isolation window upper offset"
                                            precursor.IsolationWindowUpperOffset = Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                    }
                                    innerReader.Read(); // Consume the cvParam element (no child nodes)
                                    break;
                                case "userParam":
                                    // Schema requirements: zero to many instances of this element
                                    innerReader.Skip();
                                    break;
                                default:
                                    innerReader.Skip();
                                    break;
                            }
                        }
                        innerReader.Dispose();
                        reader.Read(); // "selectedIon" might not have child nodes
                        // We will either consume the EndElement, or the same element that was passed to ReadSpectrum (in case of no child nodes)
                        break;
                    case "selectedIonList":
                        // Schema requirements: zero to one instances of this element
                        // mzML_1.0.0: one instance of this element
                        precursor.Ions = ReadSelectedIonList(reader.ReadSubtree());
                        reader.ReadEndElement(); // "selectedIonList" must have child nodes
                        break;
                    case "activation":
                        // Schema requirements: one instance of this element
                        var activationMethods = new List<ActivationMethod>();
                        innerReader = reader.ReadSubtree();
                        innerReader.MoveToContent();
                        innerReader.ReadStartElement("activation"); // Throws exception if we are not at the "activation" tag.
                        while (innerReader.ReadState == ReadState.Interactive)
                        {
                            // Handle exiting out properly at EndElement tags
                            if (innerReader.NodeType != XmlNodeType.Element)
                            {
                                innerReader.Read();
                                continue;
                            }
                            switch (innerReader.Name)
                            {
                                case "referenceableParamGroupRef":
                                    // Schema requirements: zero to many instances of this element
                                    innerReader.Skip();
                                    break;
                                case "cvParam":
                                    // Schema requirements: zero to many instances of this element
                                    /* MAY supply a *child* term of MS:1000510 (precursor activation attribute) one or more times
                                     *   e.g.: MS:1000045 (collision energy)
                                     *   e.g.: MS:1000138 (percent collision energy)
                                     *   e.g.: MS:1000245 (charge stripping)
                                     *   e.g.: MS:1000412 (buffer gas)
                                     *   e.g.: MS:1000419 (collision gas)
                                     *   e.g.: MS:1000509 (activation energy)
                                     *   e.g.: MS:1000869 (collision gas pressure)
                                     * MUST supply term MS:1000044 (dissociation method) or any of its children one or more times
                                     *   e.g.: MS:1000133 (collision-induced dissociation)
                                     *   e.g.: MS:1000134 (plasma desorption)
                                     *   e.g.: MS:1000135 (post-source decay)
                                     *   e.g.: MS:1000136 (surface-induced dissociation)
                                     *   e.g.: MS:1000242 (blackbody infrared radiative dissociation)
                                     *   e.g.: MS:1000250 (electron capture dissociation)
                                     *   e.g.: MS:1000262 (infrared multiphoton dissociation)
                                     *   e.g.: MS:1000282 (sustained off-resonance irradiation)
                                     *   e.g.: MS:1000422 (high-energy collision-induced dissociation)
                                     *   e.g.: MS:1000433 (low-energy collision-induced dissociation)
                                     *   etc.
                                     *
                                     *   e.g.: MS:1000133 "collision-induced dissociation"
                                     *   e.g.: MS:1000134 "plasma desorption"
                                     *   e.g.: MS:1000135 "post-source decay"
                                     *   e.g.: MS:1000136 "surface-induced dissociation"
                                     *   e.g.: MS:1000242 "blackbody infrared radiative dissociation"
                                     *   e.g.: MS:1000250 "electron capture dissociation"
                                     *   e.g.: MS:1000262 "infrared multiphoton dissociation"
                                     *   e.g.: MS:1000282 "sustained off-resonance irradiation"
                                     *   e.g.: MS:1000422 "beam-type collision-induced dissociation"
                                     *   e.g.: MS:1000433 "low-energy collision-induced dissociation"
                                     *   e.g.: MS:1000435 "photodissociation"
                                     *   e.g.: MS:1000598 "electron transfer dissociation"
                                     *   e.g.: MS:1000599 "pulsed q dissociation"
                                     *   e.g.: MS:1001880 "in-source collision-induced dissociation"
                                     *   e.g.: MS:1002000 "LIFT"
                                     *   e.g.: MS:1002472 "trap-type collision-induced dissociation"
                                     */
                                    switch (innerReader.GetAttribute("accession"))
                                    {
                                        case "MS:1000133":
                                            // name="collision-induced dissociation"
                                            activationMethods.Add(ActivationMethod.CID);
                                            break;
                                        case "MS:1000598":
                                            // name="electron transfer dissociation"
                                            activationMethods.Add(ActivationMethod.ETD);
                                            break;
                                        case "MS:1000422":
                                            // name="beam-type collision-induced dissociation", "high-energy collision-induced dissociation"
                                            activationMethods.Add(ActivationMethod.HCD);
                                            break;
                                        case "MS:1000250":
                                            // name="electron capture dissociation"
                                            activationMethods.Add(ActivationMethod.ECD);
                                            break;
                                        case "MS:1000599":
                                            // name="pulsed q dissociation"
                                            activationMethods.Add(ActivationMethod.PQD);
                                            break;
                                        case "MS:1000045":
                                            // name="collision energy"
                                            //energy = Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                    }
                                    innerReader.Read(); // Consume the cvParam element (no child nodes)
                                    break;
                                case "userParam":
                                    // Schema requirements: zero to many instances of this element
                                    innerReader.Skip();
                                    break;
                                default:
                                    innerReader.Skip();
                                    break;
                            }
                        }
                        innerReader.Dispose();
                        reader.Read(); // "selectedIon" might not have child nodes

                        if (activationMethods.Count > 1 && activationMethods.Contains(ActivationMethod.ETD))
                        {
                            precursor.Activation = ActivationMethod.ETD;
                        }
                        else if (activationMethods.Count > 0)
                        {
                            precursor.Activation = activationMethods[0];
                        }
                        // We will either consume the EndElement, or the same element that was passed to ReadSpectrum (in case of no child nodes)
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();
            return precursor;
        }

        /// <summary>
        /// Handle a single selectedIonList element and child nodes
        /// Called by ReadPrecursor (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single selectedIonList element</param>
        /// <returns>List of selected ions</returns>
        private List<SelectedIon> ReadSelectedIonList(XmlReader reader)
        {
            reader.MoveToContent();
            //int count = Convert.ToInt32(reader.GetAttribute("count"));
            var ions = new List<SelectedIon>();
            reader.ReadStartElement("selectedIonList"); // Throws exception if we are not at the "selectedIonList" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "selectedIon":
                        // Schema requirements: one to many instances of this element
                        var ion = new SelectedIon();
                        var innerReader = reader.ReadSubtree();
                        innerReader.MoveToContent();
                        innerReader.ReadStartElement("selectedIon"); // Throws exception if we are not at the "selectedIon" tag.
                        while (innerReader.ReadState == ReadState.Interactive)
                        {
                            // Handle exiting out properly at EndElement tags
                            if (innerReader.NodeType != XmlNodeType.Element)
                            {
                                innerReader.Read();
                                continue;
                            }
                            switch (innerReader.Name)
                            {
                                case "referenceableParamGroupRef":
                                    // Schema requirements: zero to many instances of this element
                                    innerReader.Skip();
                                    break;
                                case "cvParam":
                                    // Schema requirements: zero to many instances of this element
                                    /* MUST supply a *child* term of MS:1000455 (ion selection attribute) one or more times
                                     *   e.g.: MS:1000041 (charge state)
                                     *   e.g.: MS:1000042 (intensity)
                                     *   e.g.: MS:1000633 (possible charge state)
                                     *   e.g.: MS:1000744 (selected ion m/z)
                                     */
                                    switch (innerReader.GetAttribute("accession"))
                                    {
                                        case "MS:1000041":
                                            // name="charge state"
                                            ion.Charge = (int)Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                        case "MS:1000744":
                                            // name="selected ion m/z"
                                            ion.SelectedIonMz = Convert.ToDouble(innerReader.GetAttribute("value"));
                                            break;
                                    }
                                    innerReader.Read(); // Consume the cvParam element (no child nodes)
                                    break;
                                case "userParam":
                                    // Schema requirements: zero to many instances of this element
                                    if (innerReader.GetAttribute("name") == "old charge state")
                                    {
                                        ion.OldCharge = (int)Convert.ToDouble(innerReader.GetAttribute("value"));
                                    }
                                    innerReader.Read();
                                    break;
                                default:
                                    innerReader.Skip();
                                    break;
                            }
                        }
                        innerReader.Dispose();
                        ions.Add(ion);
                        reader.Read(); // "selectedIon" might not have child nodes
                        // We will either consume the EndElement, or the same element that was passed to ReadSpectrum (in case of no child nodes)
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();
            return ions;
        }

        /// <summary>
        /// Handle a single binaryDataArrayList element and child nodes
        /// Called by ReadSpectrum (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single binaryDataArrayList element</param>
        /// <param name="defaultArrayLength">Default array length, coming from spectrum attribute</param>
        /// <returns>Enumerable list of binary data arrays</returns>
        private IEnumerable<BinaryDataArray> ReadBinaryDataArrayList(XmlReader reader, int defaultArrayLength)
        {
            reader.MoveToContent();
            // var binaryDataArrays = Convert.ToInt32(reader.GetAttribute("count"));
            var binaryData = new List<BinaryDataArray>();
            reader.ReadStartElement("binaryDataArrayList"); // Throws exception if we are not at the "binaryDataArrayList" tag.
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "binaryDataArray":
                        // Schema requirements: two to many instances of this element
                        binaryData.Add(ReadBinaryDataArray(reader.ReadSubtree(), defaultArrayLength));
                        reader.ReadEndElement(); // "SpectrumIdentificationItem" must have child nodes
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();
            return binaryData;
        }

        /// <summary>
        /// Handle a single binaryDataArray element and child nodes
        /// Called by ReadBinaryDataArrayList (XML hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single binaryDataArray element</param>
        /// <param name="defaultLength">Default array length, coming from spectrum attribute</param>
        /// <returns>Binary data array</returns>
        private BinaryDataArray ReadBinaryDataArray(XmlReader reader, int defaultLength)
        {
            reader.MoveToContent();
            var bda = new BinaryDataArray
            {
                ArrayLength = defaultLength
            };
            // var encLength = Convert.ToInt32(reader.GetAttribute("encodedLength"));
            var arrLength = Convert.ToInt32(reader.GetAttribute("arrayLength")); // Override the default; if non-existent, should get 0
            if (arrLength > 0)
            {
                bda.ArrayLength = arrLength;
            }
            var compressed = false;
            reader.ReadStartElement("binaryDataArray"); // Throws exception if we are not at the "spectrum" tag.
            var paramList = new List<Param>();
            while (reader.ReadState == ReadState.Interactive)
            {
                // Handle exiting out properly at EndElement tags
                if (reader.NodeType != XmlNodeType.Element)
                {
                    reader.Read();
                    continue;
                }

                switch (reader.Name)
                {
                    case "referenceableParamGroupRef":
                        // Schema requirements: zero to many instances of this element
                        var rpgRef = reader.GetAttribute("ref");

                        if (string.IsNullOrWhiteSpace(rpgRef))
                        {
                            ConsoleMsgUtils.ShowWarning("Encountered referenceableParamGroupRef with null or empty ref attribute");
                        }
                        else
                        {
                            paramList.AddRange(_referenceableParamGroups[rpgRef]);
                        }

                        reader.Read();
                        break;

                    case "cvParam":
                        // Schema requirements: zero to many instances of this element
                        paramList.Add(ReadCvParam(reader.ReadSubtree()));
                        reader.Read(); // Consume the cvParam element (no child nodes)
                        break;

                    case "userParam":
                        // Schema requirements: zero to many instances of this element
                        paramList.Add(ReadUserParam(reader.ReadSubtree()));
                        reader.Read();
                        break;

                    case "binary":
                        // Schema requirements: zero to many instances of this element
                        // Process the ParamList first.
                        foreach (var param in paramList)
                        {
                            /*
                             * MUST supply a *child* term of MS:1000572 (binary data compression type) only once
                             *   e.g.: MS:1000574 (zlib compression)
                             *   e.g.: MS:1000576 (no compression)
                             * MUST supply a *child* term of MS:1000513 (binary data array) only once
                             *   e.g.: MS:1000514 (m/z array)
                             *   e.g.: MS:1000515 (intensity array)
                             *   e.g.: MS:1000516 (charge array)
                             *   e.g.: MS:1000517 (signal to noise array)
                             *   e.g.: MS:1000595 (time array)
                             *   e.g.: MS:1000617 (wavelength array)
                             *   e.g.: MS:1000786 (non-standard data array)
                             *   e.g.: MS:1000820 (flow rate array)
                             *   e.g.: MS:1000821 (pressure array)
                             *   e.g.: MS:1000822 (temperature array)
                             * MUST supply a *child* term of MS:1000518 (binary data type) only once
                             *   e.g.: MS:1000521 (32-bit float)
                             *   e.g.: MS:1000523 (64-bit float)
                             */
                            switch (param.Accession)
                            {
                                // MUST supply a *child* term of MS:1000572 (binary data compression type) only once
                                case "MS:1000574":
                                    //   e.g.: MS:1000574 (zlib compression)
                                    compressed = true;
                                    break;
                                case "MS:1000576":
                                    //   e.g.: MS:1000576 (no compression)
                                    compressed = false;
                                    break;
                                // MUST supply a *child* term of MS:1000513 (binary data array) only once
                                case "MS:1000514":
                                    //   e.g.: MS:1000514 (m/z array)
                                    bda.ArrayType = ArrayType.m_z_array;
                                    break;
                                case "MS:1000515":
                                    //   e.g.: MS:1000515 (intensity array)
                                    bda.ArrayType = ArrayType.intensity_array;
                                    break;
                                case "MS:1000516":
                                    //   e.g.: MS:1000516 (charge array)
                                    bda.ArrayType = ArrayType.charge_array;
                                    break;
                                case "MS:1000517":
                                    //   e.g.: MS:1000517 (signal to noise array)
                                    bda.ArrayType = ArrayType.signal_to_noise_array;
                                    break;
                                case "MS:1000595":
                                    //   e.g.: MS:1000595 (time array)
                                    bda.ArrayType = ArrayType.time_array;
                                    break;
                                case "MS:1000617":
                                    //   e.g.: MS:1000617 (wavelength array)
                                    bda.ArrayType = ArrayType.wavelength_array;
                                    break;
                                case "MS:1000786":
                                    //   e.g.: MS:1000786 (non-standard data array)
                                    bda.ArrayType = ArrayType.non_standard_data_array;
                                    break;
                                case "MS:1000820":
                                    //   e.g.: MS:1000820 (flow rate array)
                                    bda.ArrayType = ArrayType.flow_rate_array;
                                    break;
                                case "MS:1000821":
                                    //   e.g.: MS:1000821 (pressure array)
                                    bda.ArrayType = ArrayType.pressure_array;
                                    break;
                                case "MS:1000822":
                                    //   e.g.: MS:1000822 (temperature array)
                                    bda.ArrayType = ArrayType.temperature_array;
                                    break;
                                // MUST supply a *child* term of MS:1000518 (binary data type) only once
                                case "MS:1000521":
                                    //   e.g.: MS:1000521 (32-bit float)
                                    bda.Precision = Precision.Precision32;
                                    break;
                                case "MS:1000523":
                                    //   e.g.: MS:1000523 (64-bit float)
                                    bda.Precision = Precision.Precision64;
                                    break;
                            }
                        }

                        var dataSize = 8;
                        if (bda.Precision == Precision.Precision32)
                        {
                            dataSize = 4;
                        }

                        var bytes = Convert.FromBase64String(reader.ReadElementContentAsString()); // Consumes the start and end elements.
                        //var bytesRead = reader.ReadContentAsBase64(bytes, 0, dataSize);
                        if (compressed)
                        {
                            bytes = DecompressZLib(bytes, bda.ArrayLength * dataSize);
                        }
                        if (bytes.Length % dataSize != 0 || bytes.Length / dataSize != bda.ArrayLength)
                        {
                            // We need to fail out.
                        }
                        //byte[] oneNumber = new byte[dataSize];
                        //bool swapBytes = true;
                        bda.Data = new double[bda.ArrayLength];
                        for (var i = 0; i < bytes.Length; i += dataSize)
                        {
                            // mzML binary data should always be Little Endian. Some other data formats may use Big Endian, which would require a byte swap

                            //Array.Copy(bytes, i, oneNumber, 0, dataSize);
                            //if (swapBytes)
                            //{
                            //  Array.Reverse(oneNumber);
                            //}

                            if (dataSize == 4)
                            {
                                //bda.Data[i / dataSize] = BitConverter.ToSingle(oneNumber, 0);
                                bda.Data[i / dataSize] = BitConverter.ToSingle(bytes, i);
                            }
                            // ReSharper disable once ConditionIsAlwaysTrueOrFalse
                            else if (dataSize == 8)
                            {
                                //bda.Data[i / dataSize] = BitConverter.ToDouble(oneNumber, 0);
                                bda.Data[i / dataSize] = BitConverter.ToDouble(bytes, i);
                            }
                        }
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Dispose();
            return bda;
        }

        /*********************************************************************************************************************************************
         * TODO: Flesh out the algorithm/double check it, etc.
         * Do some more work here.
         *
         ********************************************************************************************************************************************/
        private byte[] DecompressZLib(byte[] compressedBytes, int expectedBytes)
        {
            var msCompressed = new MemoryStream(compressedBytes);
            // We must skip the first two bytes
            // See http://george.chiramattel.com/blog/2007/09/deflatestream-block-length-does-not-match.html
            // Eat the zlib headers, the rest is a normal deflated stream
            msCompressed.ReadByte();
            msCompressed.ReadByte();
            //var msInflated = new MemoryStream((int)(msCompressed.Length * 2));
            //var newBytes = new byte[msCompressed.Length * 2];
            var newBytes = new byte[expectedBytes];
            // The last 32 bits (4 bytes) are supposed to be an Adler-32 checksum. Might need to remove them as well.
            using (var inflater = new DeflateStream(msCompressed, CompressionMode.Decompress))
            {
                var bytesRead = inflater.Read(newBytes, 0, expectedBytes);
                if (bytesRead != expectedBytes)
                {
                    throw new XmlException("Fail decompressing data...");
                }
                //while (inflater.CanRead)
                //{
                //  var readBytes = new byte[4095];
                //  // Should be able to change to just this.
                //  var bytesRead = inflater.Read(readBytes, 0, readBytes.Length);
                //  if (bytesRead != 0)
                //  {
                //      msInflated.Write(readBytes, 0, bytesRead);
                //  }
                //}
            }
            //newBytes = new byte[msInflated.Length];
            //msInflated.Read(newBytes, 0, (int)msInflated.Length);
            return newBytes;
        }

        #endregion
    }
}
