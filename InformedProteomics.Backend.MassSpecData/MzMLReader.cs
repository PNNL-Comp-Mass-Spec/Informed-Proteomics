using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text;
using System.Xml;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;
using PSI_Interface.CV;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Reader for mzML files. Can handle gzipped mzML files, and read in a forward-only fashion or in a random-access fashion.
    /// </summary>
    public sealed class MzMLReader: IMassSpecDataReader
    {
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
        private long _artificialScanNum = 1;
        private long _numSpectra = -1;
        private readonly IndexList _spectrumOffsets = new IndexList() {IndexType = IndexList.IndexListType.Spectrum};
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
        /// Store the mzML version, so that we can use it to adjust how some things are processed.
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

        private abstract class Param
        {
            protected ParamType ParamType { get; set; }

            public string Name;          // Required
            public string Value;         // Optional
            public string UnitCVRef;     // Optional
            public string UnitAccession; // Optional
            public string UnitName;      // Optional

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

        private class CVParam : Param
        {
            public override string CVRef      // Required
            {
                get => _cvRef;
                set => _cvRef = value;
            }

            public override string Accession  // Required
            {
                get => _accession;
                set => _accession = value;
            }

            public CVParam()
            {
                ParamType = ParamType.cvParam;
            }
        }

        private class UserParam : Param
        {
            public override string Type       // Optional
            {
                get => _type;
                set => _type = value;
            }

            public UserParam()
            {
                ParamType = ParamType.userParam;
            }
        }

        private class IndexList
        {
            private long _artificialScanNum = 1;
            public class IndexItem // A struct would be faster, but it can also be a pain since it is a value type
            {
                public readonly string Ref;
                public readonly long Offset;
                public readonly long IdNum;

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
                _offsets.Clear();
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
            private readonly List<IndexItem> _offsets = new List<IndexItem>();

            public List<IndexItem> Offsets => _offsets;

            // Unused
            // private readonly Dictionary<string, long> OffsetsMapNative = new Dictionary<string, long>();

            public readonly Dictionary<long, long> OffsetsMapInt = new Dictionary<long, long>();

            // Unused
            // private readonly Dictionary<long, string> IdToNativeMap = new Dictionary<long, string>();

            public readonly Dictionary<string, long> NativeToIdMap = new Dictionary<string, long>();

            public void AddOffset(string idRef, string offset)
            {
                AddOffset(idRef, Int64.Parse(offset));
            }

            public void AddOffset(string idRef, long offset)
            {
                var scanNum = _artificialScanNum++;

                if (NativeIdConversion.TryGetScanNumberLong(idRef, out var num))
                {
                    scanNum = num;
                }

                var item = new IndexItem(idRef, offset, scanNum);
                AddMapForOffset(item);
                _offsets.Add(item);
            }

            private void AddMapForOffset(IndexItem offset)
            {
                if (IndexType == IndexListType.Chromatogram)
                {
                    return;
                }

                // OffsetsMapNative.Add(offset.Ref, offset.Offset);

                // This won't be sufficient until there is a valid parser for all forms of NativeID.
                // Using artificial scan number for now.
                OffsetsMapInt.Add(offset.IdNum, offset.Offset);
                // IdToNativeMap.Add(offset.IdNum, offset.Ref);
                NativeToIdMap.Add(offset.Ref, offset.IdNum);
                /*if (IndexType == IndexListType.Spectrum)
                {
                    long id = Int64.Parse(offset.Ref.Substring(offset.Ref.LastIndexOfAny(new char[] {'=', 'F'}) + 1));
                    OffsetsMapInt.Add(id, offset.Offset);
                    IdNativeMap.Add(id, offset.Ref);
                    offset.IdNum = id;
                }*/
            }

            public void RegenerateMaps()
            {
                // OffsetsMapNative.Clear();
                OffsetsMapInt.Clear();
                // IdToNativeMap.Clear();
                NativeToIdMap.Clear();
                foreach (var offset in _offsets)
                {
                    AddMapForOffset(offset);
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
                var tokens = nativeId.Split(new[] {'\t', ' ', '\n', '\r'}, StringSplitOptions.RemoveEmptyEntries);
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
            /// <returns></returns>
            public static bool TryGetScanNumberLong(string nativeId, out long num)
            {
                return long.TryParse(GetScanNumber(nativeId), out num);
            }

            /// <summary>
            /// Try to get the scan number for <paramref name="nativeId"/>, and parse it as an int
            /// </summary>
            /// <param name="nativeId"></param>
            /// <param name="num"></param>
            /// <returns></returns>
            public static bool TryGetScanNumberInt(string nativeId, out int num)
            {
                return int.TryParse(GetScanNumber(nativeId), out num);
            }

            /// <summary>
            /// For the supplied <paramref name="nativeId"/>, get the corresponding scan number
            /// </summary>
            /// <param name="nativeId"></param>
            /// <returns></returns>
            /// <remarks>Code is ported from MSData.cpp in ProteoWizard</remarks>
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
                //            return "";
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
                //        return "";
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
            //            return "";
            //    }
            //}
        }

        private readonly Dictionary<string, List<Param>> _referenceableParamGroups = new Dictionary<string, List<Param>>();

        private class SelectedIon
        {
            public double SelectedIonMz;
            public int Charge;

            // ReSharper disable once NotAccessedField.Local
            public int OldCharge;

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
        /// Initialize a MzMlReader object
        /// </summary>
        /// <param name="filePath">Path to mzML file</param>
        /// <param name="randomAccess">If mzML reader should be configured for random access</param>
        /// <param name="tryReducingMemoryUsage">If mzML reader should try to avoid reading all spectra into memory. This will reduce memory usage for a non-random access MzMLReader, as long as ReadMassSpectrum(int) isn't used.</param>
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
            _file?.Close();

            // Set a very large read buffer, it does decrease the read times for uncompressed files.
            _file = new FileStream(_filePath, FileMode.Open, FileAccess.Read, FileShare.Read, 65536);

            /*****************************************************************************************************************************************************
             * TODO: Change how the file handles are used for safety purposes - open up each time, or what?
             *****************************************************************************************************************************************************/

            if (_filePath.EndsWith(".mzML.gz"))
            {
                _isGzipped = true;
                var file = new GZipStream(_file, CompressionMode.Decompress);
                if (!_randomAccess)
                {
                    _file = file;
                }
                else
                {
                    // Unzip the file to the temp path
                    _unzippedFilePath = Path.Combine(Path.GetTempPath(), Path.GetFileNameWithoutExtension(_filePath));
                    using (_file)
                    using (file)
                    using (
                        var tempFile = new FileStream(_unzippedFilePath, FileMode.Create, FileAccess.ReadWrite,
                            FileShare.None, 65536))
                    {
                        file.CopyTo(tempFile/*, 65536*/);
                    }
                    _file = new FileStream(_unzippedFilePath, FileMode.Open, FileAccess.Read, FileShare.Read, 65536);
                }
            }
            _fileReader = new StreamReader(_file, Encoding.UTF8, true, 65536);

            if (!_isGzipped || _randomAccess) // can't reset the position on a gzipped file...
            {
                // perform a read to perform encoding autodetection
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
        /// The number of spectra in the file.
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
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000767, native spectrum identifier format
        /// </summary>
        public CV.CVID NativeIdFormat
        {
            get
            {
                RequireMetadata();
                return _nativeIdFormat;
            }
        }

        /// <summary>
        /// The Native Format of the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </summary>
        public CV.CVID NativeFormat
        {
            get
            {
                RequireMetadata();
                return _nativeFormat;
            }
        }

        /// <summary>
        /// Try to make the reader random access capable
        /// </summary>
        /// <returns>true if is random access capable, false if not</returns>
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
        /// Uses "yield return" to allow processing one spectra at a time if called from a foreach loop statement.
        /// </summary>
        /// <returns></returns>
        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            if (!_randomAccess)
            {
                return ReadAllSpectraNonRandom();
            }
            return ReadAllSpectraRandom();
        }

        /// <summary>
        /// Returns a single spectrum from the file
        /// </summary>
        /// <param name="index"></param>
        /// <param name="includePeaks"></param>
        /// <returns></returns>
        /// <remarks>If random access mode is turned on, this will respond quickly and use only as much memory as is needed to store the spectrum.
        /// If random access mode is off, this will cause the memory usage reducing mode to shut of, and all spectra will be read into memory.</remarks>
        public Spectrum ReadMassSpectrum(int index, bool includePeaks = true)
        {
            // Proper functionality when not random access
            if (!_randomAccess)
            {
                return ReadMassSpectrumNonRandom(index);
            }
            return ReadMassSpectrumRandom(index, includePeaks);
        }
        #endregion

        #region Interface Helper Functions: Non-Random Access
        /// <summary>
        /// Read all mass spectra in the file, not using random access
        /// Uses "yield return" to use less memory when called from a "foreach" statement
        /// </summary>
        /// <returns></returns>
        private IEnumerable<Spectrum> ReadAllSpectraNonRandom()
        {
            if (_reduceMemoryUsage)
            {
                _artificialScanNum = 1;
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
                        yield return ReadSpectrum(_xmlReaderForYield.ReadSubtree(), true);
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
        /// Read a single mass spectrum and return it.
        /// Causes all spectra in the file to be loaded into memory
        /// </summary>
        /// <param name="index"></param>
        private Spectrum ReadMassSpectrumNonRandom(int index)
        {
            if (!_allRead)
            {
                _artificialScanNum = 1;
                _reduceMemoryUsage = false; // They called this on a non-random access reader, now they suffer the consequences.
                ReadMzMl();
            }
            return _spectra[(int)index];
        }
        #endregion

        #region Interface Helper Functions: Random Access
        /// <summary>
        /// Read all mass spectra in the file, using random access
        /// Uses "yield return" to use less memory when called from a "foreach" statement
        /// </summary>
        /// <returns></returns>
        private IEnumerable<Spectrum> ReadAllSpectraRandom()
        {
            if (!_haveIndex || !_haveMetaData)
            {
                ReadMzMl(); // Read the index and metadata so that the offsets get populated.
            }
            foreach (var specIndex in _spectrumOffsets.Offsets)
            {
                yield return ReadMassSpectrumRandom(specIndex.IdNum);
            }
        }

        /// <summary>
        /// Read a single mass spectrum and return it.
        /// </summary>
        /// <param name="index"></param>
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
            _xmlReaderForYield?.Close();
            _fileReader?.Close();
            _file?.Close();
        }

        /// <summary>
        /// Delete unzipped file, if we had to unzip the file to read it.
        /// </summary>
        public void Cleanup()
        {
            if (_randomAccess && _isGzipped)
            {
                File.Delete(_unzippedFilePath);
            }
        }

        /// <summary>
        /// Close and cleanup file handles
        /// </summary>
        public void Dispose()
        {
            Close();
            Cleanup();
        }

        /// <summary>
        /// Close and cleanup file handles
        /// </summary>
        ~MzMLReader()
        {
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
        /// Find and read the index information, starting at the end of the file...
        /// </summary>
        private void ReadIndexFromEnd()
        {
            var stream = new FileStream(_unzippedFilePath, FileMode.Open, FileAccess.Read, FileShare.Read, 1);
            var testPos = stream.Length;
            //stream.Position = testPos; // 300 bytes from the end of the file - should be enough
            var streamReader = new StreamReader(stream, System.Text.Encoding.UTF8, true, 65536);
            streamReader.DiscardBufferedData();
            var haveOffset = false;

            while (!haveOffset)
            {
                // "<indexListOffset"
                const int bufSize = 512; //65536 (17 bits), 131072 (18 bits), 262144 (19 bits), 524288 (20 bits)
                testPos -= bufSize;
                stream.Position = testPos;
                var byteBuffer = new byte[bufSize];
                var stringBuffer = string.Empty;
                var bufStart = stream.Position;
                var bytesRead = 0;
                while (stream.Position < stream.Length && !haveOffset)
                {
                    bufStart = stream.Position;
                    bytesRead = stream.Read(byteBuffer, 0, bufSize);
                    stringBuffer = _encoding.GetString(byteBuffer, 0, bytesRead);
                    // set up the rewind to ensure full tags
                    var lastTagEnd = stringBuffer.LastIndexOf('>');
                    var lastTagStart = stringBuffer.LastIndexOf('<');
                    if (lastTagStart != -1 && lastTagEnd != -1 && lastTagStart > lastTagEnd)
                    {
                        var endOfString = lastTagStart;
                        var rewindBy = _encoding.GetByteCount(stringBuffer.Substring(endOfString));
                        stringBuffer = stringBuffer.Substring(0, endOfString);
                        stream.Seek(-rewindBy, SeekOrigin.Current);
                        //file.Position = bufEnd - rewindBy;
                    }

                    var found = stringBuffer.IndexOf("<indexListOffset");
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
                            reader2.Close();
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
                    var stringBuffer = string.Empty;
                    var bufStart = stream.Position;
                    var bytesRead = 0;
                    while (stream.Position < stream.Length && !haveOffset)
                    {
                        bufStart = stream.Position;
                        bytesRead = stream.Read(byteBuffer, 0, bufSize);
                        stringBuffer = _encoding.GetString(byteBuffer, 0, bytesRead);
                        // set up the rewind to ensure full tags
                        var lastTagEnd = stringBuffer.LastIndexOf('>');
                        var lastTagStart = stringBuffer.LastIndexOf('<');
                        if (lastTagStart != -1 && lastTagEnd != -1 && lastTagStart > lastTagEnd)
                        {
                            var endOfString = lastTagStart;
                            var rewindBy = _encoding.GetByteCount(stringBuffer.Substring(endOfString));
                            stringBuffer = stringBuffer.Substring(0, endOfString);
                            stream.Seek(-rewindBy, SeekOrigin.Current);
                            //file.Position = bufEnd - rewindBy;
                        }

                        var found = stringBuffer.IndexOf("<indexList ");
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
        /// Read the Checksum from the indexedmzML data
        /// </summary>
        private void ReadChecksum()
        {
            using (var stream = new FileStream(_unzippedFilePath, FileMode.Open, FileAccess.Read, FileShare.Read, 1))
            using (var streamReader = new StreamReader(stream, System.Text.Encoding.UTF8, true, 65536))
            {
                try
                {
                    var testPos = stream.Length - 500;
                    stream.Position = testPos;
                    streamReader.DiscardBufferedData();
                    var data = streamReader.ReadToEnd();
                    var pos = data.IndexOf("<fileChecksum", StringComparison.InvariantCultureIgnoreCase);
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
        /// Called by IndexMzMl (xml hierarchy)
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
                var stringBuffer = string.Empty;
                var bufStart = file.Position;
                var bytesRead = 0;
                var builder = string.Empty;
                while (file.Position < file.Length)
                {
                    bufStart = file.Position;
                    bytesRead = file.Read(byteBuffer, 0, maxRead);
                    stringBuffer = _encoding.GetString(byteBuffer, 0, bytesRead);
                    // set up the rewind to ensure full tags
                    var lastTagEnd = stringBuffer.LastIndexOf('>');
                    var lastTagStart = stringBuffer.LastIndexOf('<');
                    if (lastTagStart != -1 && lastTagEnd != -1 && lastTagStart > lastTagEnd)
                    {
                        var endOfString = lastTagStart;
                        var rewindBy = _encoding.GetByteCount(stringBuffer.Substring(endOfString));
                        stringBuffer = stringBuffer.Substring(0, endOfString);
                        file.Seek(-rewindBy, SeekOrigin.Current);
                        //file.Position = bufEnd - rewindBy;
                    }

                    var searchPoint = 0;
                    while (searchPoint < stringBuffer.Length)
                    {
                        var foundSpec = stringBuffer.IndexOf("<" + specTag + " ", searchPoint);
                        var foundChrom = stringBuffer.IndexOf("<" + chromTag + " ", searchPoint);
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
                        builder = stringBuffer.Substring(searchPoint + 1, (end - 1) - (searchPoint + 1));
                        // Get the ID of the tag
                        var attribName = "id";
                        if (_version == MzML_Version.mzML1_0_0)
                        {
                            attribName = "nativeID";
                        }
                        var idIndex = builder.IndexOf(attribName + "=\"");
                        var idOpenQuote = builder.IndexOf("\"", idIndex);
                        var idCloseQuote = builder.IndexOf("\"", idOpenQuote + 1);
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
                            pos = file.Position - 1; // position of caret
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
                        _spectrumOffsets.AddOffset(id, pos);
                    }
                    else if (builder.StartsWith(chromTag))
                    {
                        _chromatogramOffsets.AddOffset(id, pos);
                    }
                }*/
                _haveIndex = true;
            }
        }

        /// <summary>
        /// Handle the child nodes of the indexedmzML element
        /// Called by IndexMzMl (xml hierarchy)
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
            reader.Close();
        }

        /// <summary>
        /// Handle the child nodes of the indexList element
        /// Called by ReadIndexList (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "index" element</param>
        private void ReadIndex(XmlReader reader)
        {
            reader.MoveToContent();
            var iType = reader.GetAttribute("name");
            var eType = IndexList.IndexListType.Unknown;
            if (iType.ToLower() == "spectrum")
            {
                eType = IndexList.IndexListType.Spectrum;
            }
            else if (iType.ToLower() == "chromatogram")
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
            reader.Close();
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
                // Read to the mzML root tag, and ignore the extra indexedmzML data
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
                    _srcFileChecksum = BitConverter.ToString(hash).ToLower().Replace("-", "");
                }
            }
            var schemaName = reader.GetAttribute("xsi:schemaLocation");
            // We automatically assume it uses the mzML_1.1.0 schema. Check for the old version.
            //if (!schemaName.Contains("mzML1.1.0.xsd"))
            if (schemaName.Contains("mzML1.0.0.xsd"))
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
                reader.Close();
            } */
            if (!_reduceMemoryUsage)
            {
                // Don't worry about closing the subtree readers, just close the root reader.
                // reader is the root if it is not an indexed mzML file.
                if (indexReader == null)
                {
                    reader.Close();
                }
                else
                {
                    indexReader.Close();
                }
            }
        }
        #endregion

        #region Metadata tag readers
        /// <summary>
        /// Handle the child nodes of the fileDescription element
        /// Called by ReadMzML (xml hierarchy)
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
            reader.Close();
        }

        /// <summary>
        /// Handle a single sourceFileList element and child nodes
        /// Called by ReadMzML (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single sourceFileList element</param>
        private void ReadSourceFileList(XmlReader reader)
        {
            reader.MoveToContent();
            //int count = Convert.ToInt32(reader.GetAttribute("count"));
            reader.ReadStartElement("sourceFileList"); // Throws exception if we are not at the "sourceFileList" tag.
            //while (reader.ReadState == ReadState.Interactive && _instrument == Instrument.Unknown)
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
                                 *   et al.
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
                                 *   et al.
                                 */
                                var cv = innerReader.GetAttribute("cvRef");
                                var accession = innerReader.GetAttribute("accession");
                                if (string.IsNullOrWhiteSpace(cv))
                                {
                                    cv = "MS";
                                }
                                if (cv.ToUpper().Contains("PSI"))
                                {
                                    cv = "MS";
                                }
                                if (CV.TermAccessionLookup[cv].ContainsKey(accession))
                                {
                                    var cvid = CV.TermAccessionLookup[cv][accession];
                                    if (CV.CvidIsA(cvid, CV.CVID.MS_native_spectrum_identifier_format))
                                    {
                                        _nativeIdFormat = cvid;
                                    }
                                    else if (CV.CvidIsA(cvid, CV.CVID.MS_mass_spectrometer_file_format))
                                    {
                                        _nativeFormat = cvid;
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
                    innerReader.Close();
                    reader.Read();
                }
                else
                {
                    reader.Read();
                }
            }
            reader.Close();
        }

        /// <summary>
        /// Handle the child nodes of the referenceableParamGroupList element
        /// Called by ReadMzML (xml hierarchy)
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
                    innerReader.Close();
                    reader.Read();
                    _referenceableParamGroups.Add(id, paramList);
                }
                else
                {
                    reader.Read();
                }
            }
            reader.Close();
        }

        /// <summary>
        /// Handle the cvParam element
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "cvParam" element</param>
        private CVParam ReadCvParam(XmlReader reader)
        {
            reader.MoveToContent();
            var cvParam = new CVParam();
            cvParam.Accession = reader.GetAttribute("accession");
            cvParam.CVRef = reader.GetAttribute("cvRef");
            cvParam.Name = reader.GetAttribute("name");
            cvParam.Value = reader.GetAttribute("value");
            cvParam.UnitAccession = reader.GetAttribute("unitAccession");
            cvParam.UnitCVRef = reader.GetAttribute("unitCVRef");
            cvParam.UnitName = reader.GetAttribute("unitName");

            reader.Close();
            return cvParam;
        }

        /// <summary>
        /// Handle the userParam element
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single "userParam" element</param>
        private UserParam ReadUserParam(XmlReader reader)
        {
            reader.MoveToContent();
            var userParam = new UserParam();
            userParam.Name = reader.GetAttribute("name");
            userParam.Type = reader.GetAttribute("type");
            userParam.Value = reader.GetAttribute("value");
            userParam.UnitAccession = reader.GetAttribute("unitAccession");
            userParam.UnitCVRef = reader.GetAttribute("unitCVRef");
            userParam.UnitName = reader.GetAttribute("unitName");

            reader.Close();
            return userParam;
        }
        #endregion

        #region Run and SpectrumList Tags
        /// <summary>
        /// Handle the child nodes of the run element
        /// Called by ReadMzML (xml hierarchy)
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
            reader.Close();
        }

        /// <summary>
        /// Handle the child nodes of a spectrumList element
        /// Called by ReadRunData (xml hierarchy)
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
                //reader.Close(); // Closing can be slow for a subtree...
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
                    _spectra.Add(ReadSpectrum(reader.ReadSubtree(), true));
                    // "spectrum" might not have any child nodes
                    // We will either consume the EndElement, or the same element that was passed to ReadSpectrum (in case of no child nodes)
                    reader.Read();
                }
                else
                {
                    reader.Skip();
                }
            }
            reader.Close();
        }
        #endregion

        #region Spectrum Tag
        /// <summary>
        /// Handle a single spectrum element and child nodes
        /// Called by ReadSpectrumList (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single spectrum element</param>
        /// <param name="includePeaks">Whether to read binary data arrays</param>
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

            var scanNum = -1;
            // If a random access reader, there is already a scan number stored, based on the order of the index. Use it instead.
            if (_randomAccess)
            {
                scanNum = (int) (_spectrumOffsets.NativeToIdMap[nativeId]);
            }
            else
            {
                scanNum = (int)(_artificialScanNum++);
                // Interpret the NativeID (if the format has an interpreter) and use it instead of the artificial number.
                // TODO: Better handling than the artificial ID for other nativeIDs (ones currently not supported)
                if (NativeIdConversion.TryGetScanNumberInt(nativeId, out var num))
                {
                    scanNum = num;
                }
            }

            var defaultArraySize = Convert.ToInt32(reader.GetAttribute("defaultArrayLength"));
            reader.ReadStartElement("spectrum"); // Throws exception if we are not at the "spectrum" tag.
            var is_ms_ms = false;
            var msLevel = 0;
            var centroided = false;
            double tic = 0;
            var precursors = new List<Precursor>();
            var scans = new List<ScanData>();
            var bdas = new List<BinaryDataArray>();
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
                         *   et al.
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
                         *   et al.
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
                            bdas.AddRange(ReadBinaryDataArrayList(reader.ReadSubtree(), defaultArraySize));
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
            reader.Close();
            // Process the spectrum data
            var scan = new ScanData();
            Spectrum spectrum;
            var mzs = new BinaryDataArray();
            var intensities = new BinaryDataArray();
            foreach (var bda in bdas)
            {
                if (bda.ArrayType == ArrayType.m_z_array)
                {
                    mzs = bda;
                }
                else if (bda.ArrayType == ArrayType.intensity_array)
                {
                    intensities = bda;
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
            if (scans.Count == 1)
            {
                scan = scans[0];
            }
            else if (scans.Count > 1)
            {
                // TODO: Should do something else to appropriately handle combinations...
                scan = scans[0];
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

                var pspectrum = new ProductSpectrum(mzs.Data, intensities.Data, scanNum);
                pspectrum.ActivationMethod = precursor.Activation;
                // Select mz value to use based on presence of a Thermo-specific user param.
                // The user param has a slightly higher precision, if that matters.
                var mz = scan.MonoisotopicMz == 0.0 ? ion.SelectedIonMz : scan.MonoisotopicMz;
                pspectrum.IsolationWindow = new IsolationWindow(precursor.IsolationWindowTargetMz, precursor.IsolationWindowLowerOffset, precursor.IsolationWindowUpperOffset, mz, ion.Charge);
                //pspectrum.IsolationWindow.OldCharge = ion.OldCharge;
                //pspectrum.IsolationWindow.SelectedIonMz = ion.SelectedIonMz;
                spectrum = pspectrum;
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
        /// Called by ReadSpectrumList (xml hierarchy)
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
                        reader.ReadEndElement(); // "acquisitionList" mustt have child nodes
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
            reader.Close();
        }

        /// <summary>
        /// Handle a single scanList element and child nodes
        /// Called by ReadSpectrum (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single scanList element</param>
        /// <returns></returns>
        private List<ScanData> ReadScanList(XmlReader reader)
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
            reader.Close();
            return scans;
        }

        /// <summary>
        /// Handle a single scan element and child nodes
        /// Called by ReadSpectrum (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single scan element</param>
        /// <returns></returns>
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
            reader.Close();
            return scan;
        }

        /// <summary>
        /// Handle a single precursorList element and child nodes
        /// Called by ReadSpectrum (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single precursorList element</param>
        /// <returns></returns>
        private List<Precursor> ReadPrecursorList(XmlReader reader)
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
            reader.Close();
            return precursors;
        }

        /// <summary>
        /// Handle a single precursor element and child nodes
        /// Called by ReadPrecursorList (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single precursor element</param>
        /// <returns></returns>
        private Precursor ReadPrecursor(XmlReader reader)
        {
            reader.MoveToContent();
            reader.ReadStartElement("precursor"); // Throws exception if we are not at the "precursor" tag.
            XmlReader innerReader;
            var precursor = new Precursor();
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
                        innerReader.Close();
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
                                     *   et al.
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
                        innerReader.Close();
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
            reader.Close();
            return precursor;
        }

        /// <summary>
        /// Handle a single selectedIonList element and child nodes
        /// Called by ReadPrecursor (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single selectedIonList element</param>
        /// <returns></returns>
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
                                        ion.OldCharge = (int) Convert.ToDouble(innerReader.GetAttribute("value"));
                                    }
                                    innerReader.Read();
                                    break;
                                default:
                                    innerReader.Skip();
                                    break;
                            }
                        }
                        innerReader.Close();
                        ions.Add(ion);
                        reader.Read(); // "selectedIon" might not have child nodes
                        // We will either consume the EndElement, or the same element that was passed to ReadSpectrum (in case of no child nodes)
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
            return ions;
        }

        /// <summary>
        /// Handle a single binaryDataArrayList element and child nodes
        /// Called by ReadSpectrum (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single binaryDataArrayList element</param>
        /// <param name="defaultArrayLength">Default array length, coming from spectrum attribute</param>
        /// <returns></returns>
        private List<BinaryDataArray> ReadBinaryDataArrayList(XmlReader reader, int defaultArrayLength)
        {
            reader.MoveToContent();
            // var bdArrays = Convert.ToInt32(reader.GetAttribute("count"));
            var bdaList = new List<BinaryDataArray>();
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
                        bdaList.Add(ReadBinaryDataArray(reader.ReadSubtree(), defaultArrayLength));
                        reader.ReadEndElement(); // "SpectrumIdentificationItem" must have child nodes
                        break;
                    default:
                        reader.Skip();
                        break;
                }
            }
            reader.Close();
            return bdaList;
        }

        /// <summary>
        /// Handle a single binaryDataArray element and child nodes
        /// Called by ReadBinaryDataArrayList (xml hierarchy)
        /// </summary>
        /// <param name="reader">XmlReader that is only valid for the scope of the single binaryDataArray element</param>
        /// <param name="defaultLength">Default array length, coming from spectrum attribute</param>
        /// <returns></returns>
        private BinaryDataArray ReadBinaryDataArray(XmlReader reader, int defaultLength)
        {
            reader.MoveToContent();
            var bda = new BinaryDataArray();
            bda.ArrayLength = defaultLength;
            var encLength = Convert.ToInt32(reader.GetAttribute("encodedLength"));
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
                        paramList.AddRange(_referenceableParamGroups[rpgRef]);
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
                        //var bytesread = reader.ReadContentAsBase64(bytes, 0, dataSize);
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
            reader.Close();
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
            // EAT the zlib headers, the rest is a normal 'deflate'd stream
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
