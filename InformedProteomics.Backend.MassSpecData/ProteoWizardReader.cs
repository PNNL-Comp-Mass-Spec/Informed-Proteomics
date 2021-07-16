using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Security.Cryptography;
using InformedProteomics.Backend.Data.Spectrometry;
using pwiz.CLI.analysis;
using pwiz.CLI.cv;
using pwiz.CLI.msdata;
using Spectrum = InformedProteomics.Backend.Data.Spectrometry.Spectrum;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// ProteoWizard Reader, using the ProteoWizard pwiz_bindings_cli to utilize the ProteoWizard suite of vendor readers
    /// </summary>
    /// <remarks>
    /// <para>
    /// This class uses a custom AssemblyResolver to find an installation of ProteoWizard, specified in ProteoWizardReaderImplementation
    /// </para>
    /// <para>
    /// This class is a wrapper around ProteoWizardReaderImplementation to encapsulate the usage of the custom AssemblyResolver, which must be
    /// added to the AppDomain.CurrentDomain.AssemblyResolve event before the class is instantiated
    /// </para>
    /// </remarks>
    public sealed class ProteoWizardReader : IMassSpecDataReader
    {
        // Ignore Spelling: accessor, Apps, centroiding, cli, fid, lcd, pre, Pwiz, snr, stateful, uimf
        // Ignore Spelling: Bruker, SciEx, Shimadzu

        #region Static stateful data

        /// <summary>
        /// The path to the most recent 64-bit ProteoWizard install
        /// </summary>
        /// <remarks>
        /// If this is not null/empty, we can usually make a safe assumption that the ProteoWizard DLLs are available
        /// </remarks>
        public static string PwizPath => ProteoWizardLoader.PwizPath;

        /// <summary>
        /// Finds the path to the most recent 64-bit ProteoWizard install.
        /// PwizPath is populated from this, but only causes a single search.
        /// Paths searched, in order: "%ProteoWizard%" environment variable data, "C:\DMS_Programs\ProteoWizard", "%ProgramFiles%\ProteoWizard\(highest sorted)"
        /// </summary>
        /// <returns>ProteoWizard directory path, or null if not found</returns>
        public static string FindPwizPath()
        {
            return ProteoWizardLoader.FindPwizPath();
        }

        #endregion

        #region Constructor and Private member

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="filePath"></param>
        public ProteoWizardReader(string filePath)
        {
            ProteoWizardLoader.AddAssemblyResolver();
            ProteoWizardLoader.ValidateLoader();

            _pwizReader = new ProteoWizardReaderImplementation(filePath);
        }

        private readonly ProteoWizardReaderImplementation _pwizReader;

        #endregion

        #region IMassSpecDataReader implementation

        /// <summary>
        /// Gets all spectra
        /// </summary>
        /// <param name="includePeaks"></param>
        /// <returns>all spectra</returns>
        public IEnumerable<Spectrum> ReadAllSpectra(bool includePeaks = true)
        {
            return _pwizReader.ReadAllSpectra(includePeaks);
        }

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes
        /// </summary>
        /// <remarks>
        /// Child term of PSI-MS term MS:1000767, native spectrum identifier format
        /// </remarks>
        public PSI_Interface.CV.CV.CVID NativeIdFormat => _pwizReader.NativeIdFormat;

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes
        /// </summary>
        /// <remarks>
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </remarks>
        public PSI_Interface.CV.CV.CVID NativeFormat => _pwizReader.NativeFormat;

        /// <summary>
        /// Try to make the reader random access capable
        /// </summary>
        /// <returns>true if is random access capable, false if not</returns>
        public bool TryMakeRandomAccessCapable()
        {
            return _pwizReader.TryMakeRandomAccessCapable();
        }

        /// <summary>
        /// Returns the spectrum specified by the scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="includePeaks"></param>
        /// <returns>Spectrum object</returns>
        public Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true)
        {
            return _pwizReader.ReadMassSpectrum(scanNum, includePeaks);
        }

        /// <summary>
        /// Read the specified spectrum from the file, optionally reading only the metadata
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="includePeaks"></param>
        /// <returns>Spectrum object</returns>
        public Spectrum GetSpectrum(int scanNum, bool includePeaks = true)
        {
            return ReadMassSpectrum(scanNum, includePeaks);
        }

        /// <summary>
        /// Close the reader
        /// </summary>
        public void Close()
        {
            _pwizReader.Dispose();
        }

        /// <summary>
        /// The number of spectra in the file
        /// </summary>
        public int NumSpectra => _pwizReader.NumSpectra;

        /// <summary>
        /// Close the file
        /// </summary>
        public void Dispose()
        {
            _pwizReader.Dispose();
        }

        /// <summary>
        /// Path to the file; is <see cref="string.Empty"/> if the reader is in-memory
        /// </summary>
        public string FilePath => _pwizReader.FilePath;

        /// <summary>
        /// SHA-1 Checksum of the original input file (raw, mzML, .d folder, etc.)
        /// </summary>
        public string SrcFileChecksum => _pwizReader.SrcFileChecksum;

        /// <summary>
        /// Version of the immediate prior input file (raw, mzML, .d folder, etc.)
        /// </summary>
        public string FileFormatVersion => _pwizReader.FileFormatVersion;

        #endregion

        #region Public static utility functions (stateless)

        /// <summary>
        /// Filter string designed to be used in a file browser
        /// </summary>
        // ReSharper disable StringLiteralTypo
        public static string ProteoWizardFilterString => "All Supported|*.raw;*.mzML;*.mzML.gz;*.mzXML;*.mzXML.gz;*.mgf;*.mgf.gz;*.d;mspeak.bin;msprofile.bin;*.wiff;*.d;*.u2;FID;analysis.yep;analysis.baf;*.raw;_extern.inf;_inlet.inf;_FUNC*.DAT;*.lcd;*.uimf"
                                                         + "|Thermo .RAW|*.raw"
                                                         + "|mzML[.gz]|*.mzML;*.mzML.gz"
                                                         + "|mzXML[.gz]|*.mzXML;*.mzXML.gz"
                                                         + "|MGF[.gz]|*.mgf;*.mgf.gz"
                                                         + "|Agilent .d|*.d;mspeak.bin;msprofile.bin"
                                                         + "|AB Sciex .wiff|*.wiff"
                                                         + "|Bruker .d/FID/YEP/BAF|*.d;*.u2;FID;analysis.yep;analysis.baf"
                                                         + "|Waters .raw|*.raw;_extern.inf;_inlet.inf;_FUNC*.DAT"
                                                         + "|Shimadzu lcd|*.lcd"
                                                         + "|UIMF|*.uimf";

        // ReSharper restore StringLiteralTypo

        /// <summary>
        /// All files that ProteoWizard supports, either directly, or as part of a folder dataset
        /// </summary>
        public static List<string> SupportedFilesFilterList => new List<string>()
        {
            // ReSharper disable StringLiteralTypo
            ".raw", // Thermo (file) or Waters (folder)
            "_extern.inf", // Waters .raw folder content
            "_inlet.inf", // Waters .raw folder content
            "_FUNC*.DAT", // Waters .raw folder content
            ".d", // Agilent (folder) or Bruker (folder)
            "mspeak.bin", // Agilent .d folder content
            "msprofile.bin", // Agilent .d folder content
            ".yep", // Bruker
            ".baf", // Bruker
            "fid", // Bruker
            ".lcd", // Shimadzu
            ".wiff", // Waters
            ".mzml",
            ".mzml.gz",
            ".mzxml",
            ".mzxml.gz",
            ".mgf",
            ".mgf.gz",
            ".uimf",
            // ReSharper restore StringLiteralTypo
        };

        /// <summary>
        /// List of "folder extensions" that ProteoWizard can read
        /// </summary>
        /// <remarks>
        /// This does not include all folder type datasets - some require directory listings
        /// </remarks>
        public static List<string> SupportedDirectoryTypes => new List<string> { ".d", ".raw" };

        /// <summary>
        /// List of files that are produced by Bruker instruments that ProteoWizard can read
        /// </summary>
        public static List<string> BrukerFiles => new List<string>()
        {
            ".d",
            "analysis.yep",
            "analysis.baf",
            "fid",
        };

        /// <summary>
        /// All file extensions that ProteoWizard directly reads - i.e., we don't need to back out of a folder
        /// </summary>
        private static List<string> DirectlySupportedFilesFilterList => new List<string>()
        {
            // ReSharper disable StringLiteralTypo
            ".raw", // Thermo (file) or Waters (folder)
            ".d", // Agilent (folder) or Bruker (folder)
            ".yep", // Bruker
            ".baf", // Bruker
            "fid", // Bruker
            ".lcd", // Shimadzu
            ".wiff", // Waters
            ".mzml",
            ".mzml.gz",
            ".mzxml",
            ".mzxml.gz",
            ".mgf",
            ".mgf.gz",
            ".uimf",
            // ReSharper restore StringLiteralTypo
        };

        /// <summary>
        /// Check the file path to see if it is to files in a directory dataset type (.raw folder, or .d folder)
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns>Path to directory, if a directory, otherwise returns filePath</returns>
        public static string CheckForDirectoryDataset(string filePath)
        {
            var fullFilePath = Path.GetFullPath(filePath);
            var tempFilePath = fullFilePath.ToLower();

            // If it doesn't end with .raw, but contains .raw, and the path is not an existing file or a supported ProteoWizard file type,
            // assume that it is a Waters .raw directory.
            if (!tempFilePath.EndsWith(".raw") && tempFilePath.Contains(".raw") &&
                (!File.Exists(fullFilePath) ||
                 !DirectlySupportedFilesFilterList.Any(f => tempFilePath.EndsWith(f))))
            {
                while (!string.IsNullOrWhiteSpace(tempFilePath) && tempFilePath.Contains(".raw"))
                {
                    if (tempFilePath.EndsWith(".raw") &&
                        Directory.Exists(fullFilePath.Substring(0, tempFilePath.Length)))
                    {
                        return fullFilePath.Substring(0, tempFilePath.Length);
                    }
                    tempFilePath = Path.GetDirectoryName(tempFilePath);
                }
            }

            tempFilePath = fullFilePath.ToLower();
            // If it doesn't end with .d, but contains .d, and the path is not an existing file or it is a Bruker data type
            // Assume it is a Agilent or Bruker .d dataset, and change the path to only point to the directory.
            if (tempFilePath.EndsWith(".d") && !Directory.Exists(fullFilePath) || tempFilePath.Contains(".d") &&
                (!File.Exists(fullFilePath) || BrukerFiles.Any(f => tempFilePath.EndsWith(f)) || tempFilePath.EndsWith(".bin")))
            {
                while (!string.IsNullOrWhiteSpace(tempFilePath) && tempFilePath.Contains(".d"))
                {
                    if (tempFilePath.EndsWith(".d") &&
                        Directory.Exists(fullFilePath.Substring(0, tempFilePath.Length)))
                    {
                        return fullFilePath.Substring(0, tempFilePath.Length);
                    }
                    tempFilePath = Path.GetDirectoryName(tempFilePath);
                }
            }
            return filePath;
        }

        #endregion
    }

    /// <summary>
    /// ProteoWizardReaderImplementation, using the ProteoWizard pwiz_bindings_cli to utilize the ProteoWizard suite of vendor readers
    /// </summary>
    /// <remarks>
    /// <para>
    /// This class uses a custom AssemblyResolver in class ProteoWizardLoader to find an installation of ProteoWizard
    /// </para>
    /// <para>
    /// If there are DLL resolving problems when trying to use it, add the ProteoWizardAssemblyResolver to the
    /// AppDomain.CurrentDomain.AssemblyResolve event before first instantiating the class
    /// </para>
    /// </remarks>
    internal sealed class ProteoWizardReaderImplementation : IMassSpecDataReader
    {
        #region Constructor

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="filePath"></param>
        /// <remarks>
        /// To avoid assembly resolving errors, the ProteoWizardAssemblyResolver should be added as an AssemblyResolve event handler, as follows:
        /// AppDomain.CurrentDomain.AssemblyResolve += ProteoWizardLoader.ProteoWizardAssemblyResolver;
        /// </remarks>
        public ProteoWizardReaderImplementation(string filePath)
        {
            AppDomain.CurrentDomain.AssemblyResolve += ProteoWizardLoader.ProteoWizardAssemblyResolver;

            FilePath = ProteoWizardReader.CheckForDirectoryDataset(filePath);
        }

        #endregion

        #region Private members and functions

        private bool _loaded;
        private int _numSpectra;

        private readonly MSData _dataFile = new MSData();
        // Uses the centroiding/peak picking algorithm that the vendor libraries provide, if available; otherwise uses a low-quality centroiding algorithm
        private readonly string _vendorCentroiding = "peakPicking true 1-";
        // Continuous Wavelet Transform peak picker - high-quality peak picking, may be slow with some high-res data.
        private readonly string _cwtCentroiding = "peakPicking cwt snr=1.0 peakSpace=0.1 msLevel=1-";
        private readonly List<string> _filters = new List<string>();
        private string _fileFormatVersion = string.Empty;
        private string _srcFileChecksum = string.Empty;
        private MethodInfo _binaryDataArrayGetData;

        private void LoadPwizReader()
        {
            if (_loaded)
            {
                return;
            }

            var readers = ReaderList.FullReaderList;
            readers.read(FilePath, _dataFile);
            // ReSharper disable StringLiteralTypo
            if (new[] { ".mzml", ".mzml.gz", ".mzxml", ".mzxml.gz", ".mgf", ".mgf.gz", ".txt", "uimf", "uimf.gz" }
                .Any(ext => FilePath.ToLower().EndsWith(ext)))
            {
                // Files that do not have vendor centroiding available
                Console.WriteLine("Using cwt Centroiding");
                _filters.Add(_cwtCentroiding);
            }
            // ReSharper restore StringLiteralTypo
            else
            {
                Console.WriteLine("Using vendor Centroiding");
                _filters.Add(_vendorCentroiding);
            }
            SpectrumListFactory.wrap(_dataFile, _filters);

            _numSpectra = _dataFile.run.spectrumList.size();
            foreach (var srcFile in _dataFile.fileDescription.sourceFiles)
            {
                var cv = srcFile.cvParam(CVID.MS_SHA_1);
                if (!cv.empty())
                {
                    _srcFileChecksum = cv.value;
                    break;
                }
            }

            if (string.IsNullOrWhiteSpace(_srcFileChecksum))
            {
                string fileForChecksum;
                if (File.Exists(FilePath))
                {
                    // is file
                    fileForChecksum = FilePath;
                }
                else
                {
                    var dir = new DirectoryInfo(FilePath);
                    var files = dir.GetFiles("*", SearchOption.AllDirectories);
                    if (files.Length == 0)
                    {
                        throw new IOException("Cannot open directory with no files: " + FilePath);
                    }
                    fileForChecksum = files[0].FullName;
                }
                using (var fs = new FileStream(fileForChecksum, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
                using (var sha1 = new SHA1Managed())
                {
                    var hash = sha1.ComputeHash(fs);
                    _srcFileChecksum = BitConverter.ToString(hash).ToLower().Replace("-", string.Empty);
                }
            }

            foreach (var software in _dataFile.softwareList)
            {
                _fileFormatVersion = software.version;
            }

            // BinaryDataArray.get_data() problem fix
            // Pre-Nov. 7th, 2018 pwiz_binding_cli.dll: bda.data returns pwiz.CLI.msdata.BinaryData,  is a semi-automatic wrapper for a C++ vector, which implements IList<double>
            // Pre-Nov. 7th, 2018 pwiz_binding_cli.dll: bda.data returns pwiz.CLI.util.BinaryData implements IList<double>, but also provides other optimization functions
            // The best way to access this before was bda.data.ToArray()
            // In the future, this could be changed to bda.data.Storage.ToArray(), but that may lead to more data copying than just using the IEnumerable<double> interface
            // Both versions implement IList<double>, so I can get the object via reflection and cast it to an IList<double> (or IEnumerable<double>).

            // Get the MethodInfo for BinaryDataArray.data property accessor
            _binaryDataArrayGetData = typeof(BinaryDataArray).GetProperty("data", BindingFlags.Public | BindingFlags.Instance | BindingFlags.IgnoreCase)?.GetMethod;

            _loaded = true;
        }

        #endregion

        #region IMassSpecDataReader implementation

        /// <summary>
        /// Gets all spectra
        /// </summary>
        /// <param name="includePeaks"></param>
        /// <returns>all spectra</returns>
        public IEnumerable<Spectrum> ReadAllSpectra(bool includePeaks = true)
        {
            LoadPwizReader();
            for (var i = 1; i <= _numSpectra; i++)
            {
                yield return ReadSpectrum(i, includePeaks);
            }
        }

        /// <summary>
        /// The number of spectra in the file
        /// </summary>
        public int NumSpectra
        {
            get
            {
                LoadPwizReader();
                return _numSpectra;
            }
        }

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes
        /// </summary>
        /// <remarks>
        /// Child term of PSI-MS term MS:1000767, native spectrum identifier format
        /// </remarks>
        public PSI_Interface.CV.CV.CVID NativeIdFormat
        {
            get
            {
                LoadPwizReader();
                foreach (var file in _dataFile.fileDescription.sourceFiles)
                {
                    foreach (var cvParam in file.cvParams)
                    {
                        if (CV.cvIsA(cvParam.cvid, CVID.MS_nativeID_format))
                        {
                            var cvInt = (int)cvParam.cvid;
                            var format = (PSI_Interface.CV.CV.CVID)cvInt;
                            return format;
                        }
                    }
                }
                return PSI_Interface.CV.CV.CVID.CVID_Unknown;
            }
        }

        /// <summary>
        /// The Native Format of the source file - needed for tracking purposes
        /// </summary>
        /// <remarks>
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </remarks>
        public PSI_Interface.CV.CV.CVID NativeFormat
        {
            get
            {
                LoadPwizReader();
                foreach (var file in _dataFile.fileDescription.sourceFiles)
                {
                    foreach (var cvParam in file.cvParams)
                    {
                        if (CV.cvIsA(cvParam.cvid, CVID.MS_mass_spectrometer_file_format))
                        {
                            var cvInt = (int)cvParam.cvid;
                            var format = (PSI_Interface.CV.CV.CVID)cvInt;
                            return format;
                        }
                    }
                }
                return PSI_Interface.CV.CV.CVID.CVID_Unknown;
            }
        }

        /// <summary>
        /// Try to make the reader random access capable
        /// </summary>
        /// <returns>true if is random access capable, false if not</returns>
        public bool TryMakeRandomAccessCapable()
        {
            return true;
        }

        /// <summary>
        /// Returns the spectrum specified by the scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="includePeaks"></param>
        /// <returns>Spectrum object</returns>
        public Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true)
        {
            LoadPwizReader();
            return ReadSpectrum(scanNum, includePeaks);
        }

        /// <summary>
        /// Read the specified spectrum from the file, optionally reading only the metadata
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="includePeaks"></param>
        /// <returns>Spectrum object</returns>
        public Spectrum GetSpectrum(int scanNum, bool includePeaks = true)
        {
            return ReadMassSpectrum(scanNum, includePeaks);
        }

        /// <summary>
        /// This method exists to handle an abstraction needed for now to handle pwiz_bindings_cli.dll changes committed to ProteoWizard on November 7, 2018
        /// This abstraction may be removed at a point in the future when versions of ProteoWizard older than that date are less likely to be in use.
        /// </summary>
        /// <param name="bda"></param>
        /// <returns>Array of binary data</returns>
        private double[] GetBinaryDataAsArray(BinaryDataArray bda)
        {
            // BinaryDataArray.get_data() problem fix
            // Pre-Nov. 7th, 2018 pwiz_binding_cli.dll: bda.data returns pwiz.CLI.msdata.BinaryData,  is a semi-automatic wrapper for a C++ vector, which implements IList<double>
            // Pre-Nov. 7th, 2018 pwiz_binding_cli.dll: bda.data returns pwiz.CLI.util.BinaryData implements IList<double>, but also provides other optimization functions
            // The best way to access this before was bda.data.ToArray()
            // In the future, this could be changed to bda.data.Storage.ToArray(), but that may lead to more data copying than just using the IEnumerable<double> interface
            // Both versions implement IList<double>, so I can get the object via reflection and cast it to an IList<double> (or IEnumerable<double>).

            // Call via reflection to avoid issues of the InformedProteomics compiled reference vs. the ProteoWizard compiled DLL
            var dataObj = _binaryDataArrayGetData?.Invoke(bda, null);
            if (dataObj is IEnumerable<double> data)
            {
                return data.ToArray();
            }

            return new double[0];
        }

        /// <summary>
        /// Internal spectrum reader to eliminate excess calls to LoadPwizReader() when called from ReadAllSpectra()
        /// </summary>
        /// <param name="scanIndex"></param>
        /// <param name="includePeaks"></param>
        /// <returns>Spectrum object</returns>
        private Spectrum ReadSpectrum(int scanIndex, bool includePeaks = true)
        {
            var pwizSpec = _dataFile.run.spectrumList.spectrum(scanIndex - 1, includePeaks);

            var msLevel = (int)(pwizSpec.cvParam(CVID.MS_ms_level).value);
            var tic = (double)(pwizSpec.cvParam(CVID.MS_total_ion_current).value);
            var mzArray = new double[0];
            var intensityArray = new double[0];
            foreach (var bda in pwizSpec.binaryDataArrays)
            {
                if (bda.hasCVParam(CVID.MS_m_z_array))
                {
                    mzArray = GetBinaryDataAsArray(bda);
                }
                if (bda.hasCVParam(CVID.MS_intensity_array))
                {
                    intensityArray = GetBinaryDataAsArray(bda);
                }
            }
            double scanTime = 0;
            double driftTime = 0;
            foreach (var s in pwizSpec.scanList.scans)
            {
                if (s.hasCVParam(CVID.MS_scan_start_time))
                {
                    var timeCvParam = s.cvParam(CVID.MS_scan_start_time);
                    scanTime = (double)timeCvParam.value;
                    // Conversion: CV dictates possible units of 'minute' and 'second'
                    if (timeCvParam.units == CVID.UO_second)
                    {
                        scanTime /= 60;
                    }
                }

                if (s.hasCVParam(CVID.MS_ion_mobility_drift_time))
                {
                    var driftCvParam = s.cvParam(CVID.MS_ion_mobility_drift_time);
                    // No conversion: CV dictates only valid units of 'millisecond'
                    driftTime = (double)driftCvParam.value;
                }
            }
            if (msLevel > 1)
            {
                double? thermoMonoMass = null;
                foreach (var up in pwizSpec.userParams)
                {
                    if (up.name == "[Thermo Trailer Extra]Monoisotopic M/Z:")
                    {
                        thermoMonoMass = (double)(up.value);
                    }
                }
                var am = ActivationMethod.Unknown;
                Data.Spectrometry.IsolationWindow iw = null;
                foreach (var precursor in pwizSpec.precursors)
                {
                    var act = precursor.activation;
                    var activationMethods = new List<ActivationMethod>();
                    foreach (var param in act.cvParams)
                    {
                        switch (param.cvid)
                        {
                            case CVID.MS_collision_induced_dissociation:
                                activationMethods.Add(ActivationMethod.CID);
                                break;
                            case CVID.MS_electron_transfer_dissociation:
                                activationMethods.Add(ActivationMethod.ETD);
                                break;
                            case CVID.MS_beam_type_collision_induced_dissociation:
                                activationMethods.Add(ActivationMethod.HCD);
                                break;
                            case CVID.MS_electron_capture_dissociation:
                                activationMethods.Add(ActivationMethod.ECD);
                                break;
                            case CVID.MS_pulsed_q_dissociation:
                                activationMethods.Add(ActivationMethod.PQD);
                                break;
                        }
                    }
                    if (activationMethods.Count > 1 && activationMethods.Contains(ActivationMethod.ETD))
                    {
                        am = ActivationMethod.ETD;
                    }
                    else if (activationMethods.Count > 0)
                    {
                        am = activationMethods[0];
                    }

                    var piw = precursor.isolationWindow;
                    var target = (double)(piw.cvParam(CVID.MS_isolation_window_target_m_z).value);
                    var lowOff = (double)(piw.cvParam(CVID.MS_isolation_window_lower_offset).value);
                    var uppOff = (double)(piw.cvParam(CVID.MS_isolation_window_upper_offset).value);
                    int? charge = null;
                    double selectedIonMz = 0;
                    foreach (var si in precursor.selectedIons)
                    {
                        if (si.hasCVParam(CVID.MS_charge_state))
                        {
                            charge = (int)(si.cvParam(CVID.MS_charge_state).value);
                        }
                        selectedIonMz = (double)(si.cvParam(CVID.MS_selected_ion_m_z).value);
                    }
                    if (thermoMonoMass == null || thermoMonoMass.Value.Equals(0))
                    {
                        thermoMonoMass = selectedIonMz;
                    }
                    if (target.Equals(0))
                    {
                        target = selectedIonMz;
                    }
                    iw = new Data.Spectrometry.IsolationWindow(target, lowOff, uppOff, thermoMonoMass, charge);
                }
                return new ProductSpectrum(mzArray, intensityArray, scanIndex)
                {
                    NativeId = pwizSpec.id,
                    TotalIonCurrent = tic,
                    ActivationMethod = am,
                    IsolationWindow = iw,
                    MsLevel = msLevel,
                    ElutionTime = scanTime,
                    DriftTime = driftTime,
                };
            }
            return new Spectrum(mzArray, intensityArray, scanIndex)
            {
                NativeId = pwizSpec.id,
                TotalIonCurrent = tic,
                ElutionTime = scanTime,
                MsLevel = msLevel,
                DriftTime = driftTime,
            };
        }

        /// <summary>
        /// Close the file
        /// </summary>
        public void Close()
        {
            _dataFile.Dispose();
        }

        /// <summary>
        /// Close the file
        /// </summary>
        public void Dispose()
        {
            _dataFile.Dispose();
        }

        /// <summary>
        /// Path to the file; is <see cref="string.Empty"/> if the reader is in-memory
        /// </summary>
        public string FilePath { get; }

        /// <summary>
        /// SHA-1 Checksum of the original input file (raw, mzML, .d folder, etc.)
        /// </summary>
        public string SrcFileChecksum
        {
            get
            {
                LoadPwizReader();
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
                LoadPwizReader();
                return _fileFormatVersion;
            }
        }

        #endregion
    }
}
