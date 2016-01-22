using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Security;
using InformedProteomics.Backend.Data.Spectrometry;
using pwiz.CLI.analysis;
using pwiz.CLI.cv;
using pwiz.CLI.msdata;
using Spectrum = InformedProteomics.Backend.Data.Spectrometry.Spectrum;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// ProteoWizard Reader, using the ProteoWizard pwiz_bindings_cli to utilize the ProteoWizard suite of vendor readers.
    /// </summary>
    /// <remarks>This class uses a custom AssemblyResolver to find an installation of ProteoWizard, specified in ProteoWizardReaderImplementation.
    /// This class is a wrapper around ProteoWizardReaderImplementation to encapsulate the usage of the custom AssemblyResolver, which must be 
    /// added to the AppDomain.CurrentDomain.AssemblyResolve event before the class is instantiated.</remarks>
    public sealed class ProteoWizardReader: IMassSpecDataReader, IDisposable
    {
        #region Static stateful data

        /// <summary>
        /// The path to the most recent 64-bit ProteoWizard install
        /// If this is not null/empty, we can usually make a safe assumption that the ProteoWizard dlls are available.
        /// </summary>
        public static string PwizPath
        {
            get { return ProteoWizardReaderImplementation.PwizPath; }
        }

        /// <summary>
        /// Finds the path to the most recent 64-bit ProteoWizard install
        /// PwizPath is populated from this, but only causes a single search.
        /// Paths searched, in order: "%ProteoWizard%" environment variable data, "C:\DMS_Programs\ProteoWizard", "%ProgramFiles%\ProteoWizard\(highest sorted)"
        /// </summary>
        /// <returns></returns>
        public static string FindPwizPath()
        {
            return ProteoWizardReaderImplementation.FindPwizPath();
        }

        #endregion

        #region Constructor and Private member

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="filePath"></param>
        public ProteoWizardReader(string filePath)
        {
            ProteoWizardReaderImplementation.AddAssemblyResolver();
            ProteoWizardReaderImplementation.ValidateLoader();

            _pwizReader = new ProteoWizardReaderImplementation(filePath);
        }

        private readonly ProteoWizardReaderImplementation _pwizReader;

        #endregion

        #region IMassSpecDataReader implementation

        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            return _pwizReader.ReadAllSpectra();
        }

        public bool TryMakeRandomAccessCapable()
        {
            return _pwizReader.TryMakeRandomAccessCapable();
        }

        public Spectrum ReadMassSpectrum(int scanIndex, bool includePeaks = true)
        {
            return _pwizReader.ReadMassSpectrum(scanIndex, includePeaks);
        }

        public void Close()
        {
            _pwizReader.Dispose();
        }

        public int NumSpectra
        {
            get { return _pwizReader.NumSpectra; }
        }

        public void Dispose()
        {
            _pwizReader.Dispose();
        }

        #endregion

        #region Public static utility functions (stateless)

        /// <summary>
        /// Filter string designed to be used in a file browser
        /// </summary>
        public static string ProteoWizardFilterString
        {
            get
            {
                return "All Supported|*.raw;*.mzML;*.mzML.gz;*.mzXML;*.mzXML.gz;*.mgf;*.mgf.gz;*.d;mspeak.bin;msprofile.bin;*.wiff;*.d;*.u2;FID;analysis.yep;analysis.baf;*.raw;_extern.inf;_inlet.inf;_FUNC*.DAT;*.lcd;*.uimf"
                    + "|Thermo .RAW|*.raw"
                    + "|mzML[.gz]|*.mzML;*.mzML.gz"
                    + "|mzXML[.gz]|*.mzXML;*.mzXML.gz"
                    + "|MGF[.gz]|*.mgf;*.mgf.gz"
                    + "|Agilent .d|*.d;mspeak.bin;msprofile.bin"
                    + "|AB Sciex .wiff|*.wiff"
                    + "|Bruker .d/FID/YEP/BAF|*.d;*.u2;FID;analysis.yep;analysis.baf"
                    + "|Waters .raw|*.raw;_extern.inf;_inlet.inf;_FUNC*.DAT"
                    + "|Shimadzu lcd|*.lcd"
                    + "|UIMF|*.uimf"
                    ;
            }
        }

        /// <summary>
        /// All files that ProteoWizard supports, either directly, or as part of a folder dataset
        /// </summary>
        public static List<string> SupportedFilesFilterList
        {
            get
            {
                return new List<string>()
                {
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
                };
            }
        }

        public static List<string> SupportedDirectoryTypes
        {
            get
            {
                return new List<string>() {".d", ".raw"};
            }
        }

        public static List<string> BrukerFiles
        {
            get
            {
                return new List<string>()
                {
                    ".d",
                    "analysis.yep",
                    "analysis.baf",
                    "fid",
                };
            }
        }

        /// <summary>
        /// All file extensions that ProteoWizard directly reads - i.e., we don't need to back out of a folder
        /// </summary>
        private static List<string> DirectlySupportedFilesFilterList
        {
            get
            {
                return new List<string>()
                {
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
                };
            }
        }

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
            if ((tempFilePath.EndsWith(".d") && !Directory.Exists(fullFilePath)) || tempFilePath.Contains(".d") &&
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
    /// ProteoWizardReaderImplementation, using the ProteoWizard pwiz_bindings_cli to utilize the ProteoWizard suite of vendor readers.
    /// </summary>
    /// <remarks>This class uses a custom AssemblyResolver to find an installation of ProteoWizard. If there are DLL resolving
    /// problems when trying to use it, add the ProteoWizardAssemblyResolver to the AppDomain.CurrentDomain.AssemblyResolve
    /// event before first instantiating the class.</remarks>
    internal sealed class ProteoWizardReaderImplementation : IMassSpecDataReader, IDisposable
    {
        #region AssemblyResolverHandler for finding ProteoWizard dlls

        /// <summary>
        /// Add the Assembly Resolver to the system assembly resolver chain
        /// </summary>
        /// <remarks>This should be called early in the program, so that the ProteoWizard Assembly Resolver will 
        /// already be in the resolver chain before any other use of ProteoWizardWrapper.
        /// Also, DependencyLoader.ValidateLoader() should be used to make sure a meaningful error message is thrown if ProteoWizard is not available.</remarks>
        public static void AddAssemblyResolver()
        {
            if (!_resolverAdded)
            {
#if DEBUG
                Console.WriteLine("Adding assembly resolver...");
#endif
                AppDomain.CurrentDomain.AssemblyResolve += ProteoWizardReaderImplementation.ProteoWizardAssemblyResolver;
                _resolverAdded = true;
            }
        }

        private static bool _resolverAdded = false;

        /// <summary>
        /// On a missing DLL event, searches a path specified by FindPwizPath for the ProteoWizard dlls, and loads them
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="args"></param>
        /// <returns></returns>
        public static Assembly ProteoWizardAssemblyResolver(object sender, ResolveEventArgs args)
        {
#if DEBUG
            Console.WriteLine("Looking for: " + args.Name);
            //Console.WriteLine("Wanted by: " + args.RequestingAssembly);
#endif
            if (!args.Name.ToLower().StartsWith("pwiz_bindings_cli"))
            {
                return Assembly.LoadFrom(""); // We are not interested in searching for anything else - resolving pwiz_bindings_cli provides the hint for all of its dependencies.
                // This will actually trigger an exception, which is handled in the system code, and the dll search goes on down the chain.
                // returning null results in this code being called multiple times, for the same dependency.
            }
            Console.WriteLine("Searching for ProteoWizard files...");

            // https://support.microsoft.com/en-us/kb/837908
            //This handler is called only when the common language runtime tries to bind to the assembly and fails.
            if (string.IsNullOrWhiteSpace(PwizPath))
            {
                ValidateLoaderByPath();
                return null;
            }

            //Retrieve the list of referenced assemblies in an array of AssemblyName.
            string strTempAssmbPath = "";

            AssemblyName[] arrReferencedAssmbNames = Assembly.GetExecutingAssembly().GetReferencedAssemblies();

            //Loop through the array of referenced assembly names.
            foreach (AssemblyName strAssmbName in arrReferencedAssmbNames)
            {
                //Check for the assembly names that have raised the "AssemblyResolve" event.
                if (strAssmbName.FullName.Substring(0, strAssmbName.FullName.IndexOf(",")) == args.Name.Substring(0, args.Name.IndexOf(",")))
                {
                    //Console.WriteLine("Attempting to load DLL \"" + Path.Combine(pwizPath, args.Name.Substring(0, args.Name.IndexOf(",")) + ".dll") + "\"");
                    //Build the path of the assembly from where it has to be loaded.                
                    strTempAssmbPath = Path.Combine(PwizPath, args.Name.Substring(0, args.Name.IndexOf(",")) + ".dll");
                    break;
                }
            }
#if DEBUG
                Console.WriteLine("Loading file \"" + strTempAssmbPath + "\"");
#endif
            //Load the assembly from the specified path.  
            Assembly myAssembly = null;
            try
            {
                myAssembly = Assembly.LoadFrom(strTempAssmbPath);
            }
            catch (BadImageFormatException)
            {
                Console.WriteLine("Incompatible Assembly: \"" + strTempAssmbPath + "\"");
                throw;
            }
            catch (FileNotFoundException)
            {
                Console.WriteLine("Assembly not found: \"" + strTempAssmbPath + "\"");
                throw;
            }
            catch (FileLoadException)
            {
                Console.WriteLine("Invalid Assembly: \"" + strTempAssmbPath + "\". The assembly may be marked as \"Untrusted\" by Windows. Please unblock and try again.");
                throw;
            }
            catch (SecurityException)
            {
                Console.WriteLine("Assembly access denied: \"" + strTempAssmbPath + "\"");
                throw;
            }

            //Return the loaded assembly.
            return myAssembly;
        }

        #endregion

        #region Static stateful variable and populating functions

        /// <summary>
        /// Checks to make sure the path to ProteoWizard files is set. If not, throws an exception.
        /// </summary>
        /// <remarks>This function should generally only be called inside of a conditional statement to prevent the 
        /// exception from being thrown when the ProteoWizard dlls will not be needed.</remarks>
        public static void ValidateLoader()
        {
            try
            {
                Assembly.Load("pwiz_bindings_cli, Version=0.0.0.0, Culture=neutral, PublicKeyToken=null");
            }
            catch
            {
                var bits = Environment.Is64BitProcess ? "64" : "32";
                var message = CannotFindExceptionMessage();

                System.Console.WriteLine(message);
                throw new System.TypeLoadException(message);
            }
        }

        private static void ValidateLoaderByPath()
        {
            if (string.IsNullOrWhiteSpace(PwizPath))
            {
                var message = CannotFindExceptionMessage();

                System.Console.WriteLine(message);
                throw new System.TypeLoadException(message);
            }
        }

        private static string CannotFindExceptionMessage()
        {
            var bits = Environment.Is64BitProcess ? "64" : "32";
            var message = "Cannot load ProteoWizard dlls. Please ensure that " + bits
                + "-bit ProteoWizard is installed to its default install directory (\""
                + Environment.GetEnvironmentVariable("ProgramFiles") + "\\ProteoWizard\\ProteoWizard 3.0.[x]\").";

            return message;
        }

        /// <summary>
        /// The path to the most recent 64-bit ProteoWizard install
        /// If this is not null/empty, we can usually make a safe assumption that the ProteoWizard dlls are available.
        /// </summary>
        public static readonly string PwizPath;

        /// <summary>
        /// Finds the path to the most recent 64-bit ProteoWizard install
        /// PwizPath is populated from this, but only causes a single search.
        /// </summary>
        /// <returns></returns>
        /// <remarks>Paths searched, in order: 
        /// "%ProteoWizard%"/"%ProteoWizard%_x86" environment variable data, 
        /// "C:\DMS_Programs\ProteoWizard"/"C:\DMS_Programs\ProteoWizard_x86", 
        /// "%ProgramFiles%\ProteoWizard\(highest sorted)"</remarks>
        public static string FindPwizPath()
        {
            string pwizPath = string.Empty;

            // Set the DMS_Programs ProteoWizard path based on if the process is 32- or 64-bit.
            var dmsProgPwiz = @"C:\DMS_Programs\ProteoWizard";

            // Check for a x64 ProteoWizard environment variable
            pwizPath = Environment.GetEnvironmentVariable("ProteoWizard");

            if (!Environment.Is64BitProcess && Environment.Is64BitOperatingSystem)
            {
                // Check for a x86 ProteoWizard environment variable
                pwizPath = Environment.GetEnvironmentVariable("ProteoWizard_x86");
                dmsProgPwiz = @"C:\DMS_Programs\ProteoWizard_x86";
            }

            if (string.IsNullOrWhiteSpace(pwizPath) && Directory.Exists(dmsProgPwiz))
            {
                pwizPath = dmsProgPwiz;
            }
            if (string.IsNullOrWhiteSpace(pwizPath))
            {
                // NOTE: Should automatically function as-is to get 32-bit ProteoWizard for 32-bit process and 64-bit ProteoWizard for 64-bit process...
                var progFiles = Environment.GetEnvironmentVariable("ProgramFiles");
                if (string.IsNullOrWhiteSpace(progFiles))
                {
                    return null;
                }
                var progPwiz = Path.Combine(progFiles, "ProteoWizard");
                if (!Directory.Exists(progPwiz))
                {
                    return null;
                }
                var posPaths = Directory.GetDirectories(progPwiz, "ProteoWizard *");
                pwizPath = posPaths.Max(); // Try to get the "newest" folder
            }
            return pwizPath;
        }

        static ProteoWizardReaderImplementation()
        {
            PwizPath = FindPwizPath();
        }

        #endregion

        #region Constructor

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="filePath"></param>
        /// <remarks>To avoid assembly resolving errors, the ProteoWizardAssemblyResolver should be added as an AssemblyResolve event handler, as follows:
        /// <code>AppDomain.CurrentDomain.AssemblyResolve += ProteoWizardReader.ProteoWizardAssemblyResolver;</code>
        /// </remarks>
        public ProteoWizardReaderImplementation(string filePath)
        {
            AppDomain.CurrentDomain.AssemblyResolve += ProteoWizardAssemblyResolver;

            _filePath = ProteoWizardReader.CheckForDirectoryDataset(filePath);
        }

        #endregion

        #region Private members and functions
        
        private readonly string _filePath;
        private bool _loaded = false;
        private int _numSpectra = 0;

        private readonly MSData _dataFile = new MSData();
        // Uses the centroiding/peak picking algorithm that the vendor libraries provide, if available; otherwise uses a low-quality centroiding algorithm
        private readonly string _vendorCentroiding = "peakPicking true 1-";
        // Continuous Wavelet Transform peak picker - high-quality peak picking, may be slow with some high-res data.
        private readonly string _cwtCentroiding = "peakPicking cwt snr=1.0 peakSpace=0.1 msLevel=1-";
        private readonly List<string> _filters = new List<string>();
        private void LoadPwizReader()
        {
            if (_loaded)
            {
                return;
            }

            var readers = ReaderList.FullReaderList;
            readers.read(_filePath, _dataFile);
            if ((new string[] { ".mzml", ".mzml.gz", ".mzxml", ".mzxml.gz", ".mgf", ".mgf.gz", ".txt", "uimf", "uimf.gz" })
                .Any(ext => _filePath.ToLower().EndsWith(ext)))
            {
                // Files that do not have vendor centroiding available
                Console.WriteLine("Using cwt Centroiding");
                _filters.Add(_cwtCentroiding);
            }
            else
            {
                Console.WriteLine("Using vendor Centroiding");
                _filters.Add(_vendorCentroiding);
            }
            SpectrumListFactory.wrap(_dataFile, _filters);

            _numSpectra = _dataFile.run.spectrumList.size();
            _loaded = true;
        }

        #endregion

        #region IMassSpecDataReader implementation

        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            LoadPwizReader();
            for (int i = 1; i <= _numSpectra; i++)
            {
                yield return ReadSpectrum(i);
            }
        }

        public int NumSpectra
        {
            get
            {
                LoadPwizReader();
                return _numSpectra;
            }
        }

        public bool TryMakeRandomAccessCapable()
        {
            return true;
        }

        public Spectrum ReadMassSpectrum(int scanIndex, bool includePeaks = true)
        {
            LoadPwizReader();
            return ReadSpectrum(scanIndex, includePeaks);
        }

        /// <summary>
        /// Internal spectrum reader to eliminate excess calls to LoadPwizReader() when called from ReadAllSpectra()
        /// </summary>
        /// <param name="scanIndex"></param>
        /// <param name="includePeaks"></param>
        /// <returns></returns>
        private Spectrum ReadSpectrum(int scanIndex, bool includePeaks = true)
        {
            var pwizSpec = _dataFile.run.spectrumList.spectrum(scanIndex - 1, includePeaks);

            var msLevel = (int)(pwizSpec.cvParam(CVID.MS_ms_level).value);
            var tic = (double)(pwizSpec.cvParam(CVID.MS_total_ion_current).value);
            double[] mzArray = new double[0];
            double[] intensityArray = new double[0];
            foreach (var bda in pwizSpec.binaryDataArrays)
            {
                if (bda.hasCVParam(CVID.MS_m_z_array))
                {
                    mzArray = bda.data.ToArray();
                }
                if (bda.hasCVParam(CVID.MS_intensity_array))
                {
                    intensityArray = bda.data.ToArray();
                }
            }
            double scanTime = 0;
            foreach (var s in pwizSpec.scanList.scans)
            {
                if (s.hasCVParam(CVID.MS_scan_start_time))
                {
                    scanTime = (double)(s.cvParam(CVID.MS_scan_start_time).value);
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
                ActivationMethod am = ActivationMethod.Unknown;
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
                };
            }
            return new Spectrum(mzArray, intensityArray, scanIndex)
            {
                NativeId = pwizSpec.id,
                TotalIonCurrent = tic,
                ElutionTime = scanTime,
                MsLevel = msLevel,
            };
        }

        public void Close()
        {
            _dataFile.Dispose();
        }

        public void Dispose()
        {
            _dataFile.Dispose();
        }

        #endregion
    }
}
