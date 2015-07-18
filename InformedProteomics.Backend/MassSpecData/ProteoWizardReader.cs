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
    /// This class currently doesn't work because of an unresolved dll dependency.
    /// </summary>
    public sealed class ProteoWizardReader: IMassSpecDataReader, IDisposable
    {
        public static Assembly ProteoWizardAssemblyResolver(object sender, ResolveEventArgs args)
        {
            Console.WriteLine("Searching for ProteoWizard files...");
            // https://support.microsoft.com/en-us/kb/837908
            //This handler is called only when the common language runtime tries to bind to the assembly and fails.
            if (string.IsNullOrWhiteSpace(PwizPath))
            {
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
            Console.WriteLine("Loading file \"" + strTempAssmbPath + "\"");
            //Load the assembly from the specified path.  
            Assembly myAssembly = null;
            try
            {
                myAssembly = Assembly.LoadFrom(strTempAssmbPath);
            }
            catch (BadImageFormatException ex)
            {
                Console.WriteLine("Incompatible Assembly: \"" + strTempAssmbPath + "\"");
                throw;
            }
            catch (FileNotFoundException ex)
            {
                Console.WriteLine("Assembly not found: \"" + strTempAssmbPath + "\"");
                throw;
            }
            catch (FileLoadException ex)
            {
                Console.WriteLine("Invalid Assembly: \"" + strTempAssmbPath + "\". The assembly may be marked as \"Untrusted\" by Windows. Please unblock and try again.");
                throw;
            }
            catch (SecurityException ex)
            {
                Console.WriteLine("Assembly access denied: \"" + strTempAssmbPath + "\"");
                throw;
            }

            //Return the loaded assembly.
            return myAssembly;
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
        public static string FindPwizPath()
        {
            string pwizPath = Environment.GetEnvironmentVariable("ProteoWizard");
            if (string.IsNullOrWhiteSpace(pwizPath) && Directory.Exists(@"C:\DMS_Programs\ProteoWizard"))
            {
                pwizPath = @"C:\DMS_Programs\ProteoWizard";
            }
            if (string.IsNullOrWhiteSpace(pwizPath))
            {
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

        static ProteoWizardReader()
        {
            PwizPath = FindPwizPath();
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="filePath"></param>
        /// <remarks>To avoid assembly resolving errors, the ProteoWizardAssemblyResolver should be added as an AssemblyResolve event handler, as follows:
        /// <code>AppDomain.CurrentDomain.AssemblyResolve += ProteoWizardReader.ProteoWizardAssemblyResolver;</code>
        /// </remarks>
        public ProteoWizardReader(string filePath)
        {
            AppDomain.CurrentDomain.AssemblyResolve += ProteoWizardAssemblyResolver;

            var readers = ReaderList.FullReaderList;
            filePath = CheckForDirectoryDataset(filePath);
            readers.read(filePath, _dataFile);
            if ((new string[] {".mzml", ".mzml.gz", ".mzxml", ".mzxml.gz", ".mgf", ".mgf.gz", ".txt"})
                .Any(ext => filePath.ToLower().EndsWith(ext)))
            {
                Console.WriteLine("Using cwt Centroiding");
                _filters.Add(_cwtCentroiding);
            }
            else
            {
                Console.WriteLine("Using vendor Centroiding");
                _filters.Add(_vendorCentroiding);
            }
            SpectrumListFactory.wrap(_dataFile, _filters);
        }

        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            for (int i = 1; i <= NumSpectra; i++)
            {
                yield return ReadMassSpectrum(i);
            }
        }

        public bool TryMakeRandomAccessCapable()
        {
            return true;
        }

        public Spectrum ReadMassSpectrum(int scanIndex, bool includePeaks = true)
        {
            var pwizSpec = _dataFile.run.spectrumList.spectrum(scanIndex - 1, includePeaks);

            var msLevel = (int)(pwizSpec.cvParam(CVID.MS_ms_level).value);
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
                    scanTime = (double) (s.cvParam(CVID.MS_scan_start_time).value);
                }
            }
            if (msLevel > 1)
            {
                double? thermoMonoMass = null;
                foreach (var up in pwizSpec.userParams)
                {
                    if (up.name == "[Thermo Trailer Extra]Monoisotopic M/Z:")
                    {
                        thermoMonoMass = (double) (up.value);
                    }
                }
                ActivationMethod am = ActivationMethod.Unknown;
                Data.Spectrometry.IsolationWindow iw = null;
                foreach (var precursor in pwizSpec.precursors)
                {
                    var act = precursor.activation;
                    foreach (var param in act.cvParams)
                    {
                        switch (param.cvid)
                        {
                            case CVID.MS_collision_induced_dissociation:
                                am = ActivationMethod.CID;
                                break;
                            case CVID.MS_electron_transfer_dissociation:
                                am = ActivationMethod.ETD;
                                break;
                            case CVID.MS_beam_type_collision_induced_dissociation:
                                am = ActivationMethod.HCD;
                                break;
                            case CVID.MS_electron_capture_dissociation:
                                am = ActivationMethod.ECD;
                                break;
                            case CVID.MS_pulsed_q_dissociation:
                                am = ActivationMethod.PQD;
                                break;
                        }
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
                        selectedIonMz = (double) (si.cvParam(CVID.MS_selected_ion_m_z).value);
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
                    ActivationMethod = am,
                    IsolationWindow = iw,
                    MsLevel = msLevel,
                    ElutionTime = scanTime,
                };
            }
            return new Spectrum(mzArray, intensityArray, scanIndex)
            {
                NativeId = pwizSpec.id,
                ElutionTime = scanTime,
            };
        }

        public void Close()
        {
            _dataFile.Dispose();
        }

        public int NumSpectra
        {
            get { return _dataFile.run.spectrumList.size(); }
        }

        private readonly MSData _dataFile = new MSData();
        private readonly string _vendorCentroiding = "peakPicking true 1-";
        private readonly string _cwtCentroiding = "peakPicking cwt snr=1.0 peakSpace=0.1 msLevel=1-";
        private readonly List<string> _filters = new List<string>();

        public void Dispose()
        {
            _dataFile.Dispose();
        }

        public static string ProteoWizardFilterString
        {
            get
            {
                return "All Supported|*.raw;*.mzML;*.mzML.gz;*.mzXML;*.mzXML.gz;*.mgf;*.mgf.gz;*.d;mspeak.bin;msprofile.bin;*.wiff;*.d;*.u2;FID;analysis.yep;analysis.baf;*.raw;_extern.inf;_inlet.inf;_FUNC*.DAT;*.lcd"
                    + "|Thermo .RAW|*.raw"
                    + "|mzML[.gz]|*.mzML;*.mzML.gz"
                    + "|mzXML[.gz]|*.mzXML;*.mzXML.gz"
                    + "|MGF[.gz]|*.mgf;*.mgf.gz"
                    + "|Agilent .d|*.d;mspeak.bin;msprofile.bin"
                    + "|AB Sciex .wiff|*.wiff"
                    + "|Bruker .d/FID/YEP/BAF|*.d;*.u2;FID;analysis.yep;analysis.baf"
                    + "|Waters .raw|*.raw;_extern.inf;_inlet.inf;_FUNC*.DAT"
                    + "|Shimadzu lcd|*.lcd";
            }
        }

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
    }
}
