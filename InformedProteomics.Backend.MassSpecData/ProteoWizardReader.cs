using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Security;
using System.Security.Cryptography;
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
    public sealed class ProteoWizardReader : IMassSpecDataReader
    {
        #region Static stateful data

        /// <summary>
        /// The path to the most recent 64-bit ProteoWizard install
        /// If this is not null/empty, we can usually make a safe assumption that the ProteoWizard dlls are available.
        /// </summary>
        public static string PwizPath => ProteoWizardReaderImplementation.PwizPath;

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

        /// <summary>
        /// Gets all spectra
        /// </summary>
        /// <returns>all spectra</returns>
        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            return _pwizReader.ReadAllSpectra();
        }

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000767, native spectrum identifier format
        /// </summary>
        public PSI_Interface.CV.CV.CVID NativeIdFormat => _pwizReader.NativeIdFormat;

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </summary>
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
        /// Returns the spectrum specified by the scan number.
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="includePeaks"></param>
        /// <returns></returns>
        public Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true)
        {
            return _pwizReader.ReadMassSpectrum(scanNum, includePeaks);
        }

        /// <summary>
        /// Close the reader
        /// </summary>
        public void Close()
        {
            _pwizReader.Dispose();
        }

        /// <summary>
        /// The number of spectra in the file.
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

        /// <summary>
        /// All files that ProteoWizard supports, either directly, or as part of a folder dataset
        /// </summary>
        public static List<string> SupportedFilesFilterList => new List<string>()
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

        /// <summary>
        /// List of "folder extensions" that ProteoWizard can read. This does not include all folder type datasets - some require directory listings.
        /// </summary>
        public static List<string> SupportedDirectoryTypes => new List<string>() { ".d", ".raw" };

        /// <summary>
        /// List of files that are produced by Bruker instruments that ProteoWizard can read.
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
    internal sealed class ProteoWizardReaderImplementation : IMassSpecDataReader
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
                AppDomain.CurrentDomain.AssemblyResolve += ProteoWizardAssemblyResolver;
                _resolverAdded = true;
            }
        }

        private static bool _resolverAdded;

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

            // Retrieve the list of referenced assemblies in an array of AssemblyName.
            var strTempAssmbPath = "";

            var arrReferencedAssmbNames = Assembly.GetExecutingAssembly().GetReferencedAssemblies();

            // Loop through the array of referenced assembly names.
            foreach (var strAssmbName in arrReferencedAssmbNames)
            {
                //Check for the assembly names that have raised the "AssemblyResolve" event.
                if (strAssmbName.FullName.Substring(0, strAssmbName.FullName.IndexOf(',')) == args.Name.Substring(0, args.Name.IndexOf(',')))
                {
                    //Console.WriteLine("Attempting to load DLL \"" + Path.Combine(pwizPath, args.Name.Substring(0, args.Name.IndexOf(",")) + ".dll") + "\"");
                    //Build the path of the assembly from where it has to be loaded.
                    strTempAssmbPath = Path.Combine(PwizPath, args.Name.Substring(0, args.Name.IndexOf(',')) + ".dll");
                    break;
                }
            }
#if DEBUG
            Console.WriteLine("Loading file \"" + strTempAssmbPath + "\"");
#endif
            var assemblyFile = new FileInfo(strTempAssmbPath);

            // Load the assembly from the specified path.
            Assembly myAssembly;
            try
            {
                myAssembly = Assembly.LoadFrom(assemblyFile.FullName);
            }
            catch (BadImageFormatException)
            {
                Console.WriteLine("Incompatible Assembly: \"" + assemblyFile.FullName + "\"");
                throw;
            }
            catch (FileNotFoundException)
            {
                Console.WriteLine("Assembly not found: \"" + assemblyFile.FullName + "\"");
                throw;
            }
            catch (FileLoadException)
            {
                Console.WriteLine("Invalid Assembly: \"" + assemblyFile.FullName + "\"");
                Console.WriteLine("The assembly may be marked as \"Untrusted\" by Windows. Please unblock and try again.");
                Console.WriteLine("Use the Streams tool (https://technet.microsoft.com/en-us/sysinternals/streams.aspx) to unblock, for example");
                if (assemblyFile.DirectoryName == null)
                    Console.WriteLine("streams -d *");
                else
                    Console.WriteLine("streams -d \"" + Path.Combine(assemblyFile.DirectoryName, "*") + "\"");
                throw;
            }
            catch (SecurityException)
            {
                Console.WriteLine("Assembly access denied: \"" + assemblyFile.FullName + "\"");
                throw;
            }

            //Return the loaded assembly.
            return myAssembly;
        }

        #endregion

        #region Static stateful variable and populating functions

        /// <summary>
        /// Name of the DLL we are checking for
        /// </summary>
        public const string TargetDllName = "pwiz_bindings_cli.dll";

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
        /// "%ProteoWizard%" or "%ProteoWizard%_x86" environment variable data,
        /// "C:\DMS_Programs\ProteoWizard" or "C:\DMS_Programs\ProteoWizard_x86",
        /// "%ProgramFiles%\ProteoWizard\(highest sorted)"</remarks>
        public static string FindPwizPath()
        {
            string pwizPath;

            // Set the DMS_Programs ProteoWizard path based on if the process is 32- or 64-bit.
            string dmsProgPwiz;

            if (!Environment.Is64BitProcess)
            {
                // Check for a x86 ProteoWizard environment variable
                pwizPath = Environment.GetEnvironmentVariable("ProteoWizard_x86");

                if (string.IsNullOrEmpty(pwizPath) && !Environment.Is64BitOperatingSystem)
                {
                    pwizPath = Environment.GetEnvironmentVariable("ProteoWizard");
                }

                dmsProgPwiz = @"C:\DMS_Programs\ProteoWizard_x86";
            }
            else
            {
                // Check for a x64 ProteoWizard environment variable
                pwizPath = Environment.GetEnvironmentVariable("ProteoWizard");
                dmsProgPwiz = @"C:\DMS_Programs\ProteoWizard";
            }

            if (string.IsNullOrWhiteSpace(pwizPath) && Directory.Exists(dmsProgPwiz) &&
                new DirectoryInfo(dmsProgPwiz).GetFiles(TargetDllName).Length > 0)
            {
                return dmsProgPwiz;
            }

            if (!string.IsNullOrWhiteSpace(pwizPath) && Directory.Exists(pwizPath) && new DirectoryInfo(pwizPath).GetFiles(TargetDllName).Length > 0)
            {
                return pwizPath;
            }

            // NOTE: This call returns the 32-bit Program Files folder if the running process is 32-bit
            // or the 64-bit Program Files folder if the running process is 64-bit
            var progFiles = Environment.GetEnvironmentVariable("ProgramFiles");
            if (string.IsNullOrWhiteSpace(progFiles))
            {
                return null;
            }

            // Construct a path of the form "C:\Program Files\ProteoWizard" or "C:\Program Files (x86)\ProteoWizard"
            var progPwiz = Path.Combine(progFiles, "ProteoWizard");
            var pwizFolder = new DirectoryInfo(progPwiz);
            if (pwizFolder.Exists)
            {
                if (pwizFolder.GetFiles(TargetDllName).Length > 0)
                {
                    return progPwiz;
                }
            }
            else
            {
                // Update pwizFolder to be "C:\Program Files" or "C:\Program Files (x86)"
                pwizFolder = new DirectoryInfo(progFiles);
                if (!pwizFolder.Exists)
                {
                    return null;
                }
            }

            // Look for subfolders whose names start with ProteoWizard, for example "ProteoWizard 3.0.9490"
            var subFolders = pwizFolder.GetDirectories("ProteoWizard*").ToList();

            if (subFolders.Count <= 0)
            {
                return null;
            }

            // Try to sort by version, it properly handles the version rolling over powers of 10 (but string sorting does not)
            var byVersion = new List<Tuple<System.Version, DirectoryInfo>>();
            foreach (var folder in subFolders)
            {
                try
                {
                    // Just ignoring the directory here if it has no version
                    var version = System.Version.Parse(folder.Name.Trim().Split(' ').Last());
                    byVersion.Add(new Tuple<System.Version, DirectoryInfo>(version, folder));
                }
                catch (Exception)
                {
                    // Do nothing...
                }
            }
            if (byVersion.Count > 0)
            {
                byVersion.Sort((x, y) => x.Item1.CompareTo(y.Item1));
                byVersion.Reverse();
                var subFoldersOrig = subFolders.ToArray();
                subFolders = byVersion.Select(x => x.Item2).ToList();
                // Guarantee that any folder where we couldn't parse a version is in the list, but at the end.
                foreach (var folder in subFoldersOrig)
                {
                    if (!subFolders.Contains(folder))
                    {
                        subFolders.Add(folder);
                    }
                }
            }
            else
            {
                // Sorting by version failed, try the old method.
                subFolders.Sort((x, y) => string.Compare(x.FullName, y.FullName, StringComparison.Ordinal));
                subFolders.Reverse(); // reverse the sort order - this should give us the highest installed version of ProteoWizard first
            }

            foreach (var folder in subFolders)
            {
                if (folder.GetFiles(TargetDllName).Length > 0)
                {
                    return folder.FullName;
                }
            }
            // If the above failed, return the highest version installed
            return subFolders[0].FullName;
        }

        static ProteoWizardReaderImplementation()
        {
            PwizPath = FindPwizPath();
        }

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

                Console.WriteLine(message);
                throw new TypeLoadException(message);
            }
        }

        private static void ValidateLoaderByPath()
        {
            if (string.IsNullOrWhiteSpace(PwizPath))
            {
                var message = CannotFindExceptionMessage();

                Console.WriteLine(message);
                throw new TypeLoadException(message);
            }
        }

        private static string CannotFindExceptionMessage()
        {
            var bits = Environment.Is64BitProcess ? "64" : "32";
            var message = "Cannot load ProteoWizard dlls. Please ensure that " + bits
                + "-bit ProteoWizard is installed to its default install directory (\""
                + Environment.GetEnvironmentVariable("ProgramFiles") + "\\ProteoWizard\\ProteoWizard 3.0.[x]\")."
                + "\nCurrently trying to load ProteoWizard dlls from path \"" + PwizPath + "\".";

            return message;
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
        private void LoadPwizReader()
        {
            if (_loaded)
            {
                return;
            }

            var readers = ReaderList.FullReaderList;
            readers.read(FilePath, _dataFile);
            if (new[] { ".mzml", ".mzml.gz", ".mzxml", ".mzxml.gz", ".mgf", ".mgf.gz", ".txt", "uimf", "uimf.gz" }
                .Any(ext => FilePath.ToLower().EndsWith(ext)))
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
                    _srcFileChecksum = BitConverter.ToString(hash).ToLower().Replace("-", "");
                }
            }

            foreach (var software in _dataFile.softwareList)
            {
                _fileFormatVersion = software.version;
            }

            _loaded = true;
        }

        #endregion

        #region IMassSpecDataReader implementation

        /// <summary>
        /// Gets all spectra
        /// </summary>
        /// <returns>all spectra</returns>
        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            LoadPwizReader();
            for (var i = 1; i <= _numSpectra; i++)
            {
                yield return ReadSpectrum(i);
            }
        }

        /// <summary>
        /// The number of spectra in the file.
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
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000767, native spectrum identifier format
        /// </summary>
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
        /// The Native Format of the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </summary>
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
        /// Returns the spectrum specified by the scan number.
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="includePeaks"></param>
        /// <returns></returns>
        public Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true)
        {
            LoadPwizReader();
            return ReadSpectrum(scanNum, includePeaks);
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
            var mzArray = new double[0];
            var intensityArray = new double[0];
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
