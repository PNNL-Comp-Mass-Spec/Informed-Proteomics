using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Security;
using PRISM;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// The static methods in this class are used by ProteoWizardReader to find the ProteoWizard DLLs
    /// </summary>
    /// <remarks>
    /// These methods were originally in the ProteoWizardReader class,
    /// but were moved here to prevent mono (on Linux) from trying to resolve pwiz_bindings_cli.dll
    /// when other methods in class ProteoWizardReader are called
    /// </remarks>
    internal static class ProteoWizardLoader
    {
        // Ignore Spelling: pwiz_bindings_cli

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
        /// On a missing DLL event, searches a path specified by FindPwizPath for the ProteoWizard DLLs, and loads them
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
            if (!args.Name.StartsWith("pwiz_bindings_cli", StringComparison.OrdinalIgnoreCase))
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
            var tempAssemblyPath = string.Empty;

            var referencedAssemblyNames = Assembly.GetExecutingAssembly().GetReferencedAssemblies();

            // Loop through the array of referenced assembly names.
            foreach (var assemblyName in referencedAssemblyNames)
            {
                //Check for the assembly names that have raised the "AssemblyResolve" event.
                if (assemblyName.FullName.Substring(0, assemblyName.FullName.IndexOf(',')) == args.Name.Substring(0, args.Name.IndexOf(',')))
                {
                    //Console.WriteLine("Attempting to load DLL \"" + Path.Combine(pwizPath, args.Name.Substring(0, args.Name.IndexOf(",")) + ".dll") + "\"");
                    //Build the path of the assembly from where it has to be loaded.
                    tempAssemblyPath = Path.Combine(PwizPath, args.Name.Substring(0, args.Name.IndexOf(',')) + ".dll");
                    break;
                }
            }
#if DEBUG
            Console.WriteLine("Loading file \"" + tempAssemblyPath + "\"");
#endif
            var assemblyFile = new FileInfo(tempAssemblyPath);

            // Load the assembly from the specified path.
            Assembly myAssembly;
            try
            {
                myAssembly = Assembly.LoadFrom(assemblyFile.FullName);
            }
            catch (BadImageFormatException)
            {
                ConsoleMsgUtils.ShowError("Incompatible Assembly: \"" + assemblyFile.FullName + "\"");
                throw;
            }
            catch (FileNotFoundException)
            {
                ConsoleMsgUtils.ShowError("Assembly not found: \"" + assemblyFile.FullName + "\"");
                throw;
            }
            catch (FileLoadException)
            {
                string streamsCmdLine;
                if (assemblyFile.DirectoryName == null)
                {
                    streamsCmdLine = "streams -d *";
                }
                else
                {
                    streamsCmdLine = "streams -d \"" + Path.Combine(assemblyFile.DirectoryName, "*") + "\"";
                }

                var msg =
                    "Invalid Assembly: \"" + assemblyFile.FullName + "\"\n" +
                    "The assembly may be marked as \"Untrusted\" by Windows. Please unblock and try again.\n" +
                    "Use the Streams tool (https://technet.microsoft.com/en-us/sysinternals/streams.aspx) to unblock, for example\n" +
                    streamsCmdLine;

                ConsoleMsgUtils.ShowError(msg);
                throw;
            }
            catch (SecurityException)
            {
                ConsoleMsgUtils.ShowError("Assembly access denied: \"" + assemblyFile.FullName + "\"");
                throw;
            }

            //Return the loaded assembly.
            return myAssembly;
        }

        /// <summary>
        /// Name of the DLL we are checking for
        /// </summary>
        public const string TargetDllName = "pwiz_bindings_cli.dll";

        /// <summary>
        /// The path to the most recent 64-bit ProteoWizard install
        /// If this is not null/empty, we can usually make a safe assumption that the ProteoWizard DLLs are available.
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

            // Look for Per-user and per-machine ProteoWizard installs; use whichever install is newer.

            var possibleInstallDirs = new List<DirectoryInfo>();
            // Per-User ProteoWizard install detection
            var localAppDataDir = new DirectoryInfo(Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "Apps"));
            if (localAppDataDir.Exists)
            {
                var bitness = Environment.Is64BitProcess ? 64 : 32;
                possibleInstallDirs.AddRange(localAppDataDir.EnumerateDirectories($"ProteoWizard*{bitness}-bit"));
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

            // Look for subdirectories whose names start with ProteoWizard, for example "ProteoWizard 3.0.9490"
            possibleInstallDirs.AddRange(pwizFolder.EnumerateDirectories("ProteoWizard*"));

            if (possibleInstallDirs.Count == 0)
            {
                return null;
            }

            // Try to sort by version, it properly handles the version rolling over powers of 10 (but string sorting does not)
            var byVersion = new List<Tuple<System.Version, DirectoryInfo>>();
            foreach (var folder in possibleInstallDirs)
            {
                try
                {
                    // Just ignoring the directory here if it has no version
                    var versionString = folder.Name.Trim().Split(' ').Last();
                    if (folder.Name.EndsWith("-bit", StringComparison.OrdinalIgnoreCase))
                    {
                        var split = folder.Name.Trim().Split(' ');
                        versionString = split[split.Length - 2];
                    }

                    if (string.IsNullOrWhiteSpace(versionString) || !versionString.Contains("."))
                    {
                        continue;
                    }

                    var versionSplit = versionString.Split('.');

                    if (System.Version.TryParse(versionString, out var version))
                    {
                        // Old pre-Git SCM conversion install - only has 3 components
                        byVersion.Add(new Tuple<System.Version, DirectoryInfo>(version, folder));
                    }
                    else if (versionSplit.Length > 3 && System.Version.TryParse(string.Join(".", versionSplit.Take(versionSplit.Length - 1)), out var version2))
                    {
                        // Post-Git SCM conversion install - last section of the version is a Git hash, and will not parse
                        byVersion.Add(new Tuple<System.Version, DirectoryInfo>(version2, folder));
                    }
                }
                catch (Exception)
                {
                    // Do nothing...
                }
            }

            if (byVersion.Count > 0)
            {
                // Reverse sort the list
                byVersion.Sort((x, y) => y.Item1.CompareTo(x.Item1));
                var subFoldersOrig = possibleInstallDirs.ToArray();
                possibleInstallDirs = byVersion.ConvertAll(x => x.Item2);

                // Guarantee that any folder where we couldn't parse a version is in the list, but at the end.
                foreach (var folder in subFoldersOrig)
                {
                    if (!possibleInstallDirs.Contains(folder))
                    {
                        possibleInstallDirs.Add(folder);
                    }
                }
            }
            else
            {
                // Sorting by version failed, try the old method.
                // reverse the sort order - this should give us the highest installed version of ProteoWizard first
                possibleInstallDirs.Sort((x, y) => string.CompareOrdinal(y.FullName, x.FullName));
            }

            foreach (var folder in possibleInstallDirs)
            {
                if (folder.GetFiles(TargetDllName).Length > 0 && File.Exists(Path.Combine(folder.FullName, "pwiz_bindings_cli.dll")))
                {
                    return folder.FullName;
                }
            }
            // If the above failed, return the highest version installed
            return possibleInstallDirs[0].FullName;
        }

        static ProteoWizardLoader()
        {
            PwizPath = FindPwizPath();
        }

        /// <summary>
        /// Checks to make sure the path to ProteoWizard files is set. If not, throws an exception.
        /// </summary>
        /// <remarks>This function should generally only be called inside of a conditional statement to prevent the
        /// exception from being thrown when the ProteoWizard DLLs will not be needed.</remarks>
        public static void ValidateLoader()
        {
            try
            {
                Assembly.Load("pwiz_bindings_cli, Version=0.0.0.0, Culture=neutral, PublicKeyToken=null");
            }
            catch
            {
                var message = CannotFindExceptionMessage();

                ConsoleMsgUtils.ShowError(message);
                throw new TypeLoadException(message);
            }
        }

        private static void ValidateLoaderByPath()
        {
            if (string.IsNullOrWhiteSpace(PwizPath))
            {
                var message = CannotFindExceptionMessage();

                ConsoleMsgUtils.ShowError(message);
                throw new TypeLoadException(message);
            }
        }

        private static string CannotFindExceptionMessage()
        {
            var bits = Environment.Is64BitProcess ? "64" : "32";
            var message =
                "Cannot load ProteoWizard DLLs. Please ensure that " + bits +
                "-bit ProteoWizard is installed to its default install directory (\"" +
                Environment.GetEnvironmentVariable("ProgramFiles") + "\\ProteoWizard\\ProteoWizard 3.0.[x]\").\n" +
                "Currently trying to load ProteoWizard DLLs from path \"" + PwizPath + "\".";

            return message;
        }
    }
}
