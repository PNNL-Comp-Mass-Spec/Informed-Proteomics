using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using PRISM;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Factory class for getting the right data reader for the provided file(s)
    /// </summary>
    public static class MassSpecDataReaderFactory
    {
        /// <summary>
        /// Gets the appropriate IMassSpecDataReader for the supplied path.
        /// It is recommended that "NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        /// <remarks>It is recommended that "NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path</remarks>
        public static IMassSpecDataReader GetMassSpecDataReader(string filePath)
        {
            filePath = NormalizeDatasetPath(filePath);
            var type = GetMassSpecDataType(filePath);
            IMassSpecDataReader reader = null;
            switch (type)
            {
                case MassSpecDataType.XCaliburRun:
                    reader = new XcaliburReader(filePath);
                    break;
                case MassSpecDataType.MzMLFile:
                    reader = new MzMLReader(filePath);
                    break;
                case MassSpecDataType.PbfFile:
                    reader = new PbfLcMsRun(filePath);
                    break;
                case MassSpecDataType.DeconvolutedPbfFile:
                    reader = new DPbfLcMsRun(filePath);
                    break;
                case MassSpecDataType.Unknown:
                    if (_pwizAvailable)
                    {
                        reader = new ProteoWizardReader(filePath);
                    }
                    else
                    {
                        var arch = Environment.Is64BitProcess ? "64" : "32";
                        ConsoleMsgUtils.ShowWarning(string.Format("WARNING: Could not find a reader for file \"{0}\"." +
                                                                  " Is ProteoWizard {1}-bit installed?", filePath, arch));
                    }
                    break;
            }

            return reader;
        }

        public static ISpectrumAccessor GetMassSpecDataAccessor(string filePath, IProgress<ProgressData> progress = null)
        {
            var reader = GetMassSpecDataReader(filePath);
            if (reader == null)
            {
                return null;
            }

            if (!reader.TryMakeRandomAccessCapable())
            {
                // Convert to PBF to support random access
                reader = PbfLcMsRun.GetLcMsRun(filePath, progress);
            }

            if (reader is ISpectrumAccessor specAcc)
            {
                return specAcc;
            }

            return new SpectrumAccessorWrapper(reader, progress);
        }

        /// <summary>
        /// Gets the appropriate classifying MassSpecDataType for the supplied path
        /// </summary>
        /// <param name="filePath">Path to spec file/folder</param>
        /// <returns></returns>
        public static MassSpecDataType GetMassSpecDataType(string filePath)
        {
            // Calls "NormalizeDatasetPath" to make sure we save the file to the containing directory
            filePath = NormalizeDatasetPath(filePath);
            var pathLower = filePath.ToLower();
            if (pathLower.EndsWith("\\") || pathLower.EndsWith("/"))
            {
                pathLower = pathLower.Remove(pathLower.Length - 1);
            }

            if (_thermoRawAvailable && pathLower.EndsWith(".raw") && !Directory.Exists(filePath))
            {
                return MassSpecDataType.XCaliburRun;
            }
            // If it is an mzML file that is gzipped, use the ProteoWizardReader, since it does random access without extracting.
            if (pathLower.EndsWith(".mzml") || (pathLower.EndsWith(".mzml.gz") && !_pwizAvailable))
            {
                return MassSpecDataType.MzMLFile;
            }
            if (pathLower.EndsWith(".pbf"))
            {
                return MassSpecDataType.PbfFile;
            }
            if (pathLower.EndsWith(".dpbf"))
            {
                return MassSpecDataType.DeconvolutedPbfFile;
            }

            return MassSpecDataType.Unknown;
        }

        /// <summary>
        /// The list of all formats supported by built-in type and ProteoWizard (not all-inclusive), as an OpenFileDialog filter string.
        /// </summary>
        /// <remarks>Included are some filenames, which are inner contents for folder-type datasets.</remarks>
        public static string MassSpecDataTypeFilterString
        {
            get
            {
                if (_pwizAvailable)
                {
                    return FilterStringAll;
                }
                if (_thermoRawAvailable)
                {
                    return FilterStringBuiltIn;
                }
                return FilterStringNoExternalDll;
            }
        }

        /// <summary>
        /// The list of all formats supported by built-in type and ProteoWizard (not all-inclusive).
        /// </summary>
        /// <remarks>Included are some filenames, which are inner contents for folder-type datasets.</remarks>
        public static List<Tuple<string, string[]>> MassSpecDataTypes
        {
            get
            {
                if (_pwizAvailable)
                {
                    return SupportedTypesAll;
                }
                if (_thermoRawAvailable)
                {
                    return SupportedTypesBuiltIn;
                }
                return SupportedTypesNoExternalDll;
            }
        }

        /// <summary>
        /// The list of extensions (some file names) of all formats supported by built-in types and by ProteoWizard (not all-inclusive)
        /// </summary>
        /// <remarks>Included are some filenames, which are inner contents for folder-type datasets.</remarks>
        public static List<string> MassSpecDataTypeFilterList
        {
            get
            {
                // These file extensions need to be lowercase due to Linq queries of the form "types.Select(ext => SpecFilePath.ToLower().EndsWith(ext)).Any()"
                var internalSupported = new List<string> {".mzml", ".mzml.gz", ".pbf"};
                if (_thermoRawAvailable)
                {
                    internalSupported.Add(".raw");
                }
                if (_pwizAvailable)
                {
                    var pwizSupported = ProteoWizardReader.SupportedFilesFilterList;
                    internalSupported = internalSupported.Concat(pwizSupported.Where(ext => !internalSupported.Contains(ext))).ToList();
                }
                return internalSupported;
            }
        }

        /// <summary>
        /// The list of directory dataset type extensions that are supported by ProteoWizard (not all-inclusive)
        /// </summary>
        public static List<string> SupportedDirectoryTypes => ProteoWizardReader.SupportedDirectoryTypes;

        /// <summary>
        /// Test the supplied path to see if we can read it using available readers
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public static bool IsADirectoryDataset(string path)
        {
            path = NormalizeDatasetPath(path);
            return SupportedDirectoryTypes.Any(f => path.ToLower().EndsWith(f));
        }

        /// <summary>
        /// Gets the directory that contains the dataset; it will back out of subdirectories of folder-type datasets
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public static string GetDirectoryContainingDataset(string path)
        {
            return Path.GetDirectoryName(GetDatasetName(path));
        }

        /// <summary>
        /// Utility function: Gets the actual dataset name, needed in cases when a file in a folder-type dataset is given.
        /// </summary>
        /// <param name="path"></param>
        /// <returns>The path to the dataset, or the path to the dataset directory</returns>
        public static string GetDatasetName(string path)
        {
            return ProteoWizardReader.CheckForDirectoryDataset(path);
        }

        /// <summary>
        /// Get the normalized dataset path; only returns a different string when the given path is to a file/folder in a folder-type dataset
        /// </summary>
        /// <param name="path"></param>
        /// <returns>Normalized dataset path</returns>
        public static string NormalizeDatasetPath(string path)
        {
            return ProteoWizardReader.CheckForDirectoryDataset(path);
        }

        /// <summary>
        /// Modification of Path.ChangeExtension: properly removes multiple extensions.
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="newExt"></param>
        /// <returns></returns>
        public static string ChangeExtension(string filePath, string newExt)
        {
            var path = filePath;
            foreach (var ext in MassSpecDataTypeFilterList)
            {
                if (path.ToLower().EndsWith(ext))
                {
                    var pos = ext.IndexOf('.', 1);

                    // Remove extra extensions
                    if (pos > 0)
                    {
                        path = path.Remove(path.Length - (ext.Length - pos));
                    }
                    return Path.ChangeExtension(path, newExt);
                }
            }
            return Path.ChangeExtension(path, newExt);
        }

        /// <summary>
        /// Modification of Path.ChangeExtension(path, null), which removes the extension. This will remove multiple extensions.
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        public static string RemoveExtension(string filePath)
        {
            return ChangeExtension(filePath, null);
        }

        static MassSpecDataReaderFactory()
        {
            _pwizAvailable = IsPwizAvailable();
            _thermoRawAvailable = IsThermoRawAvailable();
        }

        private static readonly bool _pwizAvailable;
        private static readonly bool _thermoRawAvailable;

        /// <summary>
        /// Tests to see if we can load the needed ProteoWizard DLL without errors
        /// </summary>
        /// <returns>True if we can load the ProteoWizard DLLs</returns>
        public static bool IsPwizAvailable()
        {
            // TODO: Test this on a system without ProteoWizard installed...
            // Also, modify test to actually try loading the pwiz_bindings_cli
            var pwizPath = ProteoWizardReader.PwizPath;
            if (!string.IsNullOrWhiteSpace(pwizPath))
            {
                var pwizBindingsCLI = Path.Combine(pwizPath, "pwiz_bindings_cli.dll");
                try
                {
                    var pwizCLI = Assembly.LoadFrom(pwizBindingsCLI);
                    // ReSharper disable once ConditionIsAlwaysTrueOrFalse
                    if (pwizCLI != null)
                    {
                        return true;
                    }
                }
                catch (Exception)
                {
                    return false;
                }
            }
            return false;
        }

        /// <summary>
        /// Tests to see if we can can use RawFileReader (i.e. are we running as 64-bit process)
        /// </summary>
        /// <returns></returns>
        public static bool IsThermoRawAvailable()
        {
            if (Environment.Is64BitProcess)
            {
                return true;
            }

            return false;
        }

        /// <summary>
        /// Gets the list of all file types supported as a pair of description and file extensions.
        /// Contains all files supported when ProteoWizard is available.
        /// </summary>
        private static List<Tuple<string, string[]>> SupportedTypesAll => new List<Tuple<string, string[]>>
        {
            new Tuple<string, string[]>("Thermo .RAW", new[] { ".raw" }),
            new Tuple<string, string[]>("mzML", new[] { ".mzml", ".mzml.gz" }),
            new Tuple<string, string[]>("mzXML", new[] { ".mzXML", ".mzXML.gz" }),
            new Tuple<string, string[]>("Mascot Generic Format", new[] { ".mgf", ".mgf.gz" }),
            new Tuple<string, string[]>("Agilent", new[] { ".d", "mspeak.bin", "msprofile.bin" }),
            new Tuple<string, string[]>("AB Sciex", new[] { ".wiff" }),
            new Tuple<string, string[]>("Bruker", new[] { ".d", ".u2", "FID", "analysis.yep", "analysis.baf" }),
            new Tuple<string, string[]>("Waters", new[] { ".raw", "_extern.inf", "_inlet.inf", ".DAT" }),
            new Tuple<string, string[]>("Shimadzu", new[] { ".lcd" }),
            new Tuple<string, string[]>("PNNL Binary Format", new[] { ".pbf" }),
        };

        /// <summary>
        /// Gets the list of all file types supported as a pair of description and file extensions.
        /// Contains all file types supported if only ThermoRawFileReaderDLL is available.
        /// </summary>
        private static List<Tuple<string, string[]>> SupportedTypesBuiltIn => new List<Tuple<string, string[]>>
        {
            new Tuple<string, string[]>("Thermo .RAW", new[] { ".raw" }),
            new Tuple<string, string[]>("mzMl", new[] { ".mzml", ".mzml.gz" }),
            new Tuple<string, string[]>("PNNL Binary Format", new[] { ".pbf" }),
        };

        /// <summary>
        /// Gets the list of all file types supported as a pair of description and file extensions.
        /// Contains only file types that are natively supported by InformedProteomics without
        /// any external DLLs available.
        /// </summary>
        private static List<Tuple<string, string[]>> SupportedTypesNoExternalDll => new List<Tuple<string, string[]>>
        {
            new Tuple<string, string[]>("mzMl", new[] { ".mzml", ".mzml.gz" }),
            new Tuple<string, string[]>("PNNL Binary Format", new[] { ".pbf" }),
        };

        private static string FilterStringAll => "All Supported Spectrum Files|*.raw;*.mzML;*.mzML.gz;*.mzXML;*.mzXML.gz;*.mgf;*.mgf.gz;*.d;mspeak.bin;msprofile.bin;*.wiff;*.d;*.u2;FID;analysis.yep;analysis.baf;*.raw;_extern.inf;_inlet.inf;_FUNC*.DAT;*.lcd;*.uimf;*.pbf"
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
                                                 + "|PBF|*.pbf";

        private static string FilterStringBuiltIn => "All Supported Spectrum Files|*.raw;*.mzML;*.mzML.gz;*.pbf"
                                                     + "|Thermo .RAW|*.raw"
                                                     + "|mzML[.gz]|*.mzML;*.mzML.gz"
                                                     + "|PBF|*.pbf";

        private static string FilterStringNoExternalDll => "All Supported Spectrum Files|*.mzML;*.mzML.gz;*.pbf"
                                                           + "|mzML[.gz]|*.mzML;*.mzML.gz"
                                                           + "|PBF|*.pbf";
    }
}
