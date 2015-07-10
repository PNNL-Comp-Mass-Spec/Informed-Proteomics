using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace InformedProteomics.Backend.MassSpecData
{
    public class MassSpecDataReaderFactory
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
                    reader = new XCaliburReader(filePath);
                    break;
                case MassSpecDataType.MzMLFile:
                    reader = new MzMLReader(filePath);
                    break;
                case MassSpecDataType.PbfFile:
                    reader = new PbfLcMsRun(filePath);
                    break;
                case MassSpecDataType.Unknown:
                    AppDomain.CurrentDomain.AssemblyResolve += ProteoWizardReader.ProteoWizardAssemblyResolver;
                    reader = new ProteoWizardReader(filePath);
                    break;
            }

            return reader;
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
            var lower = filePath.ToLower();
            if (lower.EndsWith("\\") || lower.EndsWith("/"))
            {
                lower = lower.Remove(lower.Length - 1);
            }

            if (lower.EndsWith(".raw"))
            {
                return MassSpecDataType.XCaliburRun;
            }
            if (lower.EndsWith(".mzml") || lower.EndsWith(".mzml.gz"))
            {
                return MassSpecDataType.MzMLFile;
            }
            if (lower.EndsWith(".pbf"))
            {
                return MassSpecDataType.PbfFile;
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
                return "All Supported|*.raw;*.mzML;*.mzML.gz;*.mzXML;*.mzXML.gz;*.mgf;*.mgf.gz;*.d;mspeak.bin;msprofile.bin;*.wiff;*.d;*.u2;FID;analysis.yep;analysis.baf;*.raw;_extern.inf;_inlet.inf;_FUNC001.DAT;*.lcd;*.pbf"
                       + "|Thermo .RAW|*.raw"
                       + "|mzML[.gz]|*.mzML;*.mzML.gz"
                       + "|mzXML[.gz]|*.mzXML;*.mzXML.gz"
                       + "|MGF[.gz]|*.mgf;*.mgf.gz"
                       + "|Agilent .d|*.d;mspeak.bin;msprofile.bin"
                       + "|AB Sciex .wiff|*.wiff"
                       + "|Bruker .d/FID/YEP/BAF|*.d;*.u2;FID;analysis.yep;analysis.baf"
                       + "|Waters .raw|*.raw;_extern.inf;_inlet.inf;_FUNC001.DAT"
                       + "|Shimadzu lcd|*.lcd"
                       + "|PBF|*.pbf"
                    ;
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
                var internalSupported = new List<string>() {".raw", ".mzml", ".mzml.gz", ".pbf"};
                var pwizSupported = ProteoWizardReader.SupportedFilesFilterList;
                return internalSupported.Concat(pwizSupported.Where(ext => !internalSupported.Contains(ext))).ToList();
            }
        }

        /// <summary>
        /// The list of directory dataset type extensions that are supported by ProteoWizard (not all-inclusive)
        /// </summary>
        public static List<string> SupportedDirectoryTypes
        {
            get { return ProteoWizardReader.SupportedDirectoryTypes; }
        }

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
                    int pos = ext.IndexOf('.', 1);
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
    }
}
