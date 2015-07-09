using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace InformedProteomics.Backend.MassSpecData
{
    public class MassSpecDataReaderFactory
    {
        public static IMassSpecDataReader GetMassSpecDataReader(string filePath)
        {
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

        public static MassSpecDataType GetMassSpecDataType(string filePath)
        {
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

        public static List<string> MassSpecDataTypeFilterList
        {
            get
            {
                var internalSupported = new List<string>() {".raw", ".mzml", ".mzml.gz", ".pbf"};
                var pwizSupported = ProteoWizardReader.SupportedFilesFilterList;
                return internalSupported.Concat(pwizSupported.Where(ext => !internalSupported.Contains(ext))).ToList();
            }
        }

        public static List<string> SupportedDirectoryTypes
        {
            get { return ProteoWizardReader.SupportedDirectoryTypes; }
        }

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

        public static string RemoveExtension(string filePath)
        {
            return ChangeExtension(filePath, null);
        }
    }
}
