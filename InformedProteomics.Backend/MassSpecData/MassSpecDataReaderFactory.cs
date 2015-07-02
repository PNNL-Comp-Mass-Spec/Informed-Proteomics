
using System.Collections.Generic;

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
                    //reader = new ProteoWizardReader(filePath);
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

        public static string GetMassSpecDataTypeFilterString()
        {
            return "All Supported|*.raw;*.mzML;*.mzML.gz;*.pbf"
                   + "|Thermo .RAW|*.raw"
                   + "|mzML[.gz]|*.mzML;*.mzML.gz"
                   + "|PBF|*.pbf";
        }

        public static List<string> GetMassSpecDataTypeFilterList()
        {
            return new List<string>() { ".raw", ".mzml", ".mzml.gz", ".pbf" };
        }
    }
}
