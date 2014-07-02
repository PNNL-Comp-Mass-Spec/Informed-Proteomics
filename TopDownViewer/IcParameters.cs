using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDownViewer
{
    public class IcParameters
    {
        public LcMsRun Lcms { get; set; }
        public string DatabaseFile { get; set; }
        public int SearchMode { get; set; }
        public bool Tda { get; set; }
        public Tolerance PrecursorTolerancePpm { get; set; }
        public Tolerance ProductIonTolerancePpm { get; set; }
        public int MinSequenceLength { get; set; }
        public int MaxSequenceLength { get; set; }
        public int MinPrecursorIonCharge { get; set; }
        public int MaxPrecursorIonCharge { get; set; }
        public int MinProductIonCharge { get; set; }
        public int MaxProductIonCharge { get; set; }
        public int MinSequenceMass { get; set; }
        public int MaxSequenceMass { get; set; }
        public int MaxDynamicModificationsPerSequence { get; set; }
        public List<SearchModification> Modifications { get; set; }
        public IonTypeFactory IonTypeFactory { get; set; }

        public IcParameters(string paramFile)
        {
            Modifications = new List<SearchModification>();
            var file = File.ReadLines(paramFile);

            foreach (var line in file)
            {
                var parts = line.Split('\t');
                if (parts.Length < 2) throw new Exception("Invalid configuration file.");
                switch (parts[0])
                {
                    case "SpecFile":
                        Lcms = LcMsRun.GetLcMsRun(parts[1], MassSpecDataType.XCaliburRun, 0, 0);
                        break;
                    case "DatabaseFile":
                        DatabaseFile = parts[1];
                        break;
                    case "SearchMode":
                        SearchMode = Convert.ToInt32(parts[1]);
                        break;
                    case "Tda":
                        Tda = Convert.ToBoolean(parts[1]);
                        break;
                    case "PrecursorIonTolerancePpm":
                        PrecursorTolerancePpm = new Tolerance(Convert.ToDouble(parts[1]), ToleranceUnit.Ppm);
                        break;
                    case "ProductIonTolerancePpm":
                        ProductIonTolerancePpm = new Tolerance(Convert.ToDouble(parts[1]), ToleranceUnit.Ppm);
                        break;
                    case "MinSequenceLength":
                        MinSequenceLength = Convert.ToInt32(parts[1]);
                        break;
                    case "MaxSequenceLength":
                        MaxSequenceLength = Convert.ToInt32(parts[1]);
                        break;
                    case "MinPrecursorIonCharge":
                        MinPrecursorIonCharge = Convert.ToInt32(parts[1]);
                        break;
                    case "MaxPrecursorIonCharge":
                        MaxPrecursorIonCharge = Convert.ToInt32(parts[1]);
                        break;
                    case "MinProductIonCharge":
                        MinProductIonCharge = Convert.ToInt32(parts[1]);
                        break;
                    case "MaxProductIonCharge":
                        MaxProductIonCharge = Convert.ToInt32(parts[1]);
                        IonTypeFactory = new IonTypeFactory(MaxProductIonCharge);
                        break;
                    case "MinSequenceMass":
                        MinSequenceMass = Convert.ToInt32(parts[1]);
                        break;
                    case "MaxSequenceMass":
                        MaxSequenceMass = Convert.ToInt32(parts[1]);
                        break;
                    case "MaxDynamicModificationsPerSequence":
                        MaxDynamicModificationsPerSequence = Convert.ToInt32(parts[1]);
                        break;
                    case "Modification":
                        var modParts = parts[1].Split(',');
//                        var compositionStr = modParts[0];
                        var residue = modParts[1];
                        var isFixed = modParts[2];
                        var position = modParts[3];
                        var modName = modParts[4];
                        var modification = Modification.Get(modName);
                        SequenceLocation sequenceLocation;
                        Enum.TryParse(position, out sequenceLocation);
//                        var composition = Composition.Parse(compositionStr);
                        var isFixedModification = (isFixed == "fix");
                        Modifications.Add(new SearchModification(modification, residue.ToCharArray()[0],
                                                                 sequenceLocation, isFixedModification));
                        break;
                }
            }
        }
    }
}
