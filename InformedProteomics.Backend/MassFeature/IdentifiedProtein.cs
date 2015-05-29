using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.MassFeature
{
    public class IdentifiedProtein
    {
        public IdentifiedProtein(double mass, int charge, int scan, string sequenceStr, string modStr, int datasetId = 0, int id = 0)
        {
            ScanNum = scan;
            Mass = mass;
            Charge = charge;
            
            MinCharge = charge;
            MaxCharge = charge;
            MinScanNum = scan;
            MaxScanNum = scan;
            DataSetId = datasetId;
            Id = id;
            Sequence = sequenceStr;
            Modifications = modStr;

            var seqKeyPair = SetModifications(sequenceStr, modStr);
            SequenceText = seqKeyPair.Item2;
            _sequence = seqKeyPair.Item1;
        }

        
        public int ScanNum { get; set; }
        public int Charge { get; set; }
        public double Mass { get; set; }
        public int DataSetId { get; set; }
        public int Id { get; set; }

        public string Sequence { get; set; }
        public string Modifications { get; set; }
        
        public string SequenceText { get; set; }

        public string ProteinName { get; set; }
        public string Composition { get; set; }

        public int MinCharge { get; set; }
        public int MaxCharge { get; set; }
        public int ChargeLength { get { return (MaxCharge == 0) ? 0 : MaxCharge - MinCharge + 1; } }

        public int MinScanNum { get; set; }
        public int MaxScanNum { get; set; }
        public int ScanLength { get { return (MaxScanNum == 0) ? 0 : MaxScanNum - MinScanNum + 1; } }

        public List<string> Rows = new List<string>();
        private readonly Sequence _sequence;
        public Sequence GetSequence()
        {
            return _sequence; 
        }

        private Tuple<Sequence, string> SetModifications(string cleanSequence, string modifications)
        {
            // Build Sequence AminoAcid list
            var sequence = new Sequence(cleanSequence, new AminoAcidSet());
            var sequenceText = cleanSequence;
            var parsedModifications = ParseModifications(modifications);

            // Add modifications to sequence
            parsedModifications.Sort(new CompareModByHighestPosition());   // sort in reverse order for insertion
            foreach (var mod in parsedModifications)
            {
                var pos = mod.Item1;
                if (pos > 0)
                {
                    pos--;
                }

                var modLabel = string.Format("[{0}]", mod.Item2.Name);
                sequenceText = sequenceText.Insert(mod.Item1, modLabel);
                var aa = sequence[pos];
                var modaa = new ModifiedAminoAcid(aa, mod.Item2);
                sequence[pos] = modaa;
            }

            return new Tuple<Sequence, string>(new Sequence(sequence), sequenceText);
        }

        private List<Tuple<int, Modification>> ParseModifications(string modifications)
        {
            var mods = modifications.Split(',');
            var parsedMods = new List<Tuple<int, Modification>>();
            if (mods.Length < 1 || mods[0] == string.Empty)
            {
                return parsedMods;
            }

            foreach (var modParts in mods.Select(mod => mod.Split(' ')))
            {
                if (modParts.Length < 0)
                {
                    throw new FormatException("Unknown Modification");
                }

                var modName = modParts[0];
                var modPos = Convert.ToInt32(modParts[1]);
                var modification = Modification.Get(modName);
                parsedMods.Add(new Tuple<int, Modification>(modPos, modification));
                if (modification == null)
                {
                    throw new Exception(string.Format("Found an unrecognized modification: {0}", modName));
                }
            }

            return parsedMods;
        }

        internal class CompareModByHighestPosition : IComparer<Tuple<int, Modification>>
        {
            public int Compare(Tuple<int, Modification> x, Tuple<int, Modification> y)
            {
                return y.Item1.CompareTo(x.Item1);
            }
        }    


    }
}
