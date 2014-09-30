using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.SequenceTag
{
    public class IdentifiedSequenceTag
    {
        public int CleavageIndex1;
        public int CleavageIndex2;
        private readonly Sequence _sequence;
        public IonType[] IonTypeArray;
        public bool DerivedFromPrefix { get; private set; }

        public IdentifiedSequenceTag(Sequence seq, int c1, int c2, IonType[] iontypes)
        {
            CleavageIndex1 = c1;
            CleavageIndex2 = c2;
            DerivedFromPrefix = iontypes[0].IsPrefixIon;
            IonTypeArray = iontypes;
            _sequence = seq;
        }

        public int GetBeginResidue()
        {
            if (DerivedFromPrefix) return CleavageIndex1 + 1;
            else return _sequence.Count - CleavageIndex2 - 1;
        }

        public int GetEndResidue()
        {
            if (DerivedFromPrefix) return CleavageIndex2;
            else return _sequence.Count - CleavageIndex1 - 2;
        }

        public int GetLength()
        {
            return CleavageIndex2 - CleavageIndex1;
        }

        public string GetAnnotation(int cleavageIndex)
        {
            if (cleavageIndex < CleavageIndex1 || cleavageIndex > CleavageIndex2)
                throw new ArgumentOutOfRangeException();

            var t = IonTypeArray[cleavageIndex];

            return string.Format("{0}{1}{2}({3}+)", t.BaseIonType.Symbol, cleavageIndex, t.NeutralLoss.Name, t.Charge);
        }

        public string GetSequenceString()
        {
            var sb = new StringBuilder(GetLength());
            for (var i = GetBeginResidue(); i <= GetEndResidue(); i++) sb.Append(_sequence[i].Residue);

            return sb.ToString();
        }


        private static readonly AminoAcidSet aaSet = new AminoAcidSet();
        public static Sequence GenerateSequence(String seqStr, String modStr)
        {
            if (modStr == null || modStr.Equals(""))
                return new Sequence(seqStr, aaSet);

            int residueIndex = 0;
            string[] ptms = modStr.Split(',');
            var modificaitonInfo = new Dictionary<int, Modification>();

            foreach (var ptm in ptms)
            {
                string[] ptmPair = ptm.Split(' ');
                residueIndex = Int32.Parse(ptmPair[1]);
                if (residueIndex == 0) residueIndex++;

                modificaitonInfo.Add(residueIndex, Modification.Get(ptmPair[0]));
            }

            var aaList = new List<AminoAcid>();

            residueIndex = 1;
            foreach (var residue in seqStr)
            {
                var aa = aaSet.GetAminoAcid(residue);
                if (modificaitonInfo.ContainsKey(residueIndex))
                    aaList.Add(new ModifiedAminoAcid(aa, modificaitonInfo[residueIndex]));
                else
                    aaList.Add(aa);

                residueIndex++;
            }
            Sequence sequence = new Sequence(aaList);

            return sequence;
        }
    }

}
