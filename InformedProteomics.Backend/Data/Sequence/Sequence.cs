using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class Sequence : List<AminoAcid>, IMolecule
    {
        public Sequence(IEnumerable<AminoAcid> aaArr)
        {
            Composition composition = Composition.Zero;
            SequenceString = "";

            foreach (AminoAcid aa in aaArr)
            {
                Add(aa);
                composition += aa.Composition;
                SequenceString += aa.Residue; // added by Kyowon jeong
            }

            Composition = composition + Composition.H2O; // added by Kyowon Jeong
        }

		/*public Sequence (Composition composition, string sequence)
		{
			this.Composition = composition;
			this.SequenceString = sequence;
		}*/

        // fixed by Kyowon so that both constructors have amino acid list.
        public Sequence(Composition composition, string sequence, AminoAcidSet aminoAcidSet)
        {
            Composition = composition;
            SequenceString = sequence;
            foreach(var residue in SequenceString)
            {
                Add(aminoAcidSet.GetAminoAcid(residue));
            }
        }

        public Composition Composition { get; private set; }
    	public string SequenceString { get; set; }

        public double GetMass()
        {
            return Composition.GetMass();
        }

        public Composition GetComposition()
        {
            return Composition;
        }

        //added by kyowon jeong
        public double GetMass(int from, int to)
        {
            return GetComposition(from, to).GetMass();
        }

        //added by kyowon jeong
        public Composition GetComposition(int from, int to)
        {
            from = Math.Max(from, 0);
            to = Math.Min(to, Count);
            var composition = Composition.Zero;
            for (var i = from; i < to; i++)
                composition += this[i].Composition;
            return composition;
        }
    }
}
