using System;
using System.Collections;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class Sequence : List<AminoAcid>, IMatter, IEnumerable<Composition>
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

            this.Composition = composition;
        }

		/*public Sequence (Composition composition, string sequence)
		{
			this.Composition = composition;
			this.SequenceString = sequence;
		}*/

        // fixed by Kyowon so that both constructors have amino acid list.
        public Sequence(Composition composition, string sequence, AminoAcidSet aminoAcidSet)
        {
            this.Composition = composition;
            this.SequenceString = sequence;
            foreach(var residue in SequenceString)
            {
                Add(aminoAcidSet.GetAminoAcids(residue)[0]);
            }
        }

        public Composition Composition { get; private set; }
    	public string SequenceString { get; set; }

        public double GetMass()
        {
            return Composition.GetMass();
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

        IEnumerator<Composition> IEnumerable<Composition>.GetEnumerator()
        {
            return new CompositionEnum(this);
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return (IEnumerator) new CompositionEnum(this);
        }

        private class CompositionEnum : IEnumerator<Composition>
        {
            private Sequence _sequence;
            public CompositionEnum(Sequence sequence)
            {
                _sequence = sequence;
            }

            public Composition Current
            {
                get { throw new NotImplementedException(); }
            }

            public void Dispose()
            {
                throw new NotImplementedException();
            }

            object IEnumerator.Current
            {
                get { throw new NotImplementedException(); }
            }

            public bool MoveNext()
            {
                throw new NotImplementedException();
            }

            public void Reset()
            {
                throw new NotImplementedException();
            }
        }
    }
}
