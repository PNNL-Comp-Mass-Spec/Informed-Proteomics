using System;
using System.Collections;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class Sequence : List<AminoAcid>, IMatter, IEnumerable<Composition>
    {
        public Sequence(IEnumerable<AminoAcid> aaArr)
        {
            Composition composition = Composition.Zero;

            foreach (AminoAcid aa in aaArr)
            {
                Add(aa);
                composition += aa.Composition;
            }

            this.Composition = composition;
        }

		public Sequence (Composition composition, string sequence)
		{
			this.Composition = composition;
			this.SequenceString = sequence;
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
            from = Math.Max(from, 0);
            to = Math.Min(to, Count);
            double sum = 0;
            for (var i = from; i < to; i++)
                sum += this[i].GetMass();
            return sum;
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
