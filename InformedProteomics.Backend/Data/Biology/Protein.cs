using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Science
{
    public class Protein : Sequence.Sequence
    {
        public Protein(IEnumerable<AminoAcid> aaArr) : base(aaArr)
        {
        }
    }
}
