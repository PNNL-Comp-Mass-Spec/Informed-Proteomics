using System.Collections.Generic;

namespace InformedProteomics.Backend.Utils
{
    public class Protein : Sequence
    {
        public Protein(IEnumerable<AminoAcid> aaArr) : base(aaArr)
        {
        }
    }
}
