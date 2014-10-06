using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Sequence
{
    public interface ISequenceTagGenerator
    {
        IEnumerable<string> GetSequenceTags(ProductSpectrum spec);
    }
}
