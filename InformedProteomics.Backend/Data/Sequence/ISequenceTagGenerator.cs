using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// Interface of Sequence Tag Generation implementations
    /// </summary>
    public interface ISequenceTagGenerator
    {
        /// <summary>
        /// Get the sequence tags for ProductSpectrum <paramref name="spec"/>
        /// </summary>
        /// <param name="spec"></param>
        /// <returns></returns>
        IEnumerable<string> GetSequenceTags(ProductSpectrum spec);
    }
}
