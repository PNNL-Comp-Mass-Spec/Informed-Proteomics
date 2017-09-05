using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Interface for Scorers
    /// </summary>
    public interface IScorer
    {
        /// <summary>
        /// Get the fragment score for the provided data
        /// </summary>
        /// <param name="prefixFragmentComposition"></param>
        /// <param name="suffixFragmentComposition"></param>
        /// <param name="nTerminalResidue"></param>
        /// <param name="cTerminalResidue"></param>
        /// <returns></returns>
        double GetFragmentScore(Composition.Composition prefixFragmentComposition, Composition.Composition suffixFragmentComposition, AminoAcid nTerminalResidue = null, AminoAcid cTerminalResidue = null);
    }
}
