using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Interface for LcMsRun data objects
    /// </summary>
    public interface ILcMsRun : IChromatogramExtractor, ISpectrumAccessor, ISpectrumExtractor
    {
        /// <summary>
        /// True if the dataset is DIA data
        /// </summary>
        bool IsDia { get; }

        /// <summary>
        /// Gets scan numbers of the fragmentation spectra whose isolation window contains the precursor ion specified
        /// </summary>
        /// <param name="precursorIon"></param>
        /// <returns>scan numbers of fragmentation spectra</returns>
        int[] GetFragmentationSpectraScanNums(Ion precursorIon);

        /// <summary>
        /// Gets scan numbers of the fragmentation spectra whose isolation window contains the precursor ion specified
        /// </summary>
        /// <param name="mostAbundantIsotopeMz"></param>
        /// <returns>scan numbers of fragmentation spectra</returns>
        int[] GetFragmentationSpectraScanNums(double mostAbundantIsotopeMz);
    }
}
