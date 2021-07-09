namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Interface for Fragment scorer factories
    /// </summary>
    public interface IFragmentScorerFactory
    {
        /// <summary>
        /// Get a scorer for the provided spectrum
        /// </summary>
        /// <param name="spectrum"></param>
        /// <param name="precursorMass"></param>
        /// <param name="precursorCharge"></param>
        /// <param name="activationMethod"></param>
        /// <returns>Scorer object</returns>
        IScorer GetScorer(ProductSpectrum spectrum, double precursorMass, int precursorCharge, ActivationMethod activationMethod = ActivationMethod.Unknown);

        /// <summary>
        /// Get a scorer for the provided scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="precursorMass"></param>
        /// <param name="precursorCharge"></param>
        /// <param name="activationMethod"></param>
        /// <returns>Scorer object</returns>
        IScorer GetScorer(int scanNum, double precursorMass = 0.0, int precursorCharge = 0, ActivationMethod activationMethod = ActivationMethod.Unknown);
    }
}
