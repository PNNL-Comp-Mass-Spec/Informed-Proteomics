namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface IFragmentScorerFactory
    {
        IScorer GetScorer(ProductSpectrum spectrum, double precursorMass, int precursorCharge, ActivationMethod activationMethod = ActivationMethod.Unknown);

        IScorer GetScorer(int scanNum, double precursorMass = 0.0, int precursorCharge = 0, ActivationMethod activationMethod = ActivationMethod.Unknown);
    }
}
