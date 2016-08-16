namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface IFragmentScorerFactory
    {
        IScorer GetScorer(Spectrum spectrum, double precursorMass, int precursorCharge);
    }
}
