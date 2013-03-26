
namespace InformedProteomics.Backend.Scoring
{
    class FragmentSpectrumScorer
    {
        public FragmentSpectrum FragmentSpectrum { get; private set; }
        public FragmentSpectrumParameter Parameter { get; private set; }
        public float Score { get; private set; }
       
        public FragmentSpectrumScorer(FragmentSpectrum fragmentSpectrum, FragmentSpectrumParameter par)
        {
            FragmentSpectrum = fragmentSpectrum;
            Parameter = par;
            Score = GetScore();
        }

        private float GetScore()
        {
            return SubScoreFactory.GetProductIonSpectrumScore(FragmentSpectrum, Parameter);
        }
    }
}
