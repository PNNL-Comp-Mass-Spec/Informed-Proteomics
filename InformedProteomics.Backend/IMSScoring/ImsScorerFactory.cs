using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.IMS;

namespace InformedProteomics.Backend.IMSScoring
{
    public class ImsScorerFactory
    {
        private readonly SubScoreFactory _subScoreFactory;

        public ImsScorerFactory(string paramFilePath)
        {
            _subScoreFactory = new SubScoreFactory(paramFilePath);
        }

        public ImsScorer GetImsScorer(ImsDataCached imsData, Ion precursorIon)
        {
            return new ImsScorer(imsData, precursorIon, _subScoreFactory);
        }
    }
}
