using InformedProteomics.DIA.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestBottomUpScoring
    {
        [Test]
        public void TestParsingScoringParamFile()
        {
            const string paramFilePath = @"..\..\..\TestFiles\HCD_QExactive_Tryp.param";
            var scorer = new OldRankScorer(paramFilePath, true);
        }
    }
}
