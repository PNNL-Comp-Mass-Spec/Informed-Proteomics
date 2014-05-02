using System;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Scoring.LikelihoodScoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestSpectrumScore
    {
        [Test]
        public void SpectrumScore()
        {
            const int scanNum = 2189;
            const string trainingPath = @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\";
            const string path = @"\\protoapps\UserData\Wilkins\BottomUp\HCD_Velos_TMT10_Charles\";
            const string rawFile = path + @"raw\Leishmania_TMT_test_01_01Jan14_Frodo_13-06-25.raw";
            const string trainingSet = trainingPath + "HCD_QCShew_RankProbabilities_Charge3.txt";
            const string massErrorSet = trainingPath + "HCD_QCShew_MassErrorProbabilities_Charge3.txt";
            const string peptide = "+229.163ERPQQDATAAAAAAR";
            const string decoy = "GAVFLGEEVVLC+57.021QSQESGFEDPC+57.021KTDESGTR";
            var tolerance = new Tolerance(0.5, ToleranceUnit.Th);

            var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0);

            var spectrum = lcms.GetSpectrum(scanNum);

            var spectrumScorer = new SpectrumScore(spectrum, trainingSet, massErrorSet);

            Console.WriteLine("Target score: " + spectrumScorer.GetPeptideScore(peptide, tolerance));
            Console.WriteLine("Decoy score: " + spectrumScorer.GetPeptideScore(decoy, tolerance));
        }
    }
}
