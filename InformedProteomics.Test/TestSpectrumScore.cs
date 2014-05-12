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
            const int scanNum = 22383;
            const string trainingPath = @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QE_TMT10_Charles\";
            const string path = @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QE_TMT10_Charles\";
            const string rawFile = path + @"raw\Leishmania_TMT_NiNTA_FT_01_R2_20Jan14_Samwise_13-07-28.raw";
            const string trainingSet = trainingPath + "HCD_QE_TMT10_Charles_RankProbabilities_Charge3.txt";
            const string massErrorSet = trainingPath + "HCD_QE_TMT10_Charles_MassErrorProbabilities_Charge3.txt";
            const string peptide = "+229.163SAEALQQLNDASFNASATSAAQAAADVYK+229.163";
            const string decoy = "+229.163AC+57.021C+57.021GGDGYAVGVEPAAVEK+229.163";
            var tolerance = new Tolerance(0.5, ToleranceUnit.Th);

            var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0);

            var spectrum = lcms.GetSpectrum(scanNum);

            var spectrumScorer = new SpectrumScore(spectrum, tolerance, trainingSet, massErrorSet);

            Console.WriteLine("Target score: " + spectrumScorer.GetPeptideScore(peptide));
            Console.WriteLine("Decoy score: " + spectrumScorer.GetPeptideScore(decoy));
        }
    }
}
