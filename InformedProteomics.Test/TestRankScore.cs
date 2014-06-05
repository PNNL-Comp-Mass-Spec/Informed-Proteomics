using System;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.LikelihoodScoring.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestRankScore
    {
        [Test]
        public void RankScore()
        {
            const int ranks = 20;
            var rankScorer = new RankScore(@"\\protoapps\UserData\Wilkins\DIA\DIA.txt");
            for (int charge = 1; charge < 4; charge++)
            {
                var ionTypes = rankScorer.GetIonTypes(charge);
                foreach (var ionType in ionTypes)
                {
                    for (int r = 0; r <= ranks; r++)
                    {
                        Console.WriteLine("Charge: {0}, Ion Type: {1}, Rank: {2}, Score: {3}",
                            charge, ionType.Name, r, rankScorer.GetScore(ionType, r, charge));
                    }
                }
            }
        }

        [Test]
        public void DiaRankScore()
        {
            const string dataFile =
    @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\raw\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
            const string tsvFile =
                @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\tsv\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.tsv";
            var parser = new TsvFileParser(tsvFile);
            var sequences = parser.GetData("Peptide");
            var charges = parser.GetData("Charge");
            var scans = parser.GetData("ScanNum");

            var lcms = LcMsRun.GetLcMsRun(dataFile, MassSpecDataType.XCaliburRun, 0, 0);
            var rankScorer =
    new DiaRankScore(
        @"C:\Users\wilk011\Documents\DataFiles\TestFolder\HCD_QExactive_Tryp.txt");

            using (
                var outFile = new StreamWriter(@"C:\Users\wilk011\Documents\DataFiles\TestFolder\HCD_QCShew_Score_2.txt"))
            {
                outFile.WriteLine("Target\tDecoy");
                for (int i = 0; i < sequences.Count; i++)
                {
                    string sequenceStr = sequences[i];
                    int charge = Convert.ToInt32(charges[i]);
                    int scan = Convert.ToInt32(scans[i]);

                    var sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(sequenceStr);
                    var decoySeq = Sequence.GetSequenceFromMsGfPlusPeptideStr(sequenceStr);
                    decoySeq.Reverse();
                    var decoyStr = decoySeq.Aggregate("", (current, aa) => current + aa);
                    decoyStr = SimpleStringProcessing.Mutate(decoyStr, sequence.Count/2);
                    decoySeq = Sequence.GetSequenceFromMsGfPlusPeptideStr(decoyStr);

                    var sequenceScore = rankScorer.GetScore(sequence, charge, scan, lcms);
                    var decoyScore = rankScorer.GetScore(decoySeq, charge, scan, lcms);
                    outFile.WriteLine("{0}\t{1}", sequenceScore, decoyScore);
                }
            }
        }
    }
}
