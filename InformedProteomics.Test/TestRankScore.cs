using System;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.LikelihoodScoring.Scoring;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestRankScore
    {
        [Test]
        public void RankScoreParamResources()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const int ranks = 20;
            // Ms2DetectorType.IonTrap
            // Protocol.Standard
            var rankScorer = new RankScore(ActivationMethod.HCD, Enzyme.Trypsin);
            for (var charge = 1; charge < 4; charge++)
            {
                var ionTypes = rankScorer.GetIonTypes(charge, 0);
                foreach (var ionType in ionTypes)
                {
                    for (var r = 0; r <= ranks; r++)
                    {
                        if (r < 4 || r > ranks - 4)
                        {
                            Console.WriteLine(@"Charge: {0}, Ion Type: {1}, Rank: {2}, Score: {3:F4}",
                                charge, ionType.Name, r, rankScorer.GetScore(ionType, r, charge, 0.0));
                        }
                        else if (r == 4)
                        {
                            Console.WriteLine("  ...");
                        }
                    }
                    Console.WriteLine();
                }
                Console.WriteLine();
            }
        }

        [Test]
        [Category("Local_Testing")]
        public void RankScore()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const int ranks = 20;

            const string filePath = @"\\protoapps\UserData\Wilkins\DIA\DIA.txt";
            if (!File.Exists(filePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, filePath);
            }

            var rankScorer = new RankScore(filePath);
            for (var charge = 1; charge < 4; charge++)
            {
                var ionTypes = rankScorer.GetIonTypes(charge, 0);
                foreach (var ionType in ionTypes)
                {
                    for (var r = 0; r <= ranks; r++)
                    {
                        Console.WriteLine("Charge: {0}, Ion Type: {1}, Rank: {2}, Score: {3}",
                            charge, ionType.Name, r, rankScorer.GetScore(ionType, r, charge, 0.0));
                    }
                }
            }
        }

        [Test]
        [Category("PNL_Domain")]
        [Category("Local_Testing")]
        public void DiaRankScore()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string dataFile = @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\raw\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
            const string tsvFile = @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\tsv\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.tsv";

            if (!File.Exists(dataFile))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, dataFile);
            }

            if (!File.Exists(tsvFile))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, tsvFile);
            }

            var parser = new TsvFileParser(tsvFile);
            var sequences = parser.GetData("Peptide");
            var charges = parser.GetData("Charge");
            var scans = parser.GetData("ScanNum");

            var lcms = InMemoryLcMsRun.GetLcMsRun(dataFile, 0, 0);
            var rankScorer = new DiaRankScore(@"C:\Users\wilk011\Documents\DataFiles\TestFolder\HCD_QExactive_Tryp.txt");

            using (var outFile = new StreamWriter(@"C:\Users\wilk011\Documents\DataFiles\TestFolder\HCD_QCShew_Score_2.txt"))
            {
                outFile.WriteLine("Target\tDecoy");
                for (var i = 0; i < sequences.Count; i++)
                {
                    var sequenceStr = sequences[i];
                    var charge = Convert.ToInt32(charges[i]);
                    var scan = Convert.ToInt32(scans[i]);

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
