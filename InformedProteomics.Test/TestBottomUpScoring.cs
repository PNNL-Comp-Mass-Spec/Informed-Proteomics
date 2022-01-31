using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.MathAndStats;
using InformedProteomics.Backend.Utils;
using InformedProteomics.BottomUp.Scoring;
using InformedProteomics.Tests.Base;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    [Category("Local_Testing")]
    public class TestBottomUpScoring
    {
        [Test]
        [Category("Local_Testing")]
        public void TestPsm()
        {
            var methodName = MethodBase.GetCurrentMethod()?.Name;
            Utils.ShowStarting(methodName);

            const string specFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
            if (!File.Exists(specFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            const char pre = 'R';
            const string sequence = "LENWPPASLADDL";
            const char post = 'A';
            const string annotation = "R.LENWPPASLADDL._";
            const int charge = 2;
            const int ms2ScanNum = 25534;

            var aaSet = new AminoAcidSet();

            var run = InMemoryLcMsRun.GetLcMsRun(specFilePath, 0, 0);
            var ms2Scorer = new ProductScorerBasedOnDeconvolutedSpectra(run, 1, 2, 10, 0, 1.1);
            ms2Scorer.DeconvoluteAllProductSpectra();
            var scorer = ms2Scorer.GetMs2Scorer(ms2ScanNum);

            var graph = SequenceGraph.CreateGraph(aaSet, annotation);
            graph.SetSink(0);
            var score = graph.GetFragmentScore(scorer);
            Console.WriteLine("Fast search score: " + score);
            var composition = graph.GetSinkSequenceCompositionWithH2O();

            var informedScorer = new InformedBottomUpScorer(run, aaSet, 1, 15, new Tolerance(10));
            var refinedScore = informedScorer.GetScores(pre, sequence, post, composition, charge, ms2ScanNum);
            Console.WriteLine("RefinedScores: {0}", refinedScore);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestVennDiagram()
        {
            var methodName = MethodBase.GetCurrentMethod()?.Name;
            Utils.ShowStarting(methodName);

            const string result1Path = @"C:\cygwin\home\kims336\Data\QCShewQE\NoMod_NTT1.tsv";
            const string result2Path = @"C:\cygwin\home\kims336\Data\QCShewQE\Ic_NTT1_Test\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28_IcTda.tsv";

            if (!File.Exists(result1Path))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, result1Path);
            }

            if (!File.Exists(result2Path))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, result2Path);
            }

            const double pepQValueThreshold = 0.01;
            var result1 = new TsvFileParser(result1Path);
            var result2 = new TsvFileParser(result2Path);

            var vennDiagram = new VennDiagram<string>(result1.GetPeptides(pepQValueThreshold),
                                                      result2.GetPeptides(pepQValueThreshold));

            var intersectionPeptides = vennDiagram.Intersection;
            Console.WriteLine(vennDiagram.Set1 + " " + vennDiagram.Set2);
            Console.WriteLine(vennDiagram.Set1 + " " + vennDiagram.Intersection + " " + vennDiagram.Set2Only);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestLogLikelihoodScoring()
        {
            var methodName = MethodBase.GetCurrentMethod()?.Name;
            Utils.ShowStarting(methodName);

            const string specFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
            if (!File.Exists(specFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            //const string seqStr = "IAHESDDEKGHAAK";
            //var composition = Composition.Parse("C(62) H(98) N(20) O(24) S(0)");
            //const int charge = 4;
            //const int ms2ScanNum = 12901;

            var aaSet = new AminoAcidSet();

            var run = InMemoryLcMsRun.GetLcMsRun(specFilePath, 1.4826, 0);
            var scorer = new InformedBottomUpScorer(run, aaSet, 1, 2, new Tolerance(10));
            //var refinedScore = scorer.GetScores(AminoAcid.PeptideNTerm, seqStr, AminoAcid.PeptideCTerm, composition,
            //    charge, ms2ScanNum);
            //            Console.WriteLine("RefinedScores: {0}", refinedScore.Score);
        }

        [Test]
        [Category("Local_Testing")]
        public void TestPeptideLevelStats()
        {
            var methodName = MethodBase.GetCurrentMethod()?.Name;
            Utils.ShowStarting(methodName);

            const int topK = 10;
            const string resultDir = @"D:\Research\Data\UW\QExactive\Ic_NTT1_Rescoring";
            if (!Directory.Exists(resultDir))
            {
                Assert.Ignore("Skipping test {0} since folder not found: {1}", methodName, resultDir);
            }

            var concatenated = new List<string>();
            const string mzRange = "650to775";
            foreach (var specFilePath in Directory.GetFiles(resultDir, "*DIA*" + mzRange + "*IcTarget.tsv"))
            {
                concatenated.AddRange(File.ReadAllLines(specFilePath).Skip(1));
            }

            foreach (var specFilePath in Directory.GetFiles(resultDir, "*DIA*" + mzRange + "*IcDecoy.tsv"))
            {
                concatenated.AddRange(File.ReadAllLines(specFilePath).Skip(1));
            }

            const string headerStr =
                "Scan\tPre\tSequence\tPost\tModifications\t" +
                "Composition\tProteinName\tProteinDesc\tProteinLength\t" +
                "Start\tEnd\tCharge\tMostAbundantIsotopeMz\t" +
                "Mass\t#MatchedFragments\tIcScore";
            var header = headerStr.Split('\t').ToList();
            if (concatenated.Count <= 1)
            {
                return;
            }

            var scoreIndex = header.IndexOf("IcScore");
            var sequenceIndex = header.IndexOf("Sequence");
            var preIndex = header.IndexOf("Pre");
            var postIndex = header.IndexOf("Post");
            var proteinIndex = header.IndexOf("ProteinName");
            var scanIndex = header.IndexOf("Scan");
            if (scoreIndex < 0 || sequenceIndex < 0 || preIndex < 0 || postIndex < 0 || proteinIndex < 0)
            {
                return;
            }

            //var distinctSorted = concatenated.OrderByDescending(r => Convert.ToDouble(r.Split('\t')[scoreIndex]))
            //    .GroupBy(r => Convert.ToDouble(r.Split('\t')[scanIndex]))
            //    .Select(grp => grp.First())
            //    .ToArray();

            var distinctSorted = concatenated
                .GroupBy(r => Convert.ToInt32(r.Split('\t')[scanIndex]))
                .Select(g => new
                {
                    FirstTwo = g.OrderByDescending(r => Convert.ToDouble(r.Split('\t')[scoreIndex])).Take(topK)
                })
                .SelectMany(g => g.FirstTwo)
                .OrderByDescending(r => Convert.ToDouble(r.Split('\t')[scoreIndex]))
                .GroupBy(r => r.Split('\t')[preIndex] + r.Split('\t')[sequenceIndex] + r.Split('\t')[postIndex])
                .Select(grp => grp.First())
                .ToArray();

            // Calculate q values
            var numDecoy = 0;
            var numTarget = 0;
            var fdr = new double[distinctSorted.Length];
            for (var i = 0; i < distinctSorted.Length; i++)
            {
                var row = distinctSorted[i];
                var columns = row.Split('\t');
                var protein = columns[proteinIndex];
                if (protein.StartsWith(FastaDatabaseConstants.DecoyProteinPrefix))
                {
                    numDecoy++;
                }
                else
                {
                    numTarget++;
                }

                fdr[i] = numDecoy / (double)numTarget;
            }

            var numPeptides = 0;
            var qValue = new double[fdr.Length];
            qValue[fdr.Length - 1] = fdr[fdr.Length - 1];
            for (var i = fdr.Length - 2; i >= 0; i--)
            {
                qValue[i] = Math.Min(qValue[i + 1], fdr[i]);
                if (qValue[i] < 0.01)
                {
                    numPeptides++;
                }
            }

            Console.WriteLine("NumPeptides: {0}", numPeptides);
        }

        [Test]
        [Category("Local_Testing")]
        public void CompareIpaIc()
        {
            var methodName = MethodBase.GetCurrentMethod()?.Name;
            Utils.ShowStarting(methodName);

            const string resultDir = @"D:\Research\Data\UW\QExactive\Ic_NTT2_03";
            if (!Directory.Exists(resultDir))
            {
                Assert.Ignore("Skipping test {0} since folder not found: {1}", methodName, resultDir);
            }

            var targetPeptides = new HashSet<string>();
            foreach (var icResultFilePath in Directory.GetFiles(resultDir, "*DIA*IcTarget.tsv"))
            {
                var icParser = new TsvFileParser(icResultFilePath);
                foreach (var peptide in icParser.GetData("Sequence"))
                {
                    targetPeptides.Add(peptide);
                }
            }

            const string ipaResultPath = @"D:\Research\Data\UW\QExactive\DIA_All_Summary.tsv";
            if (!File.Exists(ipaResultPath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, methodName);
            }

            var parser = new TsvFileParser(ipaResultPath);
            var ipaPeptides = parser.GetPeptides(0.005).Select(p => p.Replace("C+57.021", "C"));
            var ipaOnly = 0;
            var both = 0;
            foreach (var ipaPeptide in ipaPeptides)
            {
                if (targetPeptides.Contains(ipaPeptide))
                {
                    ++both;
                }
                else
                {
                    ++ipaOnly;
                    Console.WriteLine(ipaPeptide);
                }
            }

            Console.WriteLine("Both: {0}, IpaOnly: {1}, Sum: {2}", both, ipaOnly, both + ipaOnly);
        }
    }
}
