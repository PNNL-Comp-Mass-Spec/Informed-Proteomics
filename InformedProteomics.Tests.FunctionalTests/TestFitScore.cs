using System;
using System.IO;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.FunctionalTests
{
    [TestFixture]
    public class TestFitScore
    {
        [Test]
        public void TestFitScoreCalculationCid()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            if (!File.Exists(FilePaths.TestTopDownRawFilePathCid))
            {
                Assert.Ignore(@"Skipping test " + methodName + @" since file not found: " + FilePaths.TestTopDownRawFilePathCid);
            }

            var run = InMemoryLcMsRun.GetLcMsRunScanRange(FilePaths.TestTopDownRawFilePathCid, 5743, 5743);
            var spec = run.GetSpectrum(5743);
            Assert.True(spec != null);

            const string protein = "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGIVVDYVLEFDVPDELIVDRIVGRRVHAASGRVYHVKFNPPKVEGKDDVTGEDLTTRKDDQEETVRKRLVEYHQMTAPLIGYYQKEAEAGNTKYAKVDGTQAVADVRAALEKILG";
            var protComp = new AminoAcidSet().GetComposition(protein) + Composition.H2O;
            Assert.True(protComp != null);
            Assert.True(protComp.C == 1035);
            Assert.True(protComp.H == 1683);
            Assert.True(protComp.N == 289);
            Assert.True(protComp.O == 318);
            Assert.True(protComp.P == 0);
            Assert.True(protComp.S == 7);
            Assert.True(Math.Abs(protComp.Mass - 23473.245267145) < 0.0000001);
            Assert.True(protComp.NominalMass == 23461);

            var ion = new Ion(protComp, 20);
//            ion.Composition.ComputeApproximateIsotopomerEnvelop();
            var isotopomerEnvelop = ion.Composition.GetIsotopomerEnvelopeRelativeIntensities();
            Console.WriteLine(@"MonoMz: {0}, MonoMass: {1}", ion.GetMonoIsotopicMz(), ion.Composition.Mass);

            var matchedPeaks = spec.GetAllIsotopePeaks(ion, new Tolerance(15), 0.1);
            for (var i = 0; i < matchedPeaks.Length; i++)
            {
                Console.WriteLine(@"{0}\t{1}\t{2}\t{3}", i, ion.GetIsotopeMz(i), isotopomerEnvelop[i], matchedPeaks[i] == null ? 0 : matchedPeaks[i].Intensity);
            }
            var fitScore = spec.GetFitScore(ion, new Tolerance(15), 0.1);
            var cosine = spec.GetConsineScore(ion, new Tolerance(15), 0.1);
            var corr = spec.GetCorrScore(ion, new Tolerance(15), 0.1);

            Console.WriteLine(@"FitScore: {0}", fitScore);
            Console.WriteLine(@"Cosine: {0}", cosine);
            Console.WriteLine(@"Corr: {0}", corr);

            Assert.True(Math.Abs(fitScore - 0.181194589537041) < 0.0001);
            Assert.True(Math.Abs(cosine - 0.917609346566222) < 0.0001);
            Assert.True(Math.Abs(corr - 0.808326778009839) < 0.0001);
        }

        [Test]
        public void TestFitScoreCalculationEtd()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            if (!File.Exists(FilePaths.TestTopDownRawFilePathEtd))
            {
                Assert.Ignore(@"Skipping test " + methodName + @" since file not found: " + FilePaths.TestTopDownRawFilePathCid);
            }

            var run = InMemoryLcMsRun.GetLcMsRunScanRange(FilePaths.TestTopDownRawFilePathEtd, 810, 810);
            var spec = run.GetSpectrum(810) as ProductSpectrum;
            Assert.True(spec != null);

            const string suf54 = "ENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK";
            var suf54Comp = new AminoAcidSet().GetComposition(suf54);
            Assert.True(suf54Comp != null);

            var ionType = new IonTypeFactory(10).GetIonType("z6");
            var ion = ionType.GetIon(suf54Comp);
            //ion.Composition.ComputeApproximateIsotopomerEnvelop();
            Console.WriteLine("MonoMz: {0}, MonoMass: {1}", ion.GetMonoIsotopicMz(), ion.Composition.Mass);

            var fitScore = spec.GetFitScore(ion, new Tolerance(15), 0.1);
            Console.WriteLine("FitScore: {0}", fitScore);
            Assert.True(fitScore < 0.15);
        }

        [Test]
        public void TestFitScoreComputationTime()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const int numTrials = 1000000;
            const int numIsotopes = 20;

            var random = new Random();
            var theoretical = new double[numTrials][];
            var observed = new double[numTrials][];
            for (var i = 0; i < numTrials; i++)
            {
                theoretical[i] = new double[numIsotopes];
                observed[i] = new double[numIsotopes];
                for (var j = 0; j < numIsotopes; j++)
                {
                    theoretical[i][j] = random.NextDouble();
                    observed[i][j] = random.NextDouble();
                }
            }

            Console.WriteLine("Calculating fit scores {0} times", numTrials * 2);

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();

            for (var trial = 0; trial < numTrials; trial++)
            {
                //FitScoreCalculator.GetDeconToolsFit(theoretical[trial], observed[trial]);
                //FitScoreCalculator.GetFitOfNormalizedVectors(theoretical[trial], observed[trial]);
                //FitScoreCalculator.GetCosine(theoretical[trial], observed[trial]);
                //FitScoreCalculator.GetFitNormalizedByTheoMaxIsotope(theoretical[trial], observed[trial]);
                FitScoreCalculator.GetPearsonCorrelation(theoretical[trial], observed[trial]);
            }

            var series1 = new double[7];
            series1[0] = 0.2;
            series1[1] = 0.5;
            series1[2] = 0.8;
            series1[3] = 0.7;
            series1[4] = 0.6;
            series1[5] = 0.3;
            series1[6] = 0.05;

            var series2 = new double[series1.Length];
            double smallestFit = 1;
            var threshold = (int)(1000000 / 10.0);

            for (var iteration = 0; iteration < 1000000; iteration++)
            {
                for (var i = 0; i < series1.Length; i++)
                {
                    series2[i] = series1[i] + random.NextDouble() / 6 - (1 / 12.0);
                }

                var result = FitScoreCalculator.GetPearsonCorrelation(series1, series2);

                if (iteration % threshold == 0)
                    Console.WriteLine(@"Fit=" + result);

                if (result < smallestFit)
                    smallestFit = result;
            }

            Console.WriteLine(@"SmallestFit=" + smallestFit);
            Assert.IsTrue(smallestFit > 0.94);

            sw.Stop();

            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sw.Elapsed.TotalSeconds);
        }
    }
}
