using System;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestFitScore
    {
        [Test]
        public void TestFitScoreCalculation()
        {
            var run = LcMsRun.GetLcMsRun(TestLcMsRun.TestTopDownRawFilePath, MassSpecDataType.XCaliburRun);
            var spec = run.GetSpectrum(810) as ProductSpectrum;
            Assert.True(spec != null);

            const string suf54 = "ENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK";
            var suf54Comp = new AminoAcidSet().GetComposition(suf54);
            Assert.True(suf54Comp != null);

            var z6Type = new IonTypeFactory(10).GetIonType("z6");
            var z54C6 = z6Type.GetIon(suf54Comp);
            Console.WriteLine("MonoMz: {0}, MonoMass: {1}", z54C6.GetMonoIsotopicMz(), z54C6.Composition.Mass);
            Console.WriteLine("FitScore: {0}", spec.GetFitScore(z54C6, new Tolerance(15), 0.1));
        }

        [Test]
        public void TestFitScoreComputationTime()
        {
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

            Console.WriteLine("Calculating fit scores {0} times", numTrials);

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();

            for (var trial = 0; trial < numTrials; trial++)
            {
                FitScoreCalculator.GetDeconToolsFit(theoretical[trial], observed[trial]);
                //FitScoreCalculator.GetFitOfNormalizedVectors(theoretical[trial], observed[trial]);
                //FitScoreCalculator.GetCosine(theoretical[trial], observed[trial]);
                //FitScoreCalculator.GetFitNormalizedByTheoMaxIsotope(theoretical[trial], observed[trial]);
            }

            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }
    }
}
