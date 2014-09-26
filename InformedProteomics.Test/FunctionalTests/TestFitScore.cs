using System;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
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
        public void TestFitScoreCalculationCid()
        {
            var run = InMemoryLcMsRun.GetLcMsRun(TestLcMsRun.TestTopDownRawFilePathCid, MassSpecDataType.XCaliburRun);
            var spec = run.GetSpectrum(5743);
            Assert.True(spec != null);

            const string protein = "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGIVVDYVLEFDVPDELIVDRIVGRRVHAASGRVYHVKFNPPKVEGKDDVTGEDLTTRKDDQEETVRKRLVEYHQMTAPLIGYYQKEAEAGNTKYAKVDGTQAVADVRAALEKILG";
            var protComp = new AminoAcidSet().GetComposition(protein) + Composition.H2O;
            Assert.True(protComp != null);

            var ion = new Ion(protComp, 20);
//            ion.Composition.ComputeApproximateIsotopomerEnvelop();
            var isotopomerEnvelop = ion.Composition.GetIsotopomerEnvelope();
            Console.WriteLine("MonoMz: {0}, MonoMass: {1}", ion.GetMonoIsotopicMz(), ion.Composition.Mass);

            var matchedPeaks = spec.GetAllIsotopePeaks(ion, new Tolerance(15), 0.1);
            for (var i = 0; i < matchedPeaks.Length; i++)
            {
                Console.WriteLine("{0}\t{1}\t{2}\t{3}", i, ion.GetIsotopeMz(i), isotopomerEnvelop[i], matchedPeaks[i] == null ? 0 : matchedPeaks[i].Intensity);
            }
            var fitScore = spec.GetFitScore(ion, new Tolerance(15), 0.1);
            Console.WriteLine("FitScore: {0}", fitScore);
            Console.WriteLine("Cosine: {0}", spec.GetConsineScore(ion, new Tolerance(15), 0.1));
            Console.WriteLine("Corr: {0}", spec.GetCorrScore(ion, new Tolerance(15), 0.1));
            //            Assert.True(fitScore < 0.15);
        }

        [Test]
        public void TestFitScoreCalculationEtd()
        {
            var run = InMemoryLcMsRun.GetLcMsRun(TestLcMsRun.TestTopDownRawFilePathEtd, MassSpecDataType.XCaliburRun);
            var spec = run.GetSpectrum(810) as ProductSpectrum;
            Assert.True(spec != null);

            const string suf54 = "ENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK";
            var suf54Comp = new AminoAcidSet().GetComposition(suf54);
            Assert.True(suf54Comp != null);

            var ionType = new IonTypeFactory(10).GetIonType("z6");
            var ion = ionType.GetIon(suf54Comp);
            ion.Composition.ComputeApproximateIsotopomerEnvelop();
            Console.WriteLine("MonoMz: {0}, MonoMass: {1}", ion.GetMonoIsotopicMz(), ion.Composition.Mass);

            var fitScore = spec.GetFitScore(ion, new Tolerance(15), 0.1);
            Console.WriteLine("FitScore: {0}", fitScore);
            Assert.True(fitScore < 0.15);
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
                //FitScoreCalculator.GetDeconToolsFit(theoretical[trial], observed[trial]);
                //FitScoreCalculator.GetFitOfNormalizedVectors(theoretical[trial], observed[trial]);
                //FitScoreCalculator.GetCosine(theoretical[trial], observed[trial]);
                //FitScoreCalculator.GetFitNormalizedByTheoMaxIsotope(theoretical[trial], observed[trial]);
                FitScoreCalculator.GetPearsonCorrelation(theoretical[trial], observed[trial]);
            }

            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }
    }
}
