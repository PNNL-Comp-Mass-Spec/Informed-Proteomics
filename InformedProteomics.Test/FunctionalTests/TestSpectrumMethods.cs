using System;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    class TestSpectrumMethods
    {
        [Test]
        public void TestGetAllIsotopePeaks()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string specFilePath = @"H:\Research\GlycoTopDown\raw\User_sample_test_02252015.raw";
            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            //const int scanNum = 17338;
            const double relativeIntensity = 0.1;
            var run = PbfLcMsRun.GetLcMsRun(specFilePath);
            var spec = run.GetSpectrum(17338);

            var comp = Composition.Parse("C(610) H(945) N(172) O(189) S(3)");
            var ion = new Ion(comp + BaseIonType.B.OffsetComposition, 9);   // b127(9+)
            Console.WriteLine("Composition: " + comp + " " + comp.Mass);
            Console.WriteLine("b127(9+): " + ion.GetMonoIsotopicMz());
            Console.WriteLine("b127(9+) 0th isotope: " + ion.GetIsotopeMz(0));
            Console.WriteLine("b127(9+) 6th isotope: " + ion.GetIsotopeMz(6));

            var peaks = spec.GetAllIsotopePeaks(ion, new Tolerance(10), relativeIntensity);
            var isotopes = ion.GetIsotopes(relativeIntensity).ToArray();

            for (var i = 0; i < isotopes.Length; i++)
            {
                if (peaks[i] == null) continue;
                var isotopeIndex = isotopes[i].Index;
                Console.WriteLine("{0}\t{1}\t{2}\t{3}", isotopeIndex, peaks[isotopeIndex].Mz, ion.GetIsotopeMz(isotopeIndex), GetPeakPpmError(peaks[isotopeIndex], ion.GetIsotopeMz(isotopeIndex)));
            }
        }

        public static double GetPeakPpmError(Peak peak, double theoMz)
        {
            return (peak.Mz - theoMz) / peak.Mz * Math.Pow(10, 6);
        }
    }
}
