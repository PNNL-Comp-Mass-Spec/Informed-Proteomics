using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestGeneratingFunction
    {
        [Test]
        public void TestGetScoreDistribution()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);
            const string rawFile = @"\\protoapps\UserData\Jungkap\TrainingSet\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            const int scanNum = 5917;
            
            if (!File.Exists(rawFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFile);
                return;
            }

            var spectralVector = GetSpectralVectorForTest(rawFile, scanNum);
            var aminoAcidsArray = AminoAcid.StandardAminoAcidArr;

            Console.WriteLine(spectralVector.Sum());

        }


        public int[] GetSpectralVectorForTest(string rawFile, int scanNum)
        {
            const string protSequence =
                "MNKSELIEKIASGADISKAAAGRALDSFIAAVTEGLKEGDKISLVGFGTFEVRERAERTGRNPQTGEEIKIAAAKIPAFKAGKALKDAVN";
            
            const int maxCharge = 15;
            const int minCharge = 1;
            const double filteringWindowSize = 1.1;
            const int isotopeOffsetTolerance = 2;
            var tolerance = new Tolerance(10);

            const string annotation = "_." + protSequence + "._";
            var seqGraph = SequenceGraph.CreateGraph(new AminoAcidSet(), AminoAcid.ProteinNTerm, protSequence,
                AminoAcid.ProteinCTerm);
            if (seqGraph == null) return null;
            seqGraph.SetSink(0);
            var neutral = seqGraph.GetSinkSequenceCompositionWithH2O() - Composition.Hydrogen;
            var proteinMass = neutral.Mass;

            var run = PbfLcMsRun.GetLcMsRun(rawFile);
            var spectrum = run.GetSpectrum(scanNum);

            // transform spectrum to spectral vector 
            var deconvolutedPeaks = Deconvoluter.GetDeconvolutedPeaks(spectrum, minCharge, maxCharge,
                isotopeOffsetTolerance, filteringWindowSize, tolerance, 0.5, 0.5);
            var nominalProteinMass = Constants.GetBinNumHighPrecision(proteinMass);
            var spectralVec = new int[nominalProteinMass];

            foreach (var peak in deconvolutedPeaks)
            {
                var binIndex = Constants.GetBinNumHighPrecision(peak.Mass);
                if (binIndex >= spectralVec.Length) break;
                spectralVec[binIndex] = 1;
            }

            return spectralVec;
        }

    }



}
