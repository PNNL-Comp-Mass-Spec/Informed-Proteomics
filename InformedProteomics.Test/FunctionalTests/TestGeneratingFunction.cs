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
using InformedProteomics.Scoring.GeneratingFunction;
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
            const int scanNum = 5927;
            
            if (!File.Exists(rawFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFile);
                return;
            }

            var spectralVector = GetSpectralVectorForTest(rawFile, scanNum);
            var aminoAcidsArray = AminoAcid.StandardAminoAcidArr;

            Console.WriteLine(spectralVector.Sum());


            var peaks = GetMassListForTest();
        }


        public List<DeconvolutedPeak> GetMassListForTest()
        {
            var peaks = new List<DeconvolutedPeak>();
            const string protSequence = "MNKSELIEKIASGADISKAAAGRALDSFIAAVTEGLKEGDKISLVGFGTFEVRERAERTGRNPQTGEEIKIAAAKIPAFKAGKALKDAVN";
            var sequence = new Sequence(protSequence, new AminoAcidSet());
            var proteinMass = sequence.Composition.Mass;

            foreach (var a in sequence.GetPrefixCompositions())
            {
                var peak = new DeconvolutedPeak(new Peak(a.Mass*0.5, 1), a.Mass, 2);
                peaks.Add(peak);
            }

            return peaks;
        }


        [Test]
        public void TestEffectiveMassBinning()
        {
            var aaSet = AminoAcid.StandardAminoAcidArr;
            var massComparer = new ProteinMassComparerWithBinning(aaSet);


            //Console.WriteLine(massComparer.GetBinNumber(50000));
            const string protSequence = "MNKSELIEKIASGADISKAAAGRALDSFIAAVTEGLKEGDKISLVGFGTFEVRERAERTGRNPQTGEEIKIAAAKIPAFKAGKALKDAVN";
            var sequence = new Sequence(protSequence, new AminoAcidSet());
            var proteinMass = sequence.Composition.Mass;

            var proteinMassIndex = massComparer.GetBinNumber(proteinMass);
            var nEdges = 0;
            for (var i = 0; i <= proteinMassIndex; i++)
            {
                var nodeMass = massComparer.GetMass(i);
                foreach (var aa in aaSet)
                {
                    var j = massComparer.GetBinNumber(nodeMass + aa.Mass);
                    //Console.WriteLine("{0}-{1} are connected by edge", i, j);
                    nEdges++;

                }
            }
            Console.WriteLine(@"# of nodes = {0}, # of edges = {1}", proteinMassIndex + 1, nEdges);
            
            /*
            const string protSequence = "MNKSELIEKIASGADISKAAAGRALDSFIAAVTEGLKEGDKISLVGFGTFEVRERAERTGRNPQTGEEIKIAAAKIPAFKAGKALKDAVN";
            
            var sequence = new Sequence(protSequence, new AminoAcidSet());
            var proteinMass = sequence.Composition.Mass;

            proteinMass = 20000;
            
            var comparer = new MzComparerWithBinning();
            var maxBinIndex = comparer.GetBinNumber(proteinMass);
            var minBinIndex = comparer.GetBinNumber(0);
            var nBins = maxBinIndex + 1;

            var aaSet = AminoAcid.StandardAminoAcidArr;
            var bins = new byte[nBins];
            bins[0] = 1;

            for (var i = minBinIndex; i < maxBinIndex; i++)
            {
                var mass = comparer.GetMzAverage(i);
                var minMass = comparer.GetMzStart(i);
                var maxMass = comparer.GetMzEnd(i);

                if (bins[i] == 1)
                {
                    foreach (var aa in aaSet)
                    {
                        //var j = comparer.GetBinNumber(mass + aa.Mass);

                        for (var j = comparer.GetBinNumber(minMass + aa.Mass);
                            j <= comparer.GetBinNumber(maxMass + aa.Mass);
                            j++)
                        {
                            if (j > maxBinIndex) continue;
                            if (bins[j] == 0) bins[j] = 1;                            
                        }
                    }
                }
            }
            Console.WriteLine("{0} -> {1}", nBins, bins.Count(b => b == 1));*/
        }


        public int[] GetSpectralVectorForTest(string rawFile, int scanNum)
        {
            const string protSequence = "MNKSELIEKIASGADISKAAAGRALDSFIAAVTEGLKEGDKISLVGFGTFEVRERAERTGRNPQTGEEIKIAAAKIPAFKAGKALKDAVN";
            //const string rawFile = @"\\protoapps\UserData\Jungkap\TrainingSet\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            //const int scanNum = 5927;

            var sequence = new Sequence(protSequence, new AminoAcidSet());
            var proteinMass = sequence.Composition.Mass;
            
            var nominalProteinMass = Constants.GetBinNum(proteinMass);
            var spectralVec = new int[nominalProteinMass + 1];

            foreach (var a in sequence.GetPrefixCompositions())
            {
                var binIndex = Constants.GetBinNum(a.Mass);
                spectralVec[binIndex] = 1;
            }
            
            
            const int maxCharge = 15;
            const int minCharge = 1;
            const double filteringWindowSize = 1.1;
            const int isotopeOffsetTolerance = 2;
            var tolerance = new Tolerance(10);


            var run = PbfLcMsRun.GetLcMsRun(rawFile);
            var spectrum = run.GetSpectrum(scanNum);

            // transform spectrum to spectral vector 
            var deconvolutedPeaks = Deconvoluter.GetDeconvolutedPeaks(spectrum, minCharge, maxCharge,
                isotopeOffsetTolerance, filteringWindowSize, tolerance, 0.5, 0.5);
            //var nominalProteinMass = Constants.GetBinNumHighPrecision(proteinMass);


            foreach (var peak in deconvolutedPeaks)
            {
                //var binIndex = Constants.GetBinNumHighPrecision(peak.Mass);
                var binIndex = Constants.GetBinNum(peak.Mass);
                if (binIndex >= spectralVec.Length) break;
                spectralVec[binIndex] = 1;
            }

            
            return spectralVec;
        }

    }



}
