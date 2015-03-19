using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Scoring;

using NUnit.Framework;


namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestProMex
    {
        [Test]
        public void TestPredictPTMfromMs1ft()
        {
            const string resultFilePath = @"\\protoapps\UserData\Jungkap\FeatureFinding\ProMex_v1.1\test\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv";
            const string ms1ftFilePath = @"\\protoapps\UserData\Jungkap\FeatureFinding\ProMex_v1.1\test\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";

            var parser = new TsvFileParser(resultFilePath);
            var sequences = parser.GetData("Sequence");
            var modifications = parser.GetData("Modifications");
            var scanNums = parser.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
            var qValues = parser.GetData("QValue").Select(s => Convert.ToDouble(s)).ToArray();
            var nMacthed = parser.GetData("#MatchedFragments");
            


            var aaSet = new AminoAcidSet();



            var ptmList = new List<Tuple<int, double, double>>();

            for (var i = 0; i < parser.NumData; i++)
            {
                if (qValues[i] > 0.01) continue;
                //var sequenceComp = aaSet.GetComposition(sequences[i]) + Composition.H2O;
                var seq = new Sequence(sequences[i], aaSet);
                var sequenceComp = seq.Composition + Composition.H2O;

                var modComposition = Composition.Zero;
                var modsStr = modifications[i];
                if (modsStr.Length == 0) continue;
                var mods = modsStr.Split(',');
                foreach (var modStr in mods.Where(str => str.Length > 0))
                {
                    var modName = modStr.Split()[0];
                    var mod = Modification.Get(modName);
                    modComposition += mod.Composition;
                }
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", scanNums[i], sequenceComp.Mass, modComposition.Mass, nMacthed[i], sequences[i], modsStr);
                //var compFromSeqAndMods = sequenceComp + modComposition;
                //Assert.True(compFromSeqAndMods.Equals(compositions[i]));

                ptmList.Add(new Tuple<int, double, double>(scanNums[i], sequenceComp.Mass, modComposition.Mass));
            }

            //var featureParser = new TsvFileParser(ms1ftFilePath);
            //var minScan = featureParser.GetData("MinScan").Select(s => Convert.ToInt32(s)).ToArray();
            //var maxScan = featureParser.GetData("MaxScan").Select(s => Convert.ToInt32(s)).ToArray();
            //var monoMass = featureParser.GetData("MonoMass").Select(s => Convert.ToDouble(s)).ToArray();
        }
        
        [Test]
        public void TestGeneratingMs1FeatureFile()
        {
            const string specFilePath = @"\\protoapps\UserData\Sangtae\TestData\SpecFiles\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.raw";
            
            const int minScanCharge = 2;
            const int maxScanCharge = 60;
            const double minScanMass = 3000;
            const double maxScanMass = 3100;
            const int maxThreads = 10;

            var run = PbfLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 0, 0);
            IMs1FeatureExtract extractor = new Ms1FeatureMatrix(run, minScanCharge, maxScanCharge, maxThreads);
            var outputFilePath = extractor.GetFeatureFile(specFilePath, minScanMass, maxScanMass);
        }

        [Test]
        public void TestProMexFilter()
        {
            const string specFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            var run = PbfLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 0, 0);
            const string ms1FtPath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\ProductionQCShew\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";
            var filter = new Ms1FtFilter(run, new Tolerance(10), ms1FtPath, 0.15);

//            Console.WriteLine("ScanNums: {0}", string.Join("\t",filter.GetMatchingMs2ScanNums(8480.327609)));
            Assert.IsTrue(filter.GetMatchingMs2ScanNums(8480.327609).Contains(5255));
        }
    }
}
