using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using DeconTools.Backend.Core;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;
using InformedProteomics.Backend.SequenceTag;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestSequenceTagFinder
    {
        [Test]
        public void TestSequenceTag()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            //const string TestRawFile = @"D:\\Vlad_TopDown\\raw\\yufeng_column_test2.raw";
            //const string TestResultFile = @"D:\\Vlad_TopDown\\results\\yufeng_column_test2_IcTda.tsv";
            const string TestRawFile = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            const string TestResultFile = @"D:\MassSpecFiles\training\IdResult\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv";
            //const string TestRawFile = @"D:\MassSpecFiles\Lewy\Lewy_intact_01.raw";
            //const string TestResultFile = @"D:\MassSpecFiles\Lewy\Lewy_intact_01_IcTda.tsv";

            if (!File.Exists(TestRawFile))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, TestRawFile);
            }

            if (!File.Exists(TestResultFile))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, TestResultFile);
            }

            //const int MaxTags = 100000;
            var tsvParser = new TsvFileParser(TestResultFile);
            var headerList = tsvParser.GetHeaders();
            var tsvData = tsvParser.GetAllData();
            var ms2ScanNumbers = tsvData["Scan"];
        
            var run = PbfLcMsRun.GetLcMsRun(TestRawFile);
            var nSpec = 0;
            var nHitSpec = 0;

            for (var i = 0; i < ms2ScanNumbers.Count; i++)
            //foreach(var scanNum in targetScans)
            {
                var scanNum = Int32.Parse(ms2ScanNumbers[i]);

                if (scanNum != 4672) continue;
                
                var spectrum = run.GetSpectrum(scanNum) as ProductSpectrum;

                int tsvIndex = ms2ScanNumbers.FindIndex(x => Int32.Parse(x) == scanNum);
                var qValue = double.Parse(tsvData["QValue"].ElementAt(tsvIndex));
                if (qValue > 0.01) break;

                var seqStr = tsvData["Sequence"].ElementAt(tsvIndex).Trim();
                var modStr = tsvData["Modifications"].ElementAt(tsvIndex).Trim();
                var tolerance = new Tolerance(5);
                var tagFinder = new SequenceTagFinder(spectrum, tolerance);
                var nTags = 0;
                var nHit = 0;

                var seqOjb = Sequence.CreateSequence(seqStr, modStr, new AminoAcidSet());
                var compWithoutH2O = seqOjb.Composition - Composition.H2O;

                Console.WriteLine(compWithoutH2O.Mass);

                foreach (var seqTagStr in tagFinder.GetAllSequenceTagString())
                {
                    if (seqStr.Contains(seqTagStr.Sequence)) //|| seqStr.Contains(Reverse(tagStr)))
                    {

                        var idx = seqStr.IndexOf(seqTagStr.Sequence);

                        //seqStr.Substring(0, idx)
                        var comp2 = seqOjb.GetComposition(0, idx);

                        Console.Write(comp2.Mass);
                        Console.Write("\t");

                        Console.Write(seqTagStr.FlankingMass);
                        Console.Write("\t");
                        Console.Write(seqTagStr.Sequence);
                        Console.Write("\t");
                        Console.Write(seqTagStr.IsPrefix);
                        Console.WriteLine("");

                        nHit++;
                    }
                    nTags++;                    
                }

                
                nSpec++;
                if (nHit > 0) nHitSpec++;

                Console.WriteLine(@"[{0}]seqLen = {1}: {2}/{3}", scanNum, seqStr.Length, nHit, nTags);
            }
            //var existingTags = tagFinder.ExtractExistingSequneceTags(sequence);
            Console.Write("{0}/{1}", nHitSpec, nSpec);
        }

        public static string Reverse(string s)
        {
            char[] charArray = s.ToCharArray();
            Array.Reverse(charArray);
            return new string(charArray);
        }
        
        [Test]
        public void TestSequenceTagGlycoData()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            //const string rawFile = @"D:\MassSpecFiles\Glyco\User_sample_test_02252015.raw";
            //const string rawFile = @"D:\MassSpecFiles\CPTAC_Intact_CR33_5_29Jun15_Bane_15-02-01RZ.raw";
            const string rawFile = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";

            if (!File.Exists(rawFile))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFile);
            }

            //var run = PbfLcMsRun.GetLcMsRun(rawFile, rawFile.EndsWith(".mzML") ? MassSpecDataType.MzMLFile : MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
            var run = PbfLcMsRun.GetLcMsRun(rawFile, 0, 0);
            var ms2ScanNums = run.GetScanNumbers(2);
            
            foreach(var scanNum in ms2ScanNums)
            {
                var spectrum = run.GetSpectrum(scanNum) as ProductSpectrum;

                Console.WriteLine(@"ScanNum = {0}; # of Peaks = {1}", scanNum, spectrum.Peaks.Length);
                Console.WriteLine(@"{0}", spectrum.ActivationMethod != ActivationMethod.ETD ? "ETD" : "HCD");
                var tolerance = new Tolerance(5);
                var tagFinder = new SequenceTagFinder(spectrum, tolerance);
                var n = 0;
                foreach (var tag in tagFinder.FindSequenceTags())
                {
                    var seqTags = tag.GetTagStrings();
                    n += seqTags.Count;
                }
                Console.WriteLine(n);
            }
            //var existingTags = tagFinder.ExtractExistingSequneceTags(sequence);
            //Console.Write(scanNum + "\t" + existingTags.Count);
        }    
    }
}
