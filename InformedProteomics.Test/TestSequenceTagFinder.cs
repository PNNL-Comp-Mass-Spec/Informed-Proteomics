using System;
using System.Linq;
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
            //const string TestRawFile = @"D:\\Vlad_TopDown\\raw\\yufeng_column_test2.raw";
            //const string TestResultFile = @"D:\\Vlad_TopDown\\results\\yufeng_column_test2_IcTda.tsv";
            const string TestRawFile = @"D:\MassSpecFiles\training\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            const string TestResultFile = @"D:\MassSpecFiles\results\ProMex\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv";
            
            //const int MaxTags = 100000;
            var tsvParser = new TsvFileParser(TestResultFile);
            var headerList = tsvParser.GetHeaders();
            var tsvData = tsvParser.GetAllData();
            var ms2ScanNumbers = tsvData["Scan"];
        
            var run = PbfLcMsRun.GetLcMsRun(TestRawFile, TestRawFile.EndsWith(".mzML") ? MassSpecDataType.MzMLFile : MassSpecDataType.XCaliburRun, 0, 0);
            var nSpec = 0;
            var nHitSpec = 0;

            for (var i = 0; i < ms2ScanNumbers.Count; i++)
            {
                var scanNum = Int32.Parse(ms2ScanNumbers[i]);
                var spectrum = run.GetSpectrum(scanNum) as ProductSpectrum;

                int tsvIndex = ms2ScanNumbers.FindIndex(x => Int32.Parse(x) == scanNum);
                var qValue = Double.Parse(tsvData["QValue"].ElementAt(tsvIndex));
                if (qValue > 0.01) continue;

                var seqStr = tsvData["Sequence"].ElementAt(tsvIndex).Trim();
                var modStr = tsvData["Modifications"].ElementAt(tsvIndex).Trim();
                //Console.WriteLine(spectrum.ScanNum);
                var tagFinder = new SequenceTagFinder(spectrum, new Tolerance(5), 5);
                var nTags = 0;
                var nHit = 0;
                foreach (var tag in tagFinder.FindSequenceTags())
                {
                    nTags++;
                    //if (tag.Count >= 5 && nTags < MaxTags)
                    //{
                        foreach (var tagStr in tag.GetTagStrings())
                        {
                            if (seqStr.Contains(tagStr) || seqStr.Contains(Reverse(tagStr))) nHit++;
                        }
                    //}
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
            const string rawFile = @"D:\MassSpecFiles\Glyco\User_sample_test_02252015.raw";


            //var run = PbfLcMsRun.GetLcMsRun(rawFile, rawFile.EndsWith(".mzML") ? MassSpecDataType.MzMLFile : MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
            var run = PbfLcMsRun.GetLcMsRun(rawFile, rawFile.EndsWith(".mzML") ? MassSpecDataType.MzMLFile : MassSpecDataType.XCaliburRun, 0, 0);
            var ms2ScanNums = run.GetScanNumbers(2);
            
            //for (var i = 12700; i < 13000; i++)
            //{
                //var scanNum = ms2ScanNums[i];
                var scanNum = 14075;

                //if (run.GetMsLevel(scanNum) != 2) continue;

                var spectrum = run.GetSpectrum(scanNum) as ProductSpectrum;

                Console.WriteLine(@"ScanNum = {0}; # of Peaks = {1}", scanNum, spectrum.Peaks.Length);
                Console.WriteLine(@"{0}", spectrum.ActivationMethod != ActivationMethod.ETD ? "ETD" : "HCD"); 

                var tagFinder = new SequenceTagFinder(spectrum, new Tolerance(5));
                var n = 0;
                foreach (var tag in tagFinder.FindSequenceTags())
                {
                    var seqTags = tag.GetTagStrings();
                    n += seqTags.Length;
                }
                Console.WriteLine(n);
            //}
            //var existingTags = tagFinder.ExtractExistingSequneceTags(sequence);
            //Console.Write(scanNum + "\t" + existingTags.Count);
        }    
    }
}
