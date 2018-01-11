using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.SequenceTag;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    class TestLongestSequence
    {
        [Test]
        public void GetLongestSequenceTest()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var rawFiles = new []
            {/*@"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_3_11Feb15_Bane_C2Column5.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep3_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep4_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep5_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep7_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep8_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_01_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_02_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_03_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_FA_01_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_FA_02_15Jan15_Bane_C2-14-08-02RZ.pbf",
             @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_FA_03_15Jan15_Bane_C2-14-08-02RZ.pbf",*/
             @"\\protoapps\UserData\Jungkap\Joshua\testData\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf"
            };
            var resultFiles = new string[]
            {/*@"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_BE100_PO4_3_11Feb15_Bane_C2Column5_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_rep3_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_rep4_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_rep5_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_rep7_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_rep8_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_SDS_T_01_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_SDS_T_02_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_SDS_T_03_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_SDS_T_FA_01_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_SDS_T_FA_02_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\CPTAC_Intact_SDS_T_FA_03_15Jan15_Bane_C2-14-08-02RZ_IcTda.tsv",*/
             @"\\protoapps\UserData\Jungkap\Joshua\IdResult\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv"
            };

            for (int i = 0; i < 1; i++)
            {
                var stringPieces = rawFiles[i].Split(new char[] {'\\', '.'});
                var fileName = stringPieces[stringPieces.Length - 2];
                Console.WriteLine(fileName);
                var csv = new StringBuilder();
                csv.Append("scanNumber,#MatchedFragments,OriginalSequence,PrefixStartIndex,PrefixEndIndex,SuffixStartIndex,SuffixEndIndex,\n");
                string outputFile = (@"\\protoapps\UserData\Jungkap\Joshua\SequenceOutFiles\longestTagOutput_"+ fileName + ".csv");
                var resultLine = ProcessFile(rawFiles[i], resultFiles[i], methodName);
                csv.Append(resultLine);
                File.WriteAllText(outputFile, csv.ToString());
            }
        }

        public string ProcessFile(string rawFile, string resultFile, string methodName)
        {
            if (!File.Exists(rawFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFile);
                return "\n";
            }

            if (!File.Exists(resultFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, resultFile);
                return "\n";
            }

            var tsvParser = new TsvFileParser(resultFile);
            var headerList = tsvParser.GetHeaders();
            var tsvData = tsvParser.GetAllData();
            var ms2ScanNumbers = tsvData["Scan"];

            var run = PbfLcMsRun.GetLcMsRun(rawFile, 0, 0);

            var resultLine = "";
            for (int i = 0; i < ms2ScanNumbers.Count; i++)
            {
                var scanNum = Int32.Parse(ms2ScanNumbers[i]);
                var spectrum = run.GetSpectrum(scanNum) as ProductSpectrum;
                int tsvIndex = ms2ScanNumbers.FindIndex(x => Int32.Parse(x) == scanNum);

                var qValue = Double.Parse(tsvData["QValue"].ElementAt(tsvIndex));
                if (qValue > 0.01) continue;

                var seqStr = tsvData["Sequence"].ElementAt(tsvIndex).Trim();
                var seqMod = tsvData["Modifications"].ElementAt(tsvIndex).Trim();
                var matchedFrags = tsvData["#MatchedFragments"].ElementAt(tsvIndex).Trim();
                var aaSet = new AminoAcidSet();
                var sequence = Sequence.CreateSequence(seqStr, seqMod, aaSet);
                var tol = new Tolerance(10);
                var sequenceFinder = new SequenceTagIndexFinder(tol, 1, 10);
                var results = sequenceFinder.GetLongestSequence(spectrum, sequence);
                resultLine += String.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},\n", scanNum, matchedFrags, seqStr, results.Item1, results.Item2,results.Item3,results.Item4,results.Item5,results.Item6);
            }
            return resultLine;
        }
    }
}
