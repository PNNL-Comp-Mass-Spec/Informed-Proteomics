using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.PostProcessing;
using InformedProteomics.TopDown.Scoring;
using InformedProteomics.TopDown.TagBasedSearch;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestSequenceTagMatching
    {


        [Test]
        public void TestSearchWithTagGeneration()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rawFilePath = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";

            if (!File.Exists(rawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFilePath);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);
            const string fastaFilePath = @"D:\MSPathFinder\Fasta\ID_002216_235ACCEA.fasta";
            //const string fastaFilePath = @"D:\MassSpecFiles\60k\ID_004973_9BA6912F.fasta";
            if (!File.Exists(fastaFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, fastaFilePath);
                return;
            }

            var fastaDb = new FastaDatabase(fastaFilePath);
            var tolerance = new Tolerance(10);
            var modsFilePath = @"D:\MSPathFinder\Fasta\Mods.txt";

            if (!File.Exists(modsFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, modsFilePath);
                return;
            }

            var aaSet = new AminoAcidSet(modsFilePath);

            //TestTagBasedSearch(run, fastaDb, tolerance, aaSet);
            var tagSearchEngine = new ScanBasedTagSearchEngine(run, new SequenceTagGenerator(run, new Tolerance(8)),  new LcMsPeakMatrix(run),  fastaDb, tolerance,aaSet);

            var matchedTags = tagSearchEngine.RunSearch(4672);
            foreach (var match in matchedTags)
            {
                Console.Write(match.Sequence);
                Console.WriteLine("\t{0}\t{1}\t{2}", match.TagMatch.StartIndex, match.TagMatch.EndIndex, match.TagMatch.Mass);
            }


        }
        
        
        [Test]
        public void TestTagBasedSearchForLewy()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string rawFilePath = @"D:\MassSpecFiles\Lewy\Lewy_AT_AD1_21May15_Bane_14-09-01RZ.pbf";
            if (!File.Exists(rawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFilePath);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);

            const int minTagLength = 4;
            var tagFilePath = MassSpecDataReaderFactory.ChangeExtension(rawFilePath, ".seqtag");
            //var tagParser = new SequenceTagParser(tagFilePath, minTagLength, 10000);

            const string fastaFilePath = @"D:\MassSpecFiles\Lewy\a4_human.fasta";
            if (!File.Exists(fastaFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, fastaFilePath);
                return;
            }

            var fastaDb = new FastaDatabase(fastaFilePath);

            var tolerance = new Tolerance(10);
            var modsFilePath = @"D:\MassSpecFiles\Lewy\Mods.txt";
            if (!File.Exists(modsFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, modsFilePath);
                return;
            }

            var aaSet = new AminoAcidSet(modsFilePath);            

            TestTagBasedSearch(run, fastaDb, tolerance, aaSet);
        }

        [Test]
        public void TestTagBasedSearchCompRef()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string dataSetPath = @"D:\MassSpecFiles\CompRef";
            const string fastaFilePath = @"D:\MassSpecFiles\CompRef\ID_003278_4B4B3CB1.fasta";
            const string modsFilePath = @"D:\MassSpecFiles\CompRef\Mods.txt";

            if (!Directory.Exists(dataSetPath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, dataSetPath);
                return;
            }
            if (!File.Exists(modsFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, modsFilePath);
                return;
            }
            if (!File.Exists(fastaFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, fastaFilePath);
                return;
            }

            var fileEntries = Directory.GetFiles(dataSetPath);

            var dataset = (from fileName in fileEntries where fileName.EndsWith("pbf") select Path.GetFileNameWithoutExtension(fileName)).ToList();
            dataset.Sort();
          
            
            var fastaDb = new FastaDatabase(fastaFilePath);
            var tolerance = new Tolerance(10);
            var aaSet = new AminoAcidSet(modsFilePath);

            for (var i = 0; i < dataset.Count; i++)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", dataSetPath, dataset[i]);
                var ms1File = string.Format(@"{0}\{1}.ms1ft", dataSetPath, dataset[i]);
                var tagFilePath = MassSpecDataReaderFactory.ChangeExtension(rawFile, ".seqtag");

                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                const int minTagLength = 5;
                //var tagParser = new SequenceTagParser(tagFilePath, minTagLength, 100);

                Console.WriteLine("-----------------{0}--------------------", rawFile);

                TestTagBasedSearch(run, fastaDb, tolerance, aaSet);

                Console.WriteLine("-----------------------------------------------------------------------");
            }

        }

        
        
        [Test]
        public void TestTagBasedSearch()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

//            const string rawFilePath = @"H:\Research\Lewy\raw\Lewy_intact_01.raw";
//            const string rawFilePath = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
//            const string rawFilePath = @"H:\Research\Yufeng\TopDownYufeng\raw\yufeng_column_test2.raw";
//            const string rawFilePath = @"H:\Research\Weijun_TopDown\raw\UC4_Intact_plasmaTest_90_6May15_Bane_14-09-01RZ.raw";
//            const string rawFilePath = @"H:\Research\Charles\TopDown\raw\SBEP_STM_001_02272012_Aragon.raw";
            const string rawFilePath = @"D:\MassSpecFiles\60k\Yufeng_SampleTest1_150614113438.pbf";
            //const string rawFilePath = @"D:\MassSpecFiles\60k\NCR_50K_Test_24Jun15_Bane_15-02-02RZ.pbf";

            if (!File.Exists(rawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFilePath);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);

            const int minTagLength = 5;
            var tagFilePath = MassSpecDataReaderFactory.ChangeExtension(rawFilePath, ".seqtag");
            //var tagParser = new SequenceTagParser(tagFilePath, minTagLength, 100);

            const string fastaFilePath = @"D:\MassSpecFiles\60k\ID_003836_DA9CC1E4.fasta";
            //const string fastaFilePath = @"D:\MassSpecFiles\60k\ID_004973_9BA6912F.fasta";
            if (!File.Exists(fastaFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, fastaFilePath);
                return;
            }

            var fastaDb = new FastaDatabase(fastaFilePath);

            var tolerance = new Tolerance(10);
//            var modsFilePath = @"H:\Research\QCShew_TopDown\Production\Mods_Methyl.txt";
            var modsFilePath = @"D:\MassSpecFiles\60k\Mods.txt";
//            var modsFilePath = "";

            if (!File.Exists(modsFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, modsFilePath);
                return;
            }

            var aaSet = new AminoAcidSet(modsFilePath);

            TestTagBasedSearch(run, fastaDb, tolerance, aaSet);
        }

        private void TestTagBasedSearch(LcMsRun run, 
            FastaDatabase fastaDb, Tolerance tolerance, AminoAcidSet aaSet)
        {
            var engine = new ScanBasedTagSearchEngine(run, new SequenceTagGenerator(run, new Tolerance(8)), new LcMsPeakMatrix(run), fastaDb, tolerance, aaSet);
//            engine.MinScan = 3400;
//            engine.MaxScan = 3900;
            engine.RunSearch();
        }

        [Test]
        public void FindProteins()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string fastaFilePath = @"H:\Research\QCShew_TopDown\Production\ID_002216_235ACCEA.fasta";
            if (!File.Exists(fastaFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, fastaFilePath);
                return;
            }

            var fastaDb = new FastaDatabase(fastaFilePath);
            var searchableDb = new SearchableDatabase(fastaDb);

            const string tagFilePath = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3_seqtag.tsv";
            if (!File.Exists(tagFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, tagFilePath);
                return;
            }

            const string outputFilePath = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3_matchedtag.tsv";
            using (var writer = new StreamWriter(outputFilePath))
            {
                var isHeader = true;
                foreach (var line in File.ReadAllLines(tagFilePath))
                {
                    if (isHeader)
                    {
                        isHeader = false;
                        writer.WriteLine(line+"\t"+"Proteins");
                        continue;
                    }

                    var token = line.Split('\t');
                    if (token.Length != 3) continue;
                    var tag = token[1];

                    var matchedProteins =
                        searchableDb.FindAllMatchedSequenceIndices(tag)
                            .Select(index => fastaDb.GetProteinName(index))
                            .Distinct().ToArray();
                    var matchedProteinStr = string.Join(",", matchedProteins);
                    writer.WriteLine("{0}\t{1}\t{2}", line, matchedProteins.Length, matchedProteinStr);
                }
            }
            Console.WriteLine(@"Done");
        }

        [Test]
        public void CountMatchedProteins()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const int minTagLength = 3;

            var scanToProtein = new Dictionary<int, string>();
            var idTag = new Dictionary<int, bool>();
            const string resultFilePath = @"H:\Research\ProMex\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv";
            if (!File.Exists(resultFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, resultFilePath);
                return;
            }

            var parser = new TsvFileParser(resultFilePath);
            var scans = parser.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
            var proteinNames = parser.GetData("ProteinName").ToArray();
            var qValues = parser.GetData("QValue").Select(Convert.ToDouble).ToArray();
            for (var i = 0; i < qValues.Length; i++)
            {
                if (qValues[i] > 0.01) break;
                scanToProtein.Add(scans[i], proteinNames[i]);
                idTag.Add(scans[i], false);
            }

            const string rawFilePath = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            if (!File.Exists(rawFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFilePath);
                return;
            }

            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);

            const string fastaFilePath = @"H:\Research\QCShew_TopDown\Production\ID_002216_235ACCEA.fasta";
            if (!File.Exists(fastaFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, fastaFilePath);
                return;
            }

//            const string fastaFilePath = @"H:\Research\QCShew_TopDown\Production\ID_002216_235ACCEA.icsfldecoy.fasta";
//            const string fastaFilePath =
//                @"D:\Research\Data\CommonContaminants\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            var fastaDb = new FastaDatabase(fastaFilePath);
            var searchableDb = new SearchableDatabase(fastaDb);
            Console.WriteLine("Sequence length: {0}", fastaDb.GetSequence().Length);

            const string tagFilePath = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3.seqtag";
            if (!File.Exists(tagFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, tagFilePath);
                return;
            }

            var hist = new Dictionary<int, int>();
            
            var scanSet = new HashSet<int>();
            HashSet<string> proteinSetForThisScan = null;
            var prevScan = -1;
            var totalNumMatches = 0L;
            var isHeader = true;
            foreach (var line in File.ReadAllLines(tagFilePath))
            {
                if (isHeader)
                {
                    isHeader = false;
                    continue;
                }

                var token = line.Split('\t');
                if (token.Length < 3) continue;
                var scan = Convert.ToInt32(token[0]);
                var proteinId = scanToProtein.ContainsKey(scan) ? scanToProtein[scan] : null;

                if (scan != prevScan)
                {
                    if (proteinSetForThisScan != null)
                    {
                        var numMatches = proteinSetForThisScan.Count;
                        int numOcc;
                        if (hist.TryGetValue(numMatches, out numOcc)) hist[numMatches] = numOcc + 1;
                        else hist.Add(numMatches, 1);
                    }

                    prevScan = scan;
                    proteinSetForThisScan = new HashSet<string>();
                }

                scanSet.Add(scan);
                var tag = token[1];
                if (tag.Length < minTagLength) continue;

                if (proteinSetForThisScan == null) continue;

                var numMatchesForThisTag = 0;
                foreach (var matchedProtein in searchableDb.FindAllMatchedSequenceIndices(tag)
                    .Select(index => fastaDb.GetProteinName(index)))
                {
                    proteinSetForThisScan.Add(matchedProtein);
                    ++numMatchesForThisTag;

                    if (proteinId != null && matchedProtein.Equals(proteinId))
                    {
                        idTag[scan] = true;
                    }
                }
                totalNumMatches += numMatchesForThisTag;
//                if (numMatchesForThisTag > 10)
//                {
//                    Console.WriteLine("{0}\t{1}", tag, numMatchesForThisTag);
//                }
            }

            if (proteinSetForThisScan != null)
            {
                var numMatches = proteinSetForThisScan.Count;
                int numOcc;
                if (hist.TryGetValue(numMatches, out numOcc)) hist[numMatches] = numOcc + 1;
                else hist.Add(numMatches, 1);
            }
            
            Console.WriteLine("AvgNumMatches: {0}", totalNumMatches/(float)scanSet.Count);
            Console.WriteLine("Histogram:");
            foreach (var entry in hist.OrderBy(e => e.Key))
            {
                Console.WriteLine("{0}\t{1}", entry.Key, entry.Value);
            }

            Console.WriteLine("NumId: {0}", idTag.Count);
            Console.WriteLine("NumIdByTag: {0}", idTag.Select(e => e.Value).Count(v => v));
        }

        [Test]
        public void CountMatchedScansPerProtein()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const int minTagLength = 6;

            var proteinToScan = new Dictionary<string, HashSet<int>>();
            const string fastaFilePath = @"D:\Research\Data\CommonContaminants\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            if (!File.Exists(fastaFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, fastaFilePath);
                return;
            }

            var fastaDb = new FastaDatabase(fastaFilePath);
            var searchableDb = new SearchableDatabase(fastaDb);
            Console.WriteLine(@"Sequence length: {0}", fastaDb.GetSequence().Length);

            //const string tagFilePath = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3_seqtag.tsv";
            //const string tagFilePath = @"\\protoapps\UserData\Jungkap\Co_culture\23B_pellet_TD_3Feb14_Bane_PL011402.seqtag";
            const string tagFilePath = @"D:\MassSpecFiles\co_culture\23A_pellet_TD_3Feb14_Bane_PL011402.seqtag";
            if (!File.Exists(tagFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, tagFilePath);
                return;
            }

            var isHeader = true;
            var numMatchedPairs = 0;
            foreach (var line in File.ReadAllLines(tagFilePath))
            {
                if (isHeader)
                {
                    isHeader = false;
                    continue;
                }

                var token = line.Split('\t');
                if (token.Length != 3) continue;
                var scan = Convert.ToInt32(token[0]);

                var tag = token[1];
                if (tag.Length < minTagLength) continue;

                foreach (var matchedProtein in searchableDb.FindAllMatchedSequenceIndices(tag)
                    .Select(index => fastaDb.GetProteinName(index)))
                {
                    ++numMatchedPairs;
                    HashSet<int> matchedScans;
                    if (proteinToScan.TryGetValue(matchedProtein, out matchedScans))
                    {
                        matchedScans.Add(scan);
                    }
                    else
                    {
                        matchedScans = new HashSet<int> {scan};
                        proteinToScan.Add(matchedProtein, matchedScans);
                    }
                }
            }

            var numMatchedProteins = proteinToScan.Keys.Count;
            var numAllProteins = fastaDb.GetNumEntries();
            Console.WriteLine("NumAllProteins: {0}", numAllProteins);
            Console.WriteLine("NumMatchedProteins: {0}", numMatchedProteins);
            Console.WriteLine("AvgMatchedScansPerProtein: {0}", numMatchedPairs / (float)numAllProteins);
        }

        [Test]
        public void FindProteinDeltaMass()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string folderPath = @"D:\MassSpecFiles\Glyco\";
            if (!Directory.Exists(folderPath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since folder not found: {1}", methodName, folderPath);
                return;
            }

            var fileSet = new string[]
            {
                "User_sample_test_02252015", "User_sample_test_MWCO_02262016", "User_sample_test_SEC_F3_03022105",
                "User_sample_test_SEC_F1_02272015", "User_sample_test_SEC_F2_02282015"
            };
            const string fastaFilePath = folderPath + "ID_003836_DA9CC1E4.fasta";

            for (var i = 0; i < fileSet.Length; i++)
            {
                var datasetName = fileSet[i];
                var tagFilePath = folderPath + datasetName + ".seqtag";
                //var outputFilePath = folderPath + datasetName + ".matchedtag";
                var outputFilePath = folderPath + datasetName + ".dmass";
                var fastaDb = new FastaDatabase(fastaFilePath);
                var searchableDb = new SearchableDatabase(fastaDb);

                using (var writer = new StreamWriter(outputFilePath))
                {
                    var isHeader = true;
                    var nReadSeqTag = 0;

                    Console.WriteLine(@"Reading {0} file", tagFilePath);

                    var nColumn = 0;
                    foreach (var line in File.ReadAllLines(tagFilePath))
                    {
                        if (isHeader)
                        {
                            isHeader = false;
                            nColumn = line.Split('\t').Length;
                            writer.WriteLine(line + "\t" + "Protein" + "\t" + "DetectedFlankingMass" + "\t" + "ExpectedFlankingMass" + "\t" + "DeltaMass");
                            continue;
                        }

                        var token = line.Split('\t');
                        if (token.Length != nColumn) continue;
                        var tag = token[1];
                        //var scan = Convert.ToInt32(token[0]);

                        if (tag.Length < 6) continue;

                        var nTerminal = token[2].Equals("1");
                        var detectedFlankingMass = Double.Parse(token[3]);

                        if (!nTerminal) detectedFlankingMass -= Composition.H2O.Mass;

                        nReadSeqTag++;

                        var matchedProteins =
                            searchableDb.FindAllMatchedSequenceIndices(tag)
                                .Select(index => fastaDb.GetProteinName(index))
                                .Distinct().ToArray();

                        if (matchedProteins.Length < 1) continue;

                        foreach (var protName in matchedProteins)
                        {
                            var seqStr = fastaDb.GetProteinSequence(protName);
                            var oriSeq = new Sequence(seqStr, AminoAcidSet.GetStandardAminoAcidSet());
                            
                            var startIdx = 0;
                            while (true)
                            {
                                var idx = seqStr.IndexOf(tag, startIdx);

                                if (idx < 0) break; //no matching

                                //var nClv = (nTerminal) ? idx : seqStr.Length - idx - tag.Length;
                                var nClv = (nTerminal) ? 2 : 1; 

                                for (var j = 0; j < nClv; j++)
                                {
                                    var flankComposition = (nTerminal)
                                        ? oriSeq.GetComposition(j, idx)
                                        : oriSeq.GetComposition(idx + tag.Length, oriSeq.Count - j);

                                    var massDiff = (detectedFlankingMass - flankComposition.Mass);
                                    if (massDiff > -500 && massDiff < 2000)
                                    {
                                        //writer.WriteLine(massDiff);
                                        writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", line, protName, detectedFlankingMass, flankComposition.Mass, massDiff);
                                    }

                                    if (massDiff > 2000) break;
                                }
                                
                                startIdx = idx + tag.Length;
                            }
                        }
                        //var matchedProteinStr = string.Join(",", matchedProteins);
                        //var massDiffStr = string.Join(",", massDiffList);
                        //writer.WriteLine("{0}\t{1}\t{2}\t{3}", line, matchedProteins.Length, matchedProteinStr, massDiffStr);
                    }

                    Console.WriteLine(@"{0} seq tags are processed", nReadSeqTag);
                }
                Console.WriteLine(@"Done");
            }
        }

    }
}
