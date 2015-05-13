using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.PostProcessing;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    class TestFeatureBasedSearch
    {
        [Test]
        public void TestFeatureId()
        {
            const string dataSet = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3";

            // Feature: 5236-5286	6-12	8480.3681	5
            const int minScanNum = 5236;
            const int maxScanNum = 5286;
            const double featureMass = 8480.3681;

            //const int minScanNum = 7251;
            //const int maxScanNum = 7326;
            //const double featureMass = 32347.18;

//            const int minScanNum = 4451;
//            const int maxScanNum = 4541;
//            const double featureMass = 31267.95;

            var tolerance = new Tolerance(10);
            var relaxedTolerance = new Tolerance(20);

            const int minTagLength = 5;
            const int minMergedTagLength = 7;
            const int minNumTagMatches = 1;

            var rawFileName = Path.ChangeExtension(dataSet, ".raw");
            var run = PbfLcMsRun.GetLcMsRun(rawFileName);

            var aminoAcidSet = AminoAcidSet.GetStandardAminoAcidSet();
            var featureFileName = Path.ChangeExtension(dataSet, ".ms1ft");
            var filter = new Ms1FtFilter(run, tolerance, featureFileName);
            var ms2ScanNums =
                filter.GetMatchingMs2ScanNums(featureMass)
                    .Where(scanNum => scanNum > minScanNum && scanNum < maxScanNum)
                    .ToArray();

            const string tagFileName = dataSet + ".seqtag"; //"_MinLength3.seqtag"; //Path.ChangeExtension(dataSet, ".seqtag");
            const string fastaFilePath = @"H:\Research\QCShew_TopDown\Production\ID_002216_235ACCEA.fasta";
            var fastaDb = new FastaDatabase(fastaFilePath);
            var searchableDb = new SearchableDatabase(fastaDb);
            var tagParser = new SequenceTagParser(tagFileName, minTagLength);

            var proteinsToTags = new Dictionary<string, IList<MatchedTag>>();
            foreach (var ms2ScanNum in ms2ScanNums)
            {
                var tags = tagParser.GetSequenceTags(ms2ScanNum);
                foreach (var tag in tags)
                {
                    var matchedIndices = searchableDb.FindAllMatchedSequenceIndices(tag.Sequence).ToArray();
                    foreach (var index in matchedIndices)
                    {
                        var protein = fastaDb.GetProteinName(index);
                        var startIndex = fastaDb.GetZeroBasedPositionInProtein(index);
                        var matchedTag = new MatchedTag(tag, startIndex, featureMass);
                        IList<MatchedTag> existingTags;
                        if (proteinsToTags.TryGetValue(protein, out existingTags))
                        {
                            existingTags.Add(matchedTag);
                        }
                        else
                        {
                            proteinsToTags.Add(protein, new List<MatchedTag> { matchedTag });
                        }
                    }
                }
            }

            foreach (var entry in proteinsToTags.OrderByDescending(e => e.Value.Count))
            {
                if (entry.Value.Count < minNumTagMatches) break;
                var proteinName = entry.Key;
                var proteinSequence = fastaDb.GetProteinSequence(proteinName);
                var protein = new Sequence(proteinSequence, aminoAcidSet);
                Console.WriteLine(proteinName + "\t" + entry.Value.Count);

                var matchedTagSet = new MatchedTagSet(proteinSequence, aminoAcidSet,
                    tolerance, relaxedTolerance);

                Console.WriteLine("********** Before merging");
                foreach (var matchedTag in entry.Value)
                {
                    var seq = proteinSequence.Substring(matchedTag.StartIndex,
                        matchedTag.EndIndex - matchedTag.StartIndex);
                    var nTermMass = protein.GetMass(0, matchedTag.StartIndex);
                    var cTermMass = protein.GetMass(matchedTag.EndIndex, protein.Count);
                    Console.WriteLine("\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
                        (matchedTag.NTermFlankingMass - nTermMass), seq, (matchedTag.CTermFlankingMass - cTermMass), matchedTag.StartIndex,
                        matchedTag.IsNTermFlankingMassReliable, matchedTag.IsCTermFlankingMassReliable);

                    matchedTagSet.Add(matchedTag);
                }

                Console.WriteLine("********** After merging");
                foreach (var matchedTag in matchedTagSet.Tags)
                {
                    if (matchedTag.Length < minMergedTagLength) continue;
                    var seq = proteinSequence.Substring(matchedTag.StartIndex,
                        matchedTag.EndIndex - matchedTag.StartIndex);
                    var nTermMass = protein.GetMass(0, matchedTag.StartIndex);
                    var cTermMass = protein.GetMass(matchedTag.EndIndex, protein.Count);
                    Console.WriteLine("\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
                        (matchedTag.NTermFlankingMass-nTermMass), seq, (matchedTag.CTermFlankingMass - cTermMass), matchedTag.StartIndex,
                        matchedTag.IsNTermFlankingMassReliable, matchedTag.IsCTermFlankingMassReliable);
                }

                break;
            }
        }

        [Test]
        public void TestGetProteinsWithTagMatchingSingleSpec()
        {
            const string dataSet = @"H:\Research\Lewy\raw\Lewy_intact_07";
            //            const int scanNum = 5158;
            const int minTagLength = 7;
            const int minNumTagMatches = 1;
            var aminoAcidSet = AminoAcidSet.GetStandardAminoAcidSet();

            const int scanNum = 2;
            // Parse sequence tags
            const string tagFileName = dataSet + ".seqtag"; //"_MinLength3.seqtag"; //Path.ChangeExtension(dataSet, ".seqtag");
            const string fastaFilePath = @"H:\Research\Lewy\ID_004858_0EE8CF61.fasta";
            var fastaDb = new FastaDatabase(fastaFilePath);
            var searchableDb = new SearchableDatabase(fastaDb);
            var tagParser = new SequenceTagParser(tagFileName, minTagLength);

            var tags = tagParser.GetSequenceTags(scanNum);
            var proteinsToTags = new Dictionary<string, IList<MatchedTag>>();

            foreach (var tag in tags)
            {
                var matchedIndices = searchableDb.FindAllMatchedSequenceIndices(tag.Sequence).ToArray();
                foreach (var index in matchedIndices)
                {
                    var protein = fastaDb.GetProteinName(index);
                    var startIndex = fastaDb.GetOneBasedPositionInProtein(index);
                    var matchedTag = new MatchedTag(tag, startIndex, 0.0);
                    IList<MatchedTag> existingTags;
                    if (proteinsToTags.TryGetValue(protein, out existingTags))
                    {
                        existingTags.Add(matchedTag);
                    }
                    else
                    {
                        proteinsToTags.Add(protein, new List<MatchedTag> { matchedTag });
                    }
                }
            }

            foreach (var entry in proteinsToTags.OrderByDescending(e => e.Value.Count))
            {
                if (entry.Value.Count < minNumTagMatches) break;
                var proteinName = entry.Key;
                var proteinSequence = fastaDb.GetProteinSequence(proteinName);
                var protein = new Sequence(proteinSequence, aminoAcidSet);
                Console.WriteLine(proteinName + "\t" + entry.Value.Count);
                foreach (var matchedTag in entry.Value)
                {
                    var seq = proteinSequence.Substring(matchedTag.StartIndex,
                        matchedTag.EndIndex - matchedTag.StartIndex);
                    var nTermMass = protein.GetMass(0, matchedTag.StartIndex);
                    var cTermMass = protein.GetMass(matchedTag.EndIndex, protein.Count);
                    Console.WriteLine("\t{0} ({1})\t{2}\t{3} ({4})\t{5}\t{6}\t{7}",
                        matchedTag.NTermFlankingMass, (matchedTag.NTermFlankingMass - nTermMass), 
                        seq, 
                        matchedTag.CTermFlankingMass, (matchedTag.CTermFlankingMass - cTermMass), 
                        matchedTag.StartIndex,
                        matchedTag.IsNTermFlankingMassReliable, matchedTag.IsCTermFlankingMassReliable);

                }
            }
        }

        [Test]
        public void TestTagMatchingSingleSpec()
        {
            const string dataSet = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3";
            const int scanNum = 4533;

            // Parse sequence tags
            var tagFileName = Path.ChangeExtension(dataSet, ".seqtag");
            const int minTagLength = 8;
            const string fastaFilePath = @"H:\Research\QCShew_TopDown\Production\ID_002216_235ACCEA.fasta";
            var fastaDb = new FastaDatabase(fastaFilePath);
            var searchableDb = new SearchableDatabase(fastaDb);
            var tagParser = new SequenceTagParser(tagFileName, minTagLength);

            var tags = tagParser.GetSequenceTags(scanNum);
            foreach (var tag in tags)
            {
                var matchedProteins = searchableDb.FindAllMatchedSequenceIndices(tag.Sequence)
                    .Select(index => fastaDb.GetProteinName(index)).ToArray();
                if (matchedProteins.Any())
                {
                    Console.WriteLine("{0}\t{1}\t{2}\t{3}", tag.Sequence, tag.IsPrefix, tag.FlankingMass, string.Join("\t", matchedProteins));
                }
            }
        }

        [Test]
        public void TestFeatureIdMatching()
        {
            const string resultFilePath = @"H:\Research\QCShew_TopDown\Production\M1_V092\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv";
            var resultParser = new MsPathFinderParser(resultFilePath);
            const double qValueThreshold = 0.01;
            const double tolerancePpm = 13;

            const string dataSet = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3";
            var rawFileName = Path.ChangeExtension(dataSet, ".raw");
            var run = PbfLcMsRun.GetLcMsRun(rawFileName);

            var idList =
                resultParser.GetIdList().TakeWhile(id => id.QValue <= qValueThreshold).OrderBy(id => id.Mass).ToList();
            var idMassList = idList.Select(id => id.Mass).ToList();
            var idFlag = new bool[idList.Count];


            // Parse sequence tags
            var tagFileName = Path.ChangeExtension(dataSet, ".seqtag");
            const int minTagLength = 6;
            const int numProtMatches = 4;
//            const string fastaFilePath = @"H:\Research\QCShew_TopDown\Production\ID_002216_235ACCEA.fasta";
            const string fastaFilePath = @"H:\Research\QCShew_TopDown\Production\ID_002216_235ACCEA.icsfldecoy.fasta";
            var fastaDb = new FastaDatabase(fastaFilePath);
            var searchableDb = new SearchableDatabase(fastaDb);

            var tagParser = new SequenceTagParser(tagFileName, minTagLength);
            
            var featureFileName = Path.ChangeExtension(dataSet, ".ms1ft");
            var featureParser = new TsvFileParser(featureFileName);

            var minScan = featureParser.GetData("MinScan").Select(s => Convert.ToInt32(s)).ToArray();
            var maxScan = featureParser.GetData("MaxScan").Select(s => Convert.ToInt32(s)).ToArray();
            var minCharge = featureParser.GetData("MinCharge").Select(s => Convert.ToInt32(s)).ToArray();
            var maxCharge = featureParser.GetData("MaxCharge").Select(s => Convert.ToInt32(s)).ToArray();
            var monoMass = featureParser.GetData("MonoMass").Select(Convert.ToDouble).ToArray();

            var numFeaturesWithId = 0;
            var numFeaturesWithMs2 = 0;
            var numFeaturesWithTags = 0;
            var numFeaturesWithMatchingTags = 0;
            var numFeaturesWithTwoOrMoreMatchingTags = 0;
            var numFeaturesWithNoIdAndMatchingTags = 0;

            for (var i = 0; i < featureParser.NumData; i++)
            {
                var mass = monoMass[i];

                // Find Id
                var tolDa = new Tolerance(tolerancePpm).GetToleranceAsDa(mass, 1);
                var minMass = mass - tolDa;
                var maxMass = mass + tolDa;
                var index = idMassList.BinarySearch(mass);
                if (index < 0) index = ~index;

                var matchedId = new List<MsPathFinderId>();
                // go down
                var curIndex = index - 1;
                while (curIndex >= 0)
                {
                    var curId = idList[curIndex];
                    if (curId.Mass < minMass) break;
                    if (curId.Scan > minScan[i] && curId.Scan < maxScan[i]
                        && curId.Charge >= minCharge[i] && curId.Charge <= maxCharge[i])
                    {
                        matchedId.Add(curId);
                        idFlag[curIndex] = true;
                    }
                    --curIndex;
                }

                // go up
                curIndex = index;
                while (curIndex < idList.Count)
                {
                    var curId = idList[curIndex];
                    if (curId.Mass > maxMass) break;
                    if (curId.Scan >= minScan[i] && curId.Scan <= maxScan[i]
                        && curId.Charge >= minCharge[i] && curId.Charge <= maxCharge[i])
                    {
                        matchedId.Add(curId);
                        idFlag[curIndex] = true;
                    }
                    ++curIndex;
                }

                var hasId = false;
                if (matchedId.Any())
                {
                    ++numFeaturesWithId;
                    hasId = true;
                }

                // Find MS2 scans
//                var numMs2Scans = 0;
                var tags = new List<SequenceTag>();
                var hasMs2 = false;
                for (var scanNum = minScan[i]; scanNum <= maxScan[i]; scanNum++)
                {
                    var isolationWindow = run.GetIsolationWindow(scanNum);
                    if (isolationWindow == null) continue;
                    var isolationWindowTargetMz = isolationWindow.IsolationWindowTargetMz;
                    var charge = (int)Math.Round(mass / isolationWindowTargetMz);
                    if (charge < minCharge[i] || charge > maxCharge[i]) continue;
                    var mz = Ion.GetIsotopeMz(mass, charge,
                        Averagine.GetIsotopomerEnvelope(mass).MostAbundantIsotopeIndex);
                    if (isolationWindow.Contains(mz))
                    {
//                        ++numMs2Scans;
                        tags.AddRange(tagParser.GetSequenceTags(scanNum));
                        hasMs2 = true;
                    }
                }
                if (hasMs2) ++numFeaturesWithMs2;
                if (tags.Any()) ++numFeaturesWithTags;
                var protHist = new Dictionary<string, int>();
                var hasMatchedTag = false;
                foreach (var tag in tags)
                {
                    var matchedProteins = searchableDb.FindAllMatchedSequenceIndices(tag.Sequence).Select(idx => fastaDb.GetProteinName(idx)).ToArray();
                    if (matchedProteins.Any())
                    {
                        hasMatchedTag = true;
                        foreach (var protein in matchedProteins)
                        {
                            int num;
                            if (protHist.TryGetValue(protein, out num)) protHist[protein] = num + 1;
                            else protHist[protein] = 1;
                        }
                    }
                }
                if (hasMatchedTag)
                {
                    ++numFeaturesWithMatchingTags;
                    if (!hasId) ++numFeaturesWithNoIdAndMatchingTags;
                }
                if (protHist.Any())
                {
                    var maxOcc = protHist.Values.Max();
                    if (maxOcc >= numProtMatches) ++numFeaturesWithTwoOrMoreMatchingTags;
                }
            }

            Console.WriteLine("NumFeatures: {0}", featureParser.NumData);
            Console.WriteLine("NumId: {0}", idList.Count);
            Console.WriteLine("NumFeaturesWithId: {0} ({1})", numFeaturesWithId, numFeaturesWithId / (float)featureParser.NumData);
            Console.WriteLine("NumFeaturesWithMs2: {0} ({1})", numFeaturesWithMs2, numFeaturesWithMs2 / (float)featureParser.NumData);
            Console.WriteLine("NumFeaturesWithTag: {0} ({1})", numFeaturesWithTags, numFeaturesWithTags / (float)featureParser.NumData);
            Console.WriteLine("NumFeaturesWithMatchedTag: {0} ({1})", numFeaturesWithMatchingTags, numFeaturesWithMatchingTags / (float)featureParser.NumData);
            Console.WriteLine("NumFeaturesWithMoreThanOneMatchedTag: {0} ({1})", numFeaturesWithTwoOrMoreMatchingTags, numFeaturesWithTwoOrMoreMatchingTags / (float)featureParser.NumData);
            Console.WriteLine("NumFeaturesWithNoIdAndMatchedTag: {0} ({1})", numFeaturesWithNoIdAndMatchingTags, numFeaturesWithNoIdAndMatchingTags / (float)featureParser.NumData);

            for (var i = 0; i < idFlag.Length; i++)
            {
                if (!idFlag[i])
                {
                    Console.WriteLine(idList[i].Scan);
                }
            }


//            Console.WriteLine(string.Join(",", filter.GetMatchingMs2ScanNums(8115.973001)));
//
//            Console.WriteLine(featureFileName);
        }

        [Test]
        public void TestTagMatching()
        {
            // Parse sequence tags
            const string dataSet = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3";
            const int minTagLength = 8;
            var tagFileName = Path.ChangeExtension(dataSet, ".seqtag");
            var tagParser = new SequenceTagParser(tagFileName, minTagLength);

            // Parse raw file
            var rawFileName = Path.ChangeExtension(dataSet, ".raw");
            var run = PbfLcMsRun.GetLcMsRun(rawFileName);

            // Parse ID file
            const string resultFilePath = @"H:\Research\QCShew_TopDown\Production\M1_V092\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv";
            var resultParser = new MsPathFinderParser(resultFilePath);
            const double qValueThreshold = 0.01;
            var idList = resultParser.GetIdWithQValuesNoLargerThan(qValueThreshold);
            var idFlag = new bool[run.MaxLcScan + 1];
            foreach (var id in idList) idFlag[id.Scan] = true;

            const string fastaFilePath = @"H:\Research\QCShew_TopDown\Production\ID_002216_235ACCEA.fasta";
//            const string fastaFilePath = @"H:\Research\QCShew_TopDown\Production\ID_002216_235ACCEA.icsfldecoy.fasta";
            var fastaDb = new FastaDatabase(fastaFilePath);
            var searchableDb = new SearchableDatabase(fastaDb);

            var numMs2Spectra = 0;
            var numSpectraWithTag = 0;
            var numSpectraWithMatchingTag = 0;
            var numSpectraWithMatchedTagNoId = 0;
            foreach (var ms2ScanNum in run.GetScanNumbers(2))
            {
                ++numMs2Spectra;
                var tags = tagParser.GetSequenceTags(ms2ScanNum);
                if (tags != null)
                {
                    ++numSpectraWithTag;
                    foreach (var tag in tags)
                    {
                        if (searchableDb.Search(tag.Sequence) >= 0)
                        {
                            //Console.WriteLine(tag.Sequence);
                            ++numSpectraWithMatchingTag;
                            if (!idFlag[ms2ScanNum]) ++numSpectraWithMatchedTagNoId;
                            break;
                        }
                    }
                }
            }
            Console.WriteLine("Tag length: {0}", minTagLength);
            Console.WriteLine("NumMs2Spectra: {0}", numMs2Spectra);
            Console.WriteLine("NumMs2SpectraWithTags: {0} ({1})", numSpectraWithTag, numSpectraWithTag/(float)numMs2Spectra);
            Console.WriteLine("NumMs2SpectraWithMatchingTags: {0} ({1})", numSpectraWithMatchingTag, numSpectraWithMatchingTag/(float)numMs2Spectra);
            Console.WriteLine("NumMs2SpectraWithMatchingTagsWithNoId: {0} ({1})", numSpectraWithMatchedTagNoId, numSpectraWithMatchedTagNoId / (float)numMs2Spectra);
        }
    }
}
