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
using InformedProteomics.TopDown.Scoring;
using MathNet.Numerics;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestIcTopDown
    {
        [Test]
        public void TestSbepSearch()
        {
            const string dirPath = @"C:\cygwin\home\kims336\Data\TopDown\raw";
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDown\databases\ID_002166_F86E3B2F.fasta";
            foreach (var specFilePath in Directory.GetFiles(dirPath))
            {
                if (!specFilePath.EndsWith(".raw")) continue;
                Console.WriteLine("Processing {0}", specFilePath);
                TestSbepSearch(specFilePath, dbFilePath);
            }
        }

        [Test]
        public void TestSbepSearch(string specFilePath, string dbFilePath)
        {
            // Search parameters
            const int maxNumNTermCleavages = 0;  // 30
            const int maxNumCTermCleavages = 0;
            const int minLength = 30;    // 7
            const int maxLength = 250; // 1000
            const int minPrecursorIonCharge = 3; // 3
            const int maxPrecursorIonCharge = 30;// 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 10;// 10
            const int numMaxModsPerProtein = 0; // 6

            var precursorTolerance = new Tolerance(10);
            var productIonTolerance = new Tolerance(10);

            // Configure amino acid set
            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var cysteinylC = new SearchModification(Modification.Cysteinyl, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                //pyroGluQ,
                //dehydroC,
                //cysteinylC,
                //deamdN,
                //deamdQ,
                //glutathioneC,
                //oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            TestTopDownSearchInternalSequences(dbFilePath, specFilePath, aaSet, minLength, maxLength, 
                minPrecursorIonCharge, maxPrecursorIonCharge,
                minProductIonCharge, maxProductIonCharge, precursorTolerance, productIonTolerance, false, false);
            //TestTopDownSearch(dbFilePath, specFilePath, aaSet, minLength, maxLength, maxNumNTermCleavages, maxNumCTermCleavages,
            //    minPrecursorIonCharge, maxPrecursorIonCharge,
            //    minProductIonCharge, maxProductIonCharge, precursorTolerance, productIonTolerance, false, false);
            //TestTopDownSearch(dbFilePath, specFilePath, aaSet, minLength, maxLength, maxNumNTermCleavages, maxNumCTermCleavages,
            //    minPrecursorIonCharge, maxPrecursorIonCharge,
            //    minProductIonCharge, maxProductIonCharge, precursorTolerance, productIonTolerance, false, true);
        }

        [Test]
        public void TestTopDownSearchInternalSequences(
            string dbFilePath, string specFilePath, AminoAcidSet aaSet,
            int minLength, int maxLength,
            int minPrecursorIonCharge, int maxPrecursorIonCharge,
            int minProductIonCharge, int maxProductIonCharge,
            Tolerance precursorTolerance, Tolerance productIonTolerance,
            bool ultraMod,
            bool isDecoy
            )
        {

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            Console.Write("Reading raw file...");
            //var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 0);
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

            //var scoringModel = new LikelihoodScoringModel(@"C:\cygwin\home\kims336\Data\TopDown\raw\CorrScores_SBEP.txt");
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);


            var targetDb = new FastaDatabase(dbFilePath);
            targetDb.Read();

            var db = !isDecoy ? targetDb : targetDb.Decoy(null, true);  // shuffled decoy

            var indexedDb = new IndexedDatabase(db);

            //var annotationsAndOffsets = indexedDb.IntactSequenceAnnotationsAndOffsets(minLength, maxLength, maxNumCTermCleavages);
            var annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(minLength, maxLength);

            var numProteins = 0;
            long totalProtCompositions = 0;
            long numPrecursorIons = 0;
            long numPrecursorIonsPassingFilter = 0;

            sw.Reset();
            sw.Start();

            var bestScorePerScan = new Dictionary<int, double>();
            var bestResultPerScan = new Dictionary<int, string>();

            var lcMsCache = new CachedLcMsRun(run, 
                minPrecursorIonCharge, maxPrecursorIonCharge, 
                minProductIonCharge, maxProductIonCharge,
                600.0, 1800.0, precursorTolerance, productIonTolerance);

            foreach (var annotationAndOffset in annotationsAndOffsets)
            {
                ++numProteins;

                var annotation = annotationAndOffset.Annotation;
                var offset = annotationAndOffset.Offset;

                //                    Console.WriteLine(annotation);
                if (numProteins % 1000 == 0)
                {
                    Console.Write("Processing {0}{1} proteins...", numProteins,
                        numProteins == 1 ? "st" : numProteins == 2 ? "nd" : numProteins == 3 ? "rd" : "th");
                    if (numProteins != 0)
                    {
                        sw.Stop();
                        sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
                        Console.WriteLine("Elapsed Time: {0:f4} sec", sec);
                        sw.Reset();
                        sw.Start();
                    }
                    //if (numProteins == 10) break;
                }

                //Console.WriteLine(protAnnotation);

                var seqGraph = SequenceGraph.CreateGraph(aaSet, annotation);
                if (seqGraph == null)
                {
                    Console.WriteLine("Ignoring illegal protein: {0}", annotation);
                    continue;
                }

                //var compSet = new HashSet<Composition>();
                var protCompositions = seqGraph.GetSequenceCompositionsWithNTermCleavage(0);
                if (ultraMod)
                    Console.WriteLine("#NTermCleavages: {0}, #ProteinCompositions: ", 0);
                for (var modIndex = 0; modIndex < protCompositions.Length; modIndex++)
                {
                    if (ultraMod)
                    {
                        if (modIndex % 100 == 0) Console.WriteLine("ModIndex: " + modIndex);
                        //                                if (modIndex >= 100) break;
                    }

                    seqGraph.SetSink(modIndex, 0);
                    var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();
                    totalProtCompositions++;
                    var matches = lcMsCache.GetMs2Matches(protCompositionWithH2O);
                    var modCombinations = seqGraph.ModificationParams.GetModificationCombination(modIndex);

                    foreach (var match in matches)
                    {
                        var charge = match.PrecursorCharge;
                        var ms2ScanNum = match.ScanNum;
                        var spec = run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                        if (spec == null) continue;
                        var scorer = lcMsCache.GetMs2Scorer(ms2ScanNum);

                        //var scorer = new LikelihoodScorer(scoringModel, spec, productIonTolerance, minProductIonCharge, maxProductIonCharge);
                        var score = seqGraph.GetScore(charge, scorer);

                        if (score <= -50) continue;

                        double existingBestScore;
                        if (bestScorePerScan.TryGetValue(ms2ScanNum, out existingBestScore) &&
                            score <= existingBestScore) continue;

                        // new best score
                        var precursorIon = new Ion(protCompositionWithH2O, charge);
                        var sequence = annotation.Substring(2,
                            annotation.Length - 4);
                        var proteinName = targetDb.GetProteinName(offset);
                        var start = targetDb.GetZeroBasedPositionInProtein(offset) + 1;
                        var end = start + sequence.Length - 1;
                        var protLength = targetDb.GetProteinLength(proteinName);
                        bestScorePerScan[ms2ScanNum] = score;
                        bestResultPerScan[ms2ScanNum] = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}",
                            annotation.Substring(2, annotation.Length - 4),   // Sequence
                            modCombinations,    // Modifications
                            protCompositionWithH2O,   // Composition
                            (isDecoy ? FastaDatabase.DecoyProteinPrefix + "_" : "") + proteinName,  // ProteinName
                            targetDb.GetProteinDescription(offset), // ProteinDescription
                            protLength, // ProteinLength
                            start,   // Start
                            end,    // End
                            charge, // precursorCharge
                            precursorIon.GetMostAbundantIsotopeMz(),    // MostAbundantIsotopeMz
                            protCompositionWithH2O.Mass,
                            score);
                    }

                    numPrecursorIonsPassingFilter++;
                }
            }

            // write results into a file
            var icExtension = !isDecoy ? ".ic2result" : "decoy.ic2result";
            var outputFilePath = Path.ChangeExtension(specFilePath, icExtension);
            using (var writer = new StreamWriter(outputFilePath))
            {
                //                writer.WriteLine("ScanNum\tAnnotation\tProtein\tProteinDesc\tComposition\tCharge\tBaseIsotopeMz\tMass\tScore");
                writer.WriteLine("ScanNum\tSequence\tModifications\tComposition\tProteinName\tProteinDesc\tProteinLength\tStart\tEnd\tCharge\tMostAbundantIsotopeMz\tScore");
                var ms2Scans = new List<int>(bestScorePerScan.Keys);
                ms2Scans.Sort();
                foreach (var ms2ScanNum in bestScorePerScan.OrderByDescending(e => e.Value).Select(scanScorePair => scanScorePair.Key))
                {
                    writer.WriteLine(ms2ScanNum + "\t" + bestResultPerScan[ms2ScanNum]);
                }
            }

            sw.Stop();
            Console.WriteLine("NumProteins: {0}", numProteins);
            Console.WriteLine("NumProteinCompositions: {0}", totalProtCompositions);
            Console.WriteLine("NumPrecursorIons: {0}", numPrecursorIons);
            Console.WriteLine("NumPrecursorIonsWithEvidence: {0}", numPrecursorIonsPassingFilter);

            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void TestTopDownSearch(
            string dbFilePath, string specFilePath, AminoAcidSet aaSet,
            int minLength, int maxLength, 
            int maxNumNTermCleavages, int maxNumCTermCleavages,
            int minPrecursorIonCharge, int maxPrecursorIonCharge,
            int minProductIonCharge, int maxProductIonCharge,
            Tolerance precursorTolerance, Tolerance productIonTolerance,
            bool ultraMod,
            bool isDecoy
            )
        {

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            Console.Write("Reading raw file...");
            //var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 0);
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

            //var scoringModel = new LikelihoodScoringModel(@"C:\cygwin\home\kims336\Data\TopDown\raw\CorrScores_SBEP.txt");
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);


            var targetDb = new FastaDatabase(dbFilePath);
            targetDb.Read();

            var db = !isDecoy ? targetDb : targetDb.Decoy(null, true);  // shuffled decoy

            var indexedDb = new IndexedDatabase(db);

            //var annotationsAndOffsets = indexedDb.IntactSequenceAnnotationsAndOffsets(minLength, maxLength, maxNumCTermCleavages);
            var annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(minLength, maxLength);

            var numProteins = 0;
            long totalProtCompositions = 0;
            long numPrecursorIons = 0;
            long numPrecursorIonsPassingFilter = 0;

            sw.Reset();
            sw.Start();

            var bestScorePerScan = new Dictionary<int, double>();
            var bestResultPerScan = new Dictionary<int, string>();

            foreach (var annotationAndOffset in annotationsAndOffsets)
            {
                ++numProteins;

                var annotation = annotationAndOffset.Annotation;
                var offset = annotationAndOffset.Offset;

                //                    Console.WriteLine(annotation);
                if (numProteins % 1000 == 0)
                {
                    Console.Write("Processing {0}{1} proteins...", numProteins,
                        numProteins == 1 ? "st" : numProteins == 2 ? "nd" : numProteins == 3 ? "rd" : "th");
                    if (numProteins != 0)
                    {
                        sw.Stop();
                        sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
                        Console.WriteLine("Elapsed Time: {0:f4} sec", sec);
                        sw.Reset();
                        sw.Start();
                    }
                    //if (numProteins == 10) break;
                }

                //Console.WriteLine(protAnnotation);

                var seqGraph = SequenceGraph.CreateGraph(aaSet, annotation);
                if (seqGraph == null)
                {
                    Console.WriteLine("Ignoring illegal protein: {0}", annotation);
                    continue;
                }

                for (var numNTermCleavage = 0; numNTermCleavage <= maxNumNTermCleavages; numNTermCleavage++)
                {
                    //var compSet = new HashSet<Composition>();
                    var protCompositions = seqGraph.GetSequenceCompositionsWithNTermCleavage(numNTermCleavage);
                    if (ultraMod)
                        Console.WriteLine("#NTermCleavages: {0}, #ProteinCompositions: ", numNTermCleavage);
                    for (var modIndex = 0; modIndex < protCompositions.Length; modIndex++)
                    {
                        if (ultraMod)
                        {
                            if (modIndex % 100 == 0) Console.WriteLine("ModIndex: " + modIndex);
                            //                                if (modIndex >= 100) break;
                        }

                        seqGraph.SetSink(modIndex, numNTermCleavage);
                        var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();
                        protCompositionWithH2O.GetIsotopomerEnvelope();

                        var modCombinations = seqGraph.ModificationParams.GetModificationCombination(modIndex);

                        totalProtCompositions++;
                        for (var charge = minPrecursorIonCharge; charge <= maxPrecursorIonCharge; charge++)
                        {
                            numPrecursorIons++;
                            var precursorIon = new Ion(protCompositionWithH2O, charge);

                            foreach (var ms2ScanNum in run.GetFragmentationSpectraScanNums(precursorIon))
                            {
                                if (run.CheckMs1Signature(precursorIon, ms2ScanNum, precursorTolerance) == false)
                                    continue;

                                numPrecursorIonsPassingFilter++;
                                var spec = run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                                if (spec == null) continue;
                                var scorer = new CorrMatchedPeakCounter(spec, productIonTolerance, minProductIonCharge, maxProductIonCharge);
                                //var scorer = new LikelihoodScorer(scoringModel, spec, productIonTolerance, minProductIonCharge, maxProductIonCharge);
                                var score = seqGraph.GetScore(charge, scorer);

                                if (score <= -50) continue;

                                double existingBestScore;
                                if (bestScorePerScan.TryGetValue(ms2ScanNum, out existingBestScore) &&
                                    score <= existingBestScore) continue;

                                // new best score
                                var sequence = annotation.Substring(numNTermCleavage + 2,
                                    annotation.Length - 4 - numNTermCleavage);
                                var proteinName = targetDb.GetProteinName(offset);
                                var start = targetDb.GetZeroBasedPositionInProtein(offset) + 1;
                                var end = start + sequence.Length - 1;
                                var protLength = targetDb.GetProteinLength(proteinName);
                                bestScorePerScan[ms2ScanNum] = score;
                                bestResultPerScan[ms2ScanNum] = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}",
                                    annotation.Substring(numNTermCleavage + 2, annotation.Length - 4 - numNTermCleavage),   // Sequence
                                    modCombinations,    // Modifications
                                    precursorIon.Composition,   // Composition
                                    (isDecoy ? FastaDatabase.DecoyProteinPrefix + "_" : "") + proteinName,  // ProteinName
                                    targetDb.GetProteinDescription(offset), // ProteinDescription
                                    protLength, // ProteinLength
                                    start,   // Start
                                    end,    // End
                                    charge, // precursorCharge
                                    precursorIon.GetMostAbundantIsotopeMz(),    // MostAbundantIsotopeMz
                                    score);
                            }
                        }
                    }
                }
            }

            // write results into a file
            var icExtension = !isDecoy ? ".ic2result" : "decoy.ic2result";
            var outputFilePath = Path.ChangeExtension(specFilePath, icExtension);
            using (var writer = new StreamWriter(outputFilePath))
            {
//                writer.WriteLine("ScanNum\tAnnotation\tProtein\tProteinDesc\tComposition\tCharge\tBaseIsotopeMz\tScore");
                writer.WriteLine("ScanNum\tSequence\tModifications\tComposition\tProteinName\tProteinDesc\tProteinLength\tStart\tEnd\tCharge\tMostAbundantIsotopeMz\tScore");
                var ms2Scans = new List<int>(bestScorePerScan.Keys);
                ms2Scans.Sort();
                foreach (var ms2ScanNum in bestScorePerScan.OrderByDescending(e => e.Value).Select(scanScorePair => scanScorePair.Key))
                {
                    writer.WriteLine(ms2ScanNum + "\t" + bestResultPerScan[ms2ScanNum]);
                }
            }

            sw.Stop();
            Console.WriteLine("NumProteins: {0}", numProteins);
            Console.WriteLine("NumProteinCompositions: {0}", totalProtCompositions);
            Console.WriteLine("NumPrecursorIons: {0}", numPrecursorIons);
            Console.WriteLine("NumPrecursorIonsWithEvidence: {0}", numPrecursorIonsPassingFilter);

            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void TestModificationSearch()
        {
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\raw\SBEP_STM_001_02272012_Aragon.raw";
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDown\databases\TestMod.fasta";

            // Search parameters
            const int maxNumNTermCleavages = 1;  // 30
            const int maxNumCTermCleavages = 0;
            const int minLength = 20;    // 7
            const int maxLength = 250; // 1000
            const int minPrecursorIonCharge = 3; // 3
            const int maxPrecursorIonCharge = 30;// 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 10;// 10
            const int numMaxModsPerProtein = 2; // 6

            var precursorTolerance = new Tolerance(10);
            var productIonTolerance = new Tolerance(10);

            // Configure amino acid set
            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var cysteinylC = new SearchModification(Modification.Cysteinyl, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                //pyroGluQ,
                dehydroC,
                //cysteinylC,
                //deamdN,
                //deamdQ
                //glutathioneC,
                //oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            TestTopDownSearchMod(dbFilePath, specFilePath, aaSet, minLength, maxLength, maxNumNTermCleavages, maxNumCTermCleavages,
                minPrecursorIonCharge, maxPrecursorIonCharge,
                minProductIonCharge, maxProductIonCharge, precursorTolerance, productIonTolerance, false, false);
        }

        [Test]
        public void TestTopDownSearchMod(
            string dbFilePath, string specFilePath, AminoAcidSet aaSet,
            int minLength, int maxLength,
            int maxNumNTermCleavages, int maxNumCTermCleavages,
            int minPrecursorIonCharge, int maxPrecursorIonCharge,
            int minProductIonCharge, int maxProductIonCharge,
            Tolerance precursorTolerance, Tolerance productIonTolerance,
            bool ultraMod,
            bool isDecoy
            )
        {

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            Console.Write("Reading raw file...");
            //var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 0);
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

            //var scoringModel = new LikelihoodScoringModel(@"C:\cygwin\home\kims336\Data\TopDown\raw\CorrScores_SBEP.txt");
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);


            var targetDb = new FastaDatabase(dbFilePath);
            targetDb.Read();

            var db = !isDecoy ? targetDb : targetDb.Decoy(null, true);  // shuffled decoy

            var indexedDb = new IndexedDatabase(db);

            var annotationsAndOffsets = indexedDb.IntactSequenceAnnotationsAndOffsets(minLength, maxLength, maxNumCTermCleavages);
            //var annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(minLength, maxLength);

            var numProteins = 0;
            long totalProtCompositions = 0;
            long numPrecursorIons = 0;
            long numPrecursorIonsPassingFilter = 0;

            sw.Reset();
            sw.Start();

            var bestScorePerScan = new Dictionary<int, double>();
            var bestResultPerScan = new Dictionary<int, string>();

            foreach (var annotationAndOffset in annotationsAndOffsets)
            {
                ++numProteins;

                var annotation = annotationAndOffset.Annotation;
                var offset = annotationAndOffset.Offset;

                //                    Console.WriteLine(annotation);
                if (numProteins % 100 == 0)
                {
                    Console.Write("Processing {0}{1} proteins...", numProteins,
                        numProteins == 1 ? "st" : numProteins == 2 ? "nd" : numProteins == 3 ? "rd" : "th");
                    if (numProteins != 0)
                    {
                        sw.Stop();
                        sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
                        Console.WriteLine("Elapsed Time: {0:f4} sec", sec);
                        sw.Reset();
                        sw.Start();
                    }
                    //if (numProteins == 10) break;
                }

                //Console.WriteLine(protAnnotation);

                var seqGraph = SequenceGraph.CreateGraph(aaSet, annotation);
                if (seqGraph == null)
                {
                    Console.WriteLine("Ignoring illegal protein: {0}", annotation);
                    continue;
                }

                for (var numNTermCleavage = 0; numNTermCleavage <= maxNumNTermCleavages; numNTermCleavage++)
                {
                    //var compSet = new HashSet<Composition>();
                    var protCompositions = seqGraph.GetSequenceCompositionsWithNTermCleavage(numNTermCleavage);
                    if (ultraMod)
                        Console.WriteLine("#NTermCleavages: {0}, #ProteinCompositions: ", numNTermCleavage);
                    for (var modIndex = 0; modIndex < protCompositions.Length; modIndex++)
                    {
                        if (ultraMod)
                        {
                            if (modIndex % 100 == 0) Console.WriteLine("ModIndex: " + modIndex);
                            //                                if (modIndex >= 100) break;
                        }

                        seqGraph.SetSink(modIndex, numNTermCleavage);
                        var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();

                        var modCombinations = seqGraph.ModificationParams.GetModificationCombination(modIndex);

                        totalProtCompositions++;
                        for (var charge = minPrecursorIonCharge; charge <= maxPrecursorIonCharge; charge++)
                        {
                            numPrecursorIons++;
                            var precursorIon = new Ion(protCompositionWithH2O, charge);

                            foreach (var ms2ScanNum in run.GetFragmentationSpectraScanNums(precursorIon))
                            {
                                if (run.CheckMs1Signature(precursorIon, ms2ScanNum, precursorTolerance) == false)
                                    continue;

                                numPrecursorIonsPassingFilter++;
                                var spec = run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                                if (spec == null) continue;
                                var scorer = new CorrMatchedPeakCounter(spec, productIonTolerance, minProductIonCharge, maxProductIonCharge);
                                //var scorer = new LikelihoodScorer(scoringModel, spec, productIonTolerance, minProductIonCharge, maxProductIonCharge);
                                var scoreAndMods = seqGraph.GetScoreAndModifications(charge, scorer);
                                var score = scoreAndMods.Item1;
                                var mods = scoreAndMods.Item2;

                                if (score <= -50) continue;

                                double existingBestScore;
                                if (bestScorePerScan.TryGetValue(ms2ScanNum, out existingBestScore) &&
                                    score <= existingBestScore) continue;

                                // new best score
                                var sequence = annotation.Substring(numNTermCleavage + 2,
                                    annotation.Length - 4 - numNTermCleavage);
                                var proteinName = targetDb.GetProteinName(offset);
                                var start = targetDb.GetZeroBasedPositionInProtein(offset) + 1;
                                var end = start + sequence.Length - 1;
                                var protLength = targetDb.GetProteinLength(proteinName);
                                bestScorePerScan[ms2ScanNum] = score;
                                bestResultPerScan[ms2ScanNum] = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}",
                                    annotation.Substring(numNTermCleavage + 2, annotation.Length - 4 - numNTermCleavage),   // Sequence
                                    //modCombinations,    // Modifications
                                    "["+mods+"]",
                                    precursorIon.Composition,   // Composition
                                    (isDecoy ? FastaDatabase.DecoyProteinPrefix + "_" : "") + proteinName,  // ProteinName
                                    targetDb.GetProteinDescription(offset), // ProteinDescription
                                    protLength, // ProteinLength
                                    start,   // Start
                                    end,    // End
                                    charge, // precursorCharge
                                    precursorIon.GetMostAbundantIsotopeMz(),    // MostAbundantIsotopeMz
                                    score);
                            }
                        }
                    }
                }
            }

            // write results into a file
            var icExtension = !isDecoy ? ".icresult" : "decoy.icresult";
            var outputFilePath = Path.ChangeExtension(specFilePath, icExtension);
            using (var writer = new StreamWriter(outputFilePath))
            {
                //                writer.WriteLine("ScanNum\tAnnotation\tProtein\tProteinDesc\tComposition\tCharge\tBaseIsotopeMz\tScore");
                writer.WriteLine("ScanNum\tSequence\tModifications\tComposition\tProteinName\tProteinDesc\tProteinLength\tStart\tEnd\tCharge\tMostAbundantIsotopeMz\tScore");
                var ms2Scans = new List<int>(bestScorePerScan.Keys);
                ms2Scans.Sort();
                foreach (var ms2ScanNum in bestScorePerScan.OrderByDescending(e => e.Value).Select(scanScorePair => scanScorePair.Key))
                {
                    writer.WriteLine(ms2ScanNum + "\t" + bestResultPerScan[ms2ScanNum]);
                }
            }

            sw.Stop();
            Console.WriteLine("NumProteins: {0}", numProteins);
            Console.WriteLine("NumProteinCompositions: {0}", totalProtCompositions);
            Console.WriteLine("NumPrecursorIons: {0}", numPrecursorIons);
            Console.WriteLine("NumPrecursorIonsWithEvidence: {0}", numPrecursorIonsPassingFilter);

            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void ExtractProteinSequences()
        {
            const string fastaFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\database\ID_003962_71E1A1D4.fasta";
            var fastaDb = new FastaDatabase(fastaFilePath);
            fastaDb.Read();

            const string proteinFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\database\Proteins.txt";
            var proteins = File.ReadAllLines(proteinFilePath);
            foreach (var protein in proteins)
            {
                var token = protein.Split();
                if (token.Length < 1) continue;
                var proteinId = protein.Split()[0];
                var proteinSequence = fastaDb.GetProteinSequence(proteinId);
                Assert.IsTrue(proteinSequence != null);
                Console.WriteLine(">"+protein);
                Console.WriteLine(proteinSequence);
            }
        }

        [Test]
        public void TestJia()
        {
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1.raw";
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\database\TargetProteins.fasta";

            // Search parameters
            const int maxNumNTermCleavages = 30;  // 30
            const int maxNumCTermCleavages = 0;
            const int minLength = 30;    // 7
            const int maxLength = 300; // 1000
            const int minPrecursorIonCharge = 5; // 3
            const int maxPrecursorIonCharge = 30;// 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 10;// 10
            const int numMaxModsPerProtein = 4; // 6

            var precursorTolerance = new Tolerance(10);
            var productIonTolerance = new Tolerance(10);

            // Configure amino acid set
            var nitrosylC = new SearchModification(Modification.Nitrosyl, 'C', SequenceLocation.Everywhere, false);
            var nethylmaleimideC = new SearchModification(Modification.Nethylmaleimide, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);

            //var cysteinylC = new SearchModification(Modification.CysteinylC, 'C', SequenceLocation.Everywhere, false);
            //var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            //var deamdN = new SearchModification(Modification.Deamidation, 'N', SequenceLocation.Everywhere, false);
            //var deamdQ = new SearchModification(Modification.Deamidation, 'Q', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                nitrosylC,
                nethylmaleimideC,
                dehydroC,
                glutathioneC,
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            TestTopDownSearchMod(dbFilePath, specFilePath, aaSet, minLength, maxLength, maxNumNTermCleavages, maxNumCTermCleavages,
                minPrecursorIonCharge, maxPrecursorIonCharge,
                minProductIonCharge, maxProductIonCharge, precursorTolerance, productIonTolerance, false, false);
        }

        [Test]
        public void TestDirectInfusion()
        {
            // Enolase
            //var protSequence = "MAVSKVYARSVYDSRGNPTVEVELTTEKGVFRSIVPSGASTGVHEALEMRDGDKSKWMGKGVLHAVKNVNDVIAPAFVKANIDVKDQKAVDDFLISLDGTANKSKLGANAILGVSLAASRAAAAEKNVPLYKHLADLSKSKTSPYVLPVPFLNVLNGGSHAGGALALQEFMIAPTGAKTFAEALRIGSEVYHNLKSLTKKRYGASAGNVGDEGGVAPNIQTAEEALDLIVDAIKAAGHDGKVKIGLDCASSEFFKDGKYDLDFKNPNSDKSKWLTGPQLADLYHSLMKRYPIVSIEDPFAEDDWEAWSHFFKTAGIQIVADDLTVTNPKRIATAIEKKAADALLLKVNQIGTLSESIKAAQDSFAAGWGVMVSHRSGETEDTFIADLVVGLRTGQIKTGAPARSERLAKLNQLLRIEEELGDNAVFAGENFHHGDKL";
            // Apo-Transferrin
            var protSequence = "MRLAVGALLVCAVLGLCLAVPDKTVRWCAVSEHEATKCQSFRDHMKSVIPSDGPSVACVKKASYLDCIRAIAANEADAVTLDAGLVYDAYLAPNNLKPVVAEFYGSKEDPQTFYYAVAVVKKDSGFQMNQLRGKKSCHTGLGRSAGWNIPIGLLYCDLPEPRKPLEKAVANFFSGSCAPCADGTDFPQLCQLCPGCGCSTLNQYFGYSGAFKCLKDGAGDVAFVKHSTIFENLANKADRDQYELLCLDNTRKPVDEYKDCHLAQVPSHTVVARSMGGKEDLIWELLNQAQEHFGKDKSKEFQLFSSPHGKDLLFKDSAHGFLKVPPRMDAKMYLGYEYVTAIRNLREGTCPEAPTDECKPVKWCALSHHERLKCDEWSVNSVGKIECVSAETTEDCIAKIMNGEADAMSLDGGFVYIAGKCGLVPVLAENYNKSDNCEDTPEAGYFAVAVVKKSASDLTWDNLKGKKSCHTAVGRTAGWNIPMGLLYNKINHCRFDEFFSEGCAPGSKKDSSLCKLCMGSGLNLCEPNNKEGYYGYTGAFRCLVEKGDVAFVKHQTVPQNTGGKNPDPWAKNLNEKDYELLCLDGTRKPVEEYANCHLARAPNHAVVTRKDKEACVHKILRQQQHLFGSNVTDCSGNFCLFRSETKDLLFRDDTVCLAKLHDRNTYEKYLGEEYVKAVGNLRKCSTSSLLEACTFRRP";
            //protSequence = SimpleStringProcessing.Mutate(SimpleStringProcessing.Shuffle(protSequence), 3);

            var protAnnotation = "_." + protSequence + "._";

            const int modIndex = 0;
            var aminoAcidSet = new AminoAcidSet();
            var seqGraph = SequenceGraph.CreateGraph(aminoAcidSet, protAnnotation);

            // Parameters
            var productIonTolerance = new Tolerance(10);

            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownDirect\ApoTransferrin_240k_N2HCD_2261_ETD_50avg.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

            var seqCompositionArr = seqGraph.GetSequenceCompositions();

            var seqComposition = seqCompositionArr[modIndex];
            Console.WriteLine("Protein mass: {0}", seqComposition.Mass);
            var peptideComposition = seqComposition + Composition.H2O;
            peptideComposition.GetIsotopomerEnvelope();

            Console.WriteLine("Composition: {0}, AveragineMass: {1}", seqComposition, seqComposition.Mass);
            seqGraph.SetSink(modIndex, 0);

            for (var ms2ScanNum = run.MinLcScan; ms2ScanNum <= run.MaxLcScan; ms2ScanNum++)
            {
                var spec = run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                Assert.True(spec != null);
                var scorer = new CorrMatchedPeakCounter(spec, productIonTolerance, 1, 20);
                var score = seqGraph.GetScore(0, scorer);
                Console.WriteLine("{0}\t{1}", ms2ScanNum, score);
            }
        }

        [Test]
        public void TestQcShewSearch()
        {
            // Search parameters
            const int maxNumNTermCleavages = 1;  // 30
            const int maxNumCTermCleavages = 0;
            const int minLength = 7;    // 7
            const int maxLength = 1000; // 1000
            const int minPrecursorIonCharge = 3; // 3
            const int maxPrecursorIonCharge = 40;// 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 10;// 10
            const int numMaxModsPerProtein = 0; // 6

            var precursorTolerance = new Tolerance(15);
            var productIonTolerance = new Tolerance(15);

            //const string dbFilePath = @"..\..\..\TestFiles\BSA.fasta";
            //const string dbFilePath =
            //    @"C:\cygwin\home\kims336\Data\TopDown\ID_003558_56D73071.fasta";
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDownSigma48\Sigma48_UPS2_2012-11-02.fasta";
            //            const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDownSigma48\P01031.fasta";

            //            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownSigma48\Ron_UPS48_2_test_1ug.raw";

            // Configure amino acid set
            //var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.PyroGluQ, 'C', SequenceLocation.Everywhere, false);
            //var cysteinylC = new SearchModification(Modification.CysteinylC, 'C', SequenceLocation.Everywhere, false);
            //var glutathioneC = new SearchModification(Modification.GlutathioneC, 'C', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                //pyroGluQ,
                dehydroC,
                //cysteinylC,
                //glutathioneC,
                oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            TestTopDownSearch(dbFilePath, specFilePath, aaSet, minLength, maxLength, maxNumNTermCleavages, maxNumCTermCleavages,
                minPrecursorIonCharge, maxPrecursorIonCharge,
                minProductIonCharge, maxProductIonCharge, precursorTolerance, productIonTolerance, false, false);
        }

        [Test]
        public void TestHistonSearch()
        {
            const bool isDecoy = true;
            // Search parameters
            const int maxNumNTermCleavages = 1;
            const int maxNumCTermCleavages = 0;
            const int minLength = 7;    // 7
            const int maxLength = 1000; // 1000
            const int minPrecursorIonCharge = 3; // 3
            const int maxPrecursorIonCharge = 40;// 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 10;// 10
            const int numMaxModsPerProtein = 11; // 6

            var precursorTolerance = new Tolerance(15);
            var productIonTolerance = new Tolerance(15);

            //const string dbFilePath = @"..\..\..\TestFiles\BSA.fasta";
            //const string dbFilePath =
            //    @"C:\cygwin\home\kims336\Data\TopDown\ID_003558_56D73071.fasta";
            //const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDownSigma48\Sigma48_UPS2_2012-11-02.fasta";
            const string dbFilePath = @"D:\Research\Data\TopDownHistone\HistoneH4.fasta";

            //            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            const string specFilePath = @"D:\Research\Data\TopDownHistone\071210_070610His0Gy070210H4_H061010A.raw";

            var acetylR = new SearchModification(Modification.Acetylation, 'R', SequenceLocation.Everywhere, false);
            var acetylK = new SearchModification(Modification.Acetylation, 'K', SequenceLocation.Everywhere, false);
            var methylR = new SearchModification(Modification.Methylation, 'R', SequenceLocation.Everywhere, false);
            var methylK = new SearchModification(Modification.Methylation, 'K', SequenceLocation.Everywhere, false);
            var diMethylR = new SearchModification(Modification.DiMethylation, 'R', SequenceLocation.Everywhere, false);
            var diMethylK = new SearchModification(Modification.DiMethylation, 'K', SequenceLocation.Everywhere, false);
            var triMethylR = new SearchModification(Modification.TriMethylation, 'R', SequenceLocation.Everywhere, false);
            var phosphoS = new SearchModification(Modification.Phosphorylation, 'S', SequenceLocation.Everywhere, false);
            var phosphoT = new SearchModification(Modification.Phosphorylation, 'T', SequenceLocation.Everywhere, false);
            var phosphoY = new SearchModification(Modification.Phosphorylation, 'Y', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
                {
                    acetylR, acetylK,
                    methylR, methylK,
                    diMethylR, diMethylK,
                    triMethylR,
                    phosphoS, phosphoT, phosphoY
                };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            TestTopDownSearch(dbFilePath, specFilePath, aaSet, minLength, maxLength, maxNumNTermCleavages, maxNumCTermCleavages,
                minPrecursorIonCharge, maxPrecursorIonCharge,
                minProductIonCharge, maxProductIonCharge, precursorTolerance, productIonTolerance, true, isDecoy);            
        }

        [Test]
        public void TestSigma48Search()
        {
            // Search parameters
            const int maxNumNTermCleavages = 1;  // 30
            const int maxNumCTermCleavages = 0;

            const int minLength = 7;    // 7
            const int maxLength = 1000; // 1000
            const int minPrecursorIonCharge = 3; // 3
            const int maxPrecursorIonCharge = 40;// 67
            const int minProductIonCharge = 1; // 1
            const int maxProductIonCharge = 10;// 10
            const int numMaxModsPerProtein = 0; // 6

            var precursorTolerance = new Tolerance(10);
            var productIonTolerance = new Tolerance(10);

            //const string dbFilePath = @"..\..\..\TestFiles\BSA.fasta";
            //const string dbFilePath =
            //    @"C:\cygwin\home\kims336\Data\TopDown\ID_003558_56D73071.fasta";
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDownSigma48\Sigma48_UPS2_2012-11-02.fasta";
//            const string dbFilePath = @"C:\cygwin\home\kims336\Data\TopDownSigma48\P01031.fasta";

            //            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDownSigma48\Ron_UPS48_2_test_1ug.raw";

            // Configure amino acid set
            //var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            //var cysteinylC = new SearchModification(Modification.CysteinylC, 'C', SequenceLocation.Everywhere, false);
            //var glutathioneC = new SearchModification(Modification.GlutathioneC, 'C', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                //pyroGluQ,
                dehydroC,
                //cysteinylC,
                //glutathioneC,
                oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            TestTopDownSearch(dbFilePath, specFilePath, aaSet, minLength, maxLength, maxNumNTermCleavages, maxNumCTermCleavages,
                minPrecursorIonCharge, maxPrecursorIonCharge,
                minProductIonCharge, maxProductIonCharge, precursorTolerance, productIonTolerance, false, false);
        }

        [Test]
        public void TestComputingIsotopomerEnvelop()
        {
            const string sequence =
                "MRLNTLSPAEGSKKAGKRLGRGIGSGLGKTGGRGHKGQKSRSGGGVRRGFEGGQMPLYRRLPKFGFTSRKAAITAEIRLSDLAKVEGGVVDLNTLKAANIIGIQIEFAKVILAGEVTTPVTVRGLRVTKGARAAIEAAGGKIEE";
            const int charge = 17;

            var aaSet = new AminoAcidSet();
            var composition = aaSet.GetComposition(sequence);

            var precursorIon = new Ion(composition + Composition.H2O, charge);

            Console.WriteLine("{0}\t{1}\t{2}", precursorIon.Composition, charge, precursorIon.GetMonoIsotopicMz());
            foreach (var isotope in precursorIon.GetIsotopes(relativeIntensityThreshold: 0.05))
            {
                Console.WriteLine("{0}: {1}\t{2}", isotope.Index, precursorIon.GetIsotopeMz(isotope.Index), isotope.Ratio);
            }
        }

        [Test]
        public void TestTopDownSearchOneProtein()
        {
            // Parameters
            const int minCharge = 3;    // 3
            const int maxCharge = 67;   // 67
            var precursorIonTolerance = new Tolerance(10);
            var productIonTolerance = new Tolerance(10);

            var sw = new System.Diagnostics.Stopwatch();

            // Configure amino acids
            //var aaSet = new AminoAcidSet();


            // Configure amino acid set
            const int numMaxModsPerProtein = 6;
            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            var dehydro = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var cysteinylC = new SearchModification(Modification.Cysteinyl, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                pyroGluQ,
                //dehydro,
                //cysteinylC,
                //glutathioneC,
                //oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            const string protAnnotation = "A.HAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK._";
            //const string protAnnotation =
            //    "_.MKLYNLKDHNEQVSFAQAVTQGLGKNQGLFFPHDLPEFSLTEIDEMLKLDFVTRSAKILSAFIGDEIPQEILEERVRAAFAFPAPVANVESDVGCLELFHGPTLAFKDFGGRFMAQMLTHIAGDKPVTILTATSGDTGAAVAHAFYGLPNVKVVILYPRGKISPLQEKLFCTLGGNIETVAIDGDFDACQALVKQAFDDEELKVALGLNSANSINISRLLAQICYYFEAVAQLPQETRNQLVVSVPSGNFGDLTAGLLAKSLGLPVKRFIAATNVNDTVPRFLHDGQWSPKATQATLSNAMDVSQPNNWPRVEELFRRKIWQLKELGYAAVDDETTQQTMRELKELGYTSEPHAAVAYRALRDQLNPGEYGLFLGTAHPAKFKESVEAILGETLDLPKELAERADLPLLSHNLPADFAALRKLMMNHQ._";

            // Create a sequence graph
            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            if (seqGraph == null)
            {
                Console.WriteLine("Invalid sequence: {0}", protAnnotation);
                return;
            }

            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

            sw.Start();
            var precursorFilter = new PrecursorFilter(run, precursorIonTolerance);

            var seqCompositionArr = seqGraph.GetSequenceCompositions();
            Console.WriteLine("Length: {0}\tNumCompositions: {1}", protAnnotation.Length-4, seqCompositionArr.Length);

            for (var modIndex = 0; modIndex < seqCompositionArr.Length; modIndex++)
            {
                var seqComposition = seqCompositionArr[modIndex];
                var peptideComposition = seqComposition + Composition.H2O;
                peptideComposition.GetIsotopomerEnvelope();

                Console.WriteLine("Composition: {0}, AveragineMass: {1}", seqComposition, seqComposition.Mass);

                //if (Math.Abs(seqComposition.GetMass() - 47162.1844822) > 0.001) continue;

                for (var charge = minCharge; charge <= maxCharge; charge++)
                {
                    var precursorIon = new Ion(peptideComposition, charge);

                    var bestScore = Double.NegativeInfinity;
                    var bestScanNum = -1;
                    foreach (var ms2ScanNum in run.GetFragmentationSpectraScanNums(precursorIon))
                    {
                        if (!precursorFilter.IsValid(precursorIon, ms2ScanNum)) continue;

                        var spec = run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                        if (spec == null) continue;
                        var scorer = new MatchedPeakCounter(spec, productIonTolerance, 1, 6);
                        var score = seqGraph.GetScore(charge, scorer);

                        //Console.WriteLine("{0}\t{1}\t{2}\t{3}", precursorIon.GetMostAbundantIsotopeMz(), precursorIon.precursorCharge, ms2ScanNum, score);
                        if (score > bestScore)
                        {
                            bestScore = score;
                            bestScanNum = ms2ScanNum;
                        }
                    }
                    if (bestScore > 10)
                    {
                        Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", protAnnotation, charge, precursorIon.GetMostAbundantIsotopeMz(), bestScanNum, bestScore);
                    }
                }
            }
            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

        }


        [Test]
        public void CountProteins()
        {
            const string filePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.icresult";

            var protSet = new HashSet<string>();
            foreach (var line in File.ReadLines(filePath))
            {
                var token = line.Split('\t');
                if (token.Length != 4) continue;
                protSet.Add(token[0]);
            }
            Console.WriteLine("NumProteins: {0}", protSet.Count);
        }

        [Test]
        public void TestGeneratingAllXics()
        {
            // Search parameters
            const int numNTermCleavages = 1;  // 30
            const int minLength = 7;
            const int maxLength = 1000;
            const int minCharge = 3; // 3
            const int maxCharge = 67; // 67
            const int numMaxModsPerProtein = 0; // 6
            var precursorTolerance = new Tolerance(20);
            //const string dbFilePath = @"..\..\..\TestFiles\BSA.fasta";
            const string dbFilePath =
                @"C:\cygwin\home\kims336\Data\TopDown\ID_003558_56D73071.fasta";

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            Console.Write("Reading raw file...");
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

            // Configure amino acid set
            //            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.ProteinNTerm, false);
            var dehydro = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var cysteinylC = new SearchModification(Modification.Cysteinyl, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            //            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
                {
                    //pyroGluQ,
                    dehydro,
                    cysteinylC,
                    glutathioneC,
                    //oxM
                };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var targetDb = new FastaDatabase(dbFilePath);
            var indexedDb = new IndexedDatabase(targetDb);

            var numProteins = 0;
            long totalProtCompositions = 0;
            long numXics = 0;

            sw.Reset();
            sw.Start();
            Console.WriteLine("Generating XICs...");

            foreach (var protAnnotationAndOffset in indexedDb.IntactSequenceAnnotationsAndOffsets(minLength, maxLength))
            {
                ++numProteins;

                if (numProteins % 1000 == 0)
                {
                    Console.WriteLine("Processed {0} proteins", numProteins);
                }

                //Console.WriteLine(protAnnotation);

                var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotationAndOffset.Annotation);
                if (seqGraph == null) continue;

                for (var nTermCleavages = 0; nTermCleavages <= numNTermCleavages; nTermCleavages++)
                {
                    var protCompositions = seqGraph.GetSequenceCompositionsWithNTermCleavage(nTermCleavages);
                    foreach (var protComposition in protCompositions)
                    {
                        totalProtCompositions++;
                        var mostAbundantIsotopeIndex = protComposition.GetMostAbundantIsotopeZeroBasedIndex();
                        for (var charge = minCharge; charge <= maxCharge; charge++)
                        {
                            numXics++;
                            var precursorIon = new Ion(protComposition + Composition.H2O, charge);
                            run.GetExtractedIonChromatogram(precursorIon.GetIsotopeMz(mostAbundantIsotopeIndex), precursorTolerance);
                        }
                    }
                }
            }

            sw.Stop();
            Console.WriteLine("NumProteins: {0}", numProteins);
            Console.WriteLine("NumProteinCompositions: {0}", totalProtCompositions);
            Console.WriteLine("NumXics: {0}", numXics);
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }
        
        [Test]
        public void TestGeneratingXics()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            const string specFilePath = @"C:\workspace\TopDown\E_coli_iscU_60_mock.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);
            const string protAnnotation = "A.HAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK.-";
            var aaSet = new AminoAcidSet();

            var precursorTolerance = new Tolerance(10);

            // Create a sequence graph
            var protSeq = protAnnotation.Substring(2, protAnnotation.Length - 4);
            var seqGraph = SequenceGraph.CreateGraph(aaSet, protSeq);
            foreach (var protComposition in seqGraph.GetSequenceCompositions())
            {
                var mostAbundantIsotopeIndex = protComposition.GetMostAbundantIsotopeZeroBasedIndex();
                Console.WriteLine("Composition\t{0}", protComposition);
                Console.WriteLine("MostAbundantIsotopeIndex\t{0}", mostAbundantIsotopeIndex);
                Console.WriteLine();

                for (var charge = 10; charge <= 14; charge++)
                {
                    var precursorIon = new Ion(protComposition+Composition.H2O, charge);
                    var xic = run.GetExtractedIonChromatogram(precursorIon.GetIsotopeMz(mostAbundantIsotopeIndex), precursorTolerance);
                    //Console.WriteLine(xic[0].ScanNum + " " + xic[1].ScanNum);
                    
                    Console.WriteLine("ScanNum\t{0}", string.Join("\t", xic.Select(p => p.ScanNum.ToString())));
                    Console.WriteLine("precursorCharge " + charge + "\t" + string.Join("\t", xic.Select(p => p.Intensity.ToString())));
                }

                Console.WriteLine("\nCharge\tm/z");
                for (var charge = 9; charge <= 18; charge++)
                {
                    var precursorIon = new Ion(protComposition + Composition.H2O, charge);
                    Console.WriteLine("{0}\t{1}", charge, precursorIon.GetIsotopeMz(mostAbundantIsotopeIndex));
                }
            }

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }



        [Test]
        public void TestHistonEnumeration()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            const int numNTermCleavages = 0;
            const int numMaxModsPerProtein = 11;
            const string protAnnotation = "_.MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG._";  // Histone H4
            //const string protAnnotation =
            //    "_.MARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA._"; // Histone H3.2

            var acetylR = new SearchModification(Modification.Acetylation, 'R', SequenceLocation.Everywhere, false);
            var acetylK = new SearchModification(Modification.Acetylation, 'K', SequenceLocation.Everywhere, false);
            var methylR = new SearchModification(Modification.Methylation, 'R', SequenceLocation.Everywhere, false);
            var methylK = new SearchModification(Modification.Methylation, 'K', SequenceLocation.Everywhere, false);
            var diMethylR = new SearchModification(Modification.DiMethylation, 'R', SequenceLocation.Everywhere, false);
            var diMethylK = new SearchModification(Modification.DiMethylation, 'K', SequenceLocation.Everywhere, false);
            var triMethylR = new SearchModification(Modification.TriMethylation, 'R', SequenceLocation.Everywhere, false);
            var phosphoS = new SearchModification(Modification.Phosphorylation, 'S', SequenceLocation.Everywhere, false);
            var phosphoT = new SearchModification(Modification.Phosphorylation, 'T', SequenceLocation.Everywhere, false);
            var phosphoY = new SearchModification(Modification.Phosphorylation, 'Y', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
                {
                    acetylR, acetylK,
                    methylR, methylK,
                    diMethylR, diMethylK,
                    triMethylR,
                    phosphoS, phosphoT, phosphoY
                };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var numPossiblyModifiedResidues = 0;
            var numR = 0;
            var numK = 0;
            var numSTY = 0;
            foreach (var aa in protAnnotation.Substring(2, protAnnotation.Length - 4))
            {
                var numMods = aaSet.GetModificationIndices(aa).Length;
                if (aa == 'S' || aa == 'T' || aa == 'Y') numSTY++;
                if (aa == 'R') numR++;
                if (aa == 'K') numK++;
                if (numMods >= 1)
                {
                    numPossiblyModifiedResidues += 1;
                }
            }

            var numProteoforms = 0.0;
            for (var numMods = 0; numMods <= numMaxModsPerProtein; numMods++)
            {
                for (var i = 0; i <= numMods; i++)
                {
                    for (var j = 0; i + j <= numMods; j++)
                    {
                        var k = numMods - i - j;
                        numProteoforms += SpecialFunctions.Binomial(numR, i)*Math.Pow(4, i)
                                          *SpecialFunctions.Binomial(numK, j)*Math.Pow(3, j)
                                          *SpecialFunctions.Binomial(numSTY, k);
                    }
                }
            }
            Console.WriteLine("#Proteoforms: {0:E2}", numProteoforms);
            Console.WriteLine("#PossiblyModifiedResidues: {0}", numPossiblyModifiedResidues);
            Console.WriteLine("#STY: {0}", numSTY);
            Console.WriteLine("#K: {0}", numK);
            Console.WriteLine("#R: {0}", numR);
            Console.WriteLine("5 choose 2: {0}", SpecialFunctions.Binomial(5, 2));

            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            if (seqGraph == null)
            {
                Console.WriteLine("Invalid sequence: {0}", protAnnotation);
                return;
            }

            Console.WriteLine("Num sequence compositions: {0}, {1}", seqGraph.GetNumSequenceCompositions(), seqGraph.GetNumDistinctSequenceCompositions()
                );

            Console.WriteLine("Num product compositions: {0}", seqGraph.GetNumFragmentCompositions()
                );

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void TestTopDownEnumeration()
        {
            // Search parameters
            const int numNTermCleavages = 30;
            const int minLength = 7;
            const int maxLength = 1000;
            const int numMaxModsPerProtein = 6;
            //const string dbFilePath = @"..\..\..\TestFiles\BSA.fasta";
            const string dbFilePath =
                @"..\..\..\TestFiles\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            // Configure amino acid set
//            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.ProteinNTerm, false);
            var dehydro = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var cysteinylC = new SearchModification(Modification.Cysteinyl, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
//            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
                {
                    //pyroGluQ,
                    dehydro,
                    cysteinylC,
                    glutathioneC,
                    //oxM
                };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var targetDb = new FastaDatabase(dbFilePath);
            var indexedDb = new IndexedDatabase(targetDb);

            var numProteins = 0;
            long totalProtCompositions = 0;
            foreach (var protAnnotationAndOffset in indexedDb.IntactSequenceAnnotationsAndOffsets(minLength, maxLength))
            {
                ++numProteins;

                if (numProteins % 1000 == 0)
                {
                    Console.WriteLine("Processed {0} proteins", numProteins);
                }

                var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotationAndOffset.Annotation);
                if (seqGraph == null) continue;

                for (var nTermCleavage = 0; nTermCleavage <= numNTermCleavages; nTermCleavage++)
                {
                    totalProtCompositions += seqGraph.GetNumSequenceCompositionsWithNTermCleavage(nTermCleavage);
                }
            }

            sw.Stop();
            Console.WriteLine("NumProteins: {0}", numProteins);
            Console.WriteLine("NumProteinCompositions: {0}", totalProtCompositions);
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }
    }
}
