using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.Scoring;

namespace InformedProteomics.TopDown.Execution
{
    public class IcTopDown
    {
        public IcTopDown(
            string dbFilePath,
            string specFilePath,
            AminoAcidSet aaSet)
        {

        }

        public IcTopDown(
            string dbFilePath, 
            string specFilePath, 
            AminoAcidSet aaSet,
            int minProteinLength, 
            int maxProteinLength, 
            int maxNumNTermCleavages,
            int maxNumCTermCleavages,
            int minPrecursorIonCharge, 
            int maxPrecursorIonCharge,
            int minProductIonCharge, 
            int maxProductIonCharge,
            Tolerance precursorIonTolerance, 
            Tolerance productIonTolerance,
            bool runTargetDecoyAnalysis)
        {
            
        }

        public void TestTopDownSearch(
            string dbFilePath, string specFilePath, AminoAcidSet aaSet,
            int minLength, int maxLength, int maxNumNTermCleavages,
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
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 0);
            //var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

            var scoringModel = new LikelihoodScoringModel(@"C:\cygwin\home\kims336\Data\TopDown\raw\cidCosineMatches.txt");
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);


            var targetDb = new FastaDatabase(dbFilePath);
            targetDb.Read();

            var db = !isDecoy ? targetDb : targetDb.Decoy(null, true);  // shuffled decoy

            var indexedDb = new IndexedDatabase(db);

            var annotationsAndOffsets = indexedDb.IntactSequenceAnnotationsAndOffsets(minLength, maxLength);

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
                        sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
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
                        Console.WriteLine("#NTermCleavages: {0}, #ProteinCompositions: ", numNTermCleavage,
                            protCompositions.Length);
                    for (var modIndex = 0; modIndex < protCompositions.Length; modIndex++)
                    {
                        if (ultraMod)
                        {
                            if (modIndex % 100 == 0) Console.WriteLine("ModIndex: " + modIndex);
                            //                                if (modIndex >= 100) break;
                        }

                        seqGraph.SetSink(modIndex, numNTermCleavage);
                        var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();
                        //if (compSet.Contains(protCompositionWithH2O)) continue;
                        //compSet.Add(protCompositionWithH2O);

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
                                //var scorer = new MatchedPeakCounter(spec, productIonTolerance, minProductIonCharge, maxProductIonCharge);
                                var scorer = new LikelihoodScorer(scoringModel, spec, productIonTolerance,
                                    minProductIonCharge, maxProductIonCharge);
                                var score = seqGraph.GetScore(charge, scorer);

                                if (score <= 0) continue;

                                double existingBestScore;
                                if (bestScorePerScan.TryGetValue(ms2ScanNum, out existingBestScore) &&
                                    score <= existingBestScore) continue;

                                // new best score
                                bestScorePerScan[ms2ScanNum] = score;
                                bestResultPerScan[ms2ScanNum] = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}",
                                    annotation.Substring(numNTermCleavage + 2,
                                        annotation.Length - 4 - numNTermCleavage),
                                    (isDecoy ? "XXX_" : "") + targetDb.GetProteinName(offset),
                                    targetDb.GetProteinDescription(offset),
                                    precursorIon.Composition,
                                    charge,
                                    precursorIon.GetMostAbundantIsotopeMz(),
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
                writer.WriteLine("ScanNum\tAnnotation\tProtein\tProteinDesc\tComposition\tCharge\tBaseIsotopeMz\tScore");
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
    }
}
