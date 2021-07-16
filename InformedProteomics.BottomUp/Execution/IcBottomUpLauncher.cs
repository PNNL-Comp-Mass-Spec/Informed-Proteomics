using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.SearchResults;
using InformedProteomics.BottomUp.Scoring;
using InformedProteomics.TopDown.Execution;
using InformedProteomics.TopDown.Scoring;
using PRISM;

namespace InformedProteomics.BottomUp.Execution
{
    public class IcBottomUpLauncher
    {
        public const string TargetFileExtension = "_IcTarget.tsv";
        public const string DecoyFileExtension = "_IcDecoy.tsv";
        public const string TdaFileExtension = "_IcTda.tsv";

        public IcBottomUpLauncher(
            string specFilePath,
            string dbFilePath,
            string outputDir,
            AminoAcidSet aaSet,
            Enzyme enzyme)
        {
            ErrorMessage = string.Empty;

            SpecFilePath = specFilePath;
            DatabaseFilePath = dbFilePath;
            AminoAcidSet = aaSet;
            Enzyme = enzyme;

            if (outputDir == null)
            {
                OutputDir = Path.GetDirectoryName(SpecFilePath);
            }
            else
            {
                if (!Directory.Exists(outputDir))
                {
                    if (File.Exists(outputDir) && !File.GetAttributes(outputDir).HasFlag(FileAttributes.Directory))
                    {
                        throw new Exception(outputDir + " is not a directory!");
                    }
                    Directory.CreateDirectory(outputDir);
                }
                OutputDir = outputDir;
            }

            OutputDir = outputDir;
            MinSequenceLength = 6;
            MaxSequenceLength = 30;
            MinPrecursorIonCharge = 1;
            MaxPrecursorIonCharge = 4;
            MinProductIonCharge = 1;
            MaxProductIonCharge = 3;
            PrecursorIonTolerance = new Tolerance(10);
            ProductIonTolerance = new Tolerance(10);
            RunTargetDecoyAnalysis = DatabaseSearchMode.Both;
            NumTolerableTermini = 1;
            MatchesPerSpectrumToReport = 10;
        }

        public string ErrorMessage { get; private set; }
        public string SpecFilePath { get; }
        public string DatabaseFilePath { get; }
        public string OutputDir { get; }
        public AminoAcidSet AminoAcidSet { get; }
        public Enzyme Enzyme { get; }

        /// <summary>
        /// Min Sequence Length
        /// </summary>
        /// <remarks>default 6</remarks>
        public int MinSequenceLength { get; set; }

        /// <summary>
        /// Max Sequence Length
        /// </summary>
        /// <remarks>default 30</remarks>
        public int MaxSequenceLength { get; set; }

        /// <summary>
        /// Min Precursor Ion Charge
        /// </summary>
        /// <remarks>default 1</remarks>
        public int MinPrecursorIonCharge { get; set; }

        /// <summary>
        /// Max Precursor Ion Charge
        /// </summary>
        /// <remarks>default 4</remarks>
        public int MaxPrecursorIonCharge { get; set; }

        /// <summary>
        /// Min Product Ion Charge
        /// </summary>
        /// <remarks>default 1</remarks>
        public int MinProductIonCharge { get; set; }

        /// <summary>
        /// Max Product Ion Charge
        /// </summary>
        /// <remarks>default 3</remarks>
        public int MaxProductIonCharge { get; set; }

        /// <summary>
        /// Precursor Ion Tolerance
        /// </summary>
        /// <remarks>default 10 ppm</remarks>
        public Tolerance PrecursorIonTolerance { get; set; }

        /// <summary>
        /// Product Ion Tolerance
        /// </summary>
        /// <remarks>default 10 ppm</remarks>
        public Tolerance ProductIonTolerance { get; set; }

        /// <summary>
        /// Run Target Decoy Analysis boolean (reads/updates RunTargetDecoyAnalysis)
        /// </summary>
        /// <remarks>
        /// default true
        /// true: target and decoy, false: target only, null: decoy only
        /// </remarks>
        public bool? RunTargetDecoyAnalysisBool
        {
            get
            {
                if (RunTargetDecoyAnalysis == DatabaseSearchMode.Both)
                {
                    return true;
                }

                if (RunTargetDecoyAnalysis == DatabaseSearchMode.Decoy)
                {
                    return null;
                }
                //(Tda2 == DatabaseSearchMode.Target)
                return false;
            }
            set
            {
                if (value == null)
                {
                    RunTargetDecoyAnalysis = DatabaseSearchMode.Decoy;
                }
                else if (value.Value)
                {
                    RunTargetDecoyAnalysis = DatabaseSearchMode.Both;
                }
                else
                {
                    RunTargetDecoyAnalysis = DatabaseSearchMode.Target;
                }
            }
        }

        /// <summary>
        /// Run Target Decoy Analysis
        /// </summary>
        /// <remarks>default Both</remarks>
        public DatabaseSearchMode RunTargetDecoyAnalysis { get; set; }

        /// <summary>
        /// Num Tolerable Termini
        /// </summary>
        /// <remarks>default 1</remarks>
        public int NumTolerableTermini { get; set; }

        /// <summary>
        /// Matches per spectrum to report in the results
        /// </summary>
        /// <remarks>default 10</remarks>
        public int MatchesPerSpectrumToReport { get; set; }

        /// <summary>
        /// Max Number of Threads
        /// </summary>
        /// <remarks>default 4</remarks>
        public int MaxNumThreads { get; set; } // TODO: Actually implement multithreading....

        /// <summary>
        /// Precursor Ion Tolerance, in ppm
        /// </summary>
        /// <remarks>default 10 ppm</remarks>
        public double PrecursorIonTolerancePpm
        {
            get => PrecursorIonTolerance.GetValue();
            set => PrecursorIonTolerance = new Tolerance(value);
        }

        /// <summary>
        /// Product Ion Tolerance, in ppm
        /// </summary>
        /// <remarks>default 10 ppm</remarks>
        public double ProductIonTolerancePpm
        {
            get => ProductIonTolerance.GetValue();
            set => ProductIonTolerance = new Tolerance(value);
        }

        private LcMsRun _run;
        private ProductScorerBasedOnDeconvolutedSpectra _ms2ScorerFactory;
        private InformedBottomUpScorer _bottomUpScorer;

        public bool RunSearch(double corrThreshold)
        {
            var sw = new Stopwatch();
            ErrorMessage = string.Empty;

            Console.Write("Reading raw file...");
            sw.Start();
            _run = InMemoryLcMsRun.GetLcMsRun(SpecFilePath, 1.4826, 1.4826);
            _bottomUpScorer = new InformedBottomUpScorer(_run, AminoAcidSet, MinProductIonCharge, MaxProductIonCharge, ProductIonTolerance);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
            Console.WriteLine("Elapsed Time: {0:f4} sec", sec);

            sw.Reset();
            Console.Write("Determining precursor masses...");
            sw.Start();
            var ms1Filter = new Ms1IsotopeAndChargeCorrFilter(_run, PrecursorIonTolerance, MinPrecursorIonCharge, MaxPrecursorIonCharge,
                400, 5000, corrThreshold, 0, 0);
            sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
            Console.WriteLine("Elapsed Time: {0:f4} sec", sec);

            sw.Reset();
            Console.Write("Deconvoluting MS2 spectra...");
            sw.Start();
            _ms2ScorerFactory = new ProductScorerBasedOnDeconvolutedSpectra(
                _run,
                MinProductIonCharge, MaxProductIonCharge,
                new Tolerance(10),
                0
                );
            _ms2ScorerFactory.DeconvoluteAllProductSpectra();
            sw.Stop();
            sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
            Console.WriteLine("Elapsed Time: {0:f4} sec", sec);

            // Target database
            var targetDb = new FastaDatabase(DatabaseFilePath);

            //            string dirName = OutputDir ?? Path.GetDirectoryName(SpecFilePath);

            var baseName = Path.GetFileNameWithoutExtension(SpecFilePath);
            var targetOutputFilePath = Path.Combine(OutputDir, baseName + TargetFileExtension);
            var decoyOutputFilePath = Path.Combine(OutputDir, baseName + DecoyFileExtension);
            var tdaOutputFilePath = Path.Combine(OutputDir, baseName + TdaFileExtension);

            if (RunTargetDecoyAnalysis.HasFlag(DatabaseSearchMode.Target))
            {
                sw.Reset();
                Console.Write("Reading the target database...");
                sw.Start();
                targetDb.Read();
                sw.Stop();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine("Elapsed Time: {0:f4} sec", sec);

                sw.Reset();
                Console.WriteLine("Searching the target database");
                sw.Start();
                var targetMatches = RunSearch(GetAnnotationsAndOffsets(targetDb), ms1Filter, false);
                sw.Stop();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine("Target database search elapsed time: {0:f4} sec", sec);

                sw.Reset();
                Console.Write("Rescoring and writing target results...");
                sw.Start();
                WriteResultsToFile(targetMatches, targetOutputFilePath, targetDb);
                sw.Stop();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine("Elapsed time: {0:f4} sec", sec);
            }

            if (RunTargetDecoyAnalysis.HasFlag(DatabaseSearchMode.Decoy))
            {
                // Decoy database
                sw.Reset();
                Console.Write("Reading the decoy database...");
                sw.Start();
                var decoyDb = targetDb.Decoy(Enzyme);
                decoyDb.Read();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine("Elapsed Time: {0:f4} sec", sec);

                sw.Reset();
                Console.WriteLine("Searching the decoy database");
                sw.Start();
                var decoyMatches = RunSearch(GetAnnotationsAndOffsets(decoyDb), ms1Filter, true);
                sw.Stop();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine("Decoy database search elapsed Time: {0:f4} sec", sec);

                sw.Reset();
                Console.Write("Rescoring and writing decoy results...");
                sw.Start();
                WriteResultsToFile(decoyMatches, decoyOutputFilePath, decoyDb);
                sw.Stop();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine("Elapsed time: {0:f4} sec", sec);
            }

            if (RunTargetDecoyAnalysis.HasFlag(DatabaseSearchMode.Both))
            {
                var fdrCalculator = new FdrCalculator(targetOutputFilePath, decoyOutputFilePath);
                if (fdrCalculator.HasError())
                {
                    ErrorMessage = fdrCalculator.ErrorMessage;
                    Console.WriteLine("Error computing FDR: " + fdrCalculator.ErrorMessage);
                    return false;
                }

                fdrCalculator.WriteTo(tdaOutputFilePath);
            }

            Console.WriteLine("Done");
            return true;
        }

        private IEnumerable<AnnotationAndOffset> GetAnnotationsAndOffsets(FastaDatabase database)
        {
            var indexedDbTarget = new IndexedDatabase(database);

            if (NumTolerableTermini == 0)
            {
                return indexedDbTarget.AnnotationsAndOffsetsNoEnzyme(MinSequenceLength, MaxSequenceLength);
            }

            return indexedDbTarget.AnnotationsAndOffsets(MinSequenceLength, MaxSequenceLength, NumTolerableTermini, 2, Enzyme);
        }

        private SortedSet<DatabaseSequenceSpectrumMatch>[] RunSearch(IEnumerable<AnnotationAndOffset> annotationsAndOffsets, ISequenceFilter ms1Filter, bool isDecoy)
        {
            var sw = new Stopwatch();
            var numPeptides = 0;
            sw.Reset();
            sw.Start();

            var matches = new SortedSet<DatabaseSequenceSpectrumMatch>[_run.MaxLcScan + 1];

            // TODO: N-term Met cleavage
            foreach (var annotationAndOffset in annotationsAndOffsets)
            {
                ++numPeptides;

                var annotation = annotationAndOffset.Annotation;
                var offset = annotationAndOffset.Offset;

                if (numPeptides % 100000 == 0)
                {
                    Console.Write("Processing {0}{1} peptides...", numPeptides,
                        numPeptides == 1 ? "st" : numPeptides == 2 ? "nd" : numPeptides == 3 ? "rd" : "th");
                    if (numPeptides != 0)
                    {
                        sw.Stop();
                        var sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                        Console.WriteLine("Elapsed Time: {0:f4} sec", sec);
                        sw.Reset();
                        sw.Start();
                    }
                }

                var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, annotation);
                if (seqGraph == null)
                {
                    // ConsoleMsgUtils.ShowWarning(string.Format("Ignoring illegal protein: {0}", annotation));
                    continue;
                }

                //var protCompositions = seqGraph.GetSequenceCompositions();
                var numProteoforms = seqGraph.GetNumProteoformCompositions();
                var modCombs = seqGraph.GetModificationCombinations();
                for (var modIndex = 0; modIndex < numProteoforms; modIndex++)
                {
                    seqGraph.SetSink(modIndex);
                    var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();
                    var sequenceMass = protCompositionWithH2O.Mass;
                    var modCombinations = modCombs[modIndex];

                    foreach (var ms2ScanNum in ms1Filter.GetMatchingMs2ScanNums(sequenceMass))
                    {
                        if (_run.GetSpectrum(ms2ScanNum) is not ProductSpectrum spec)
                        {
                            continue;
                        }

                        var charge =
                            (int)Math.Round(sequenceMass / (spec.IsolationWindow.IsolationWindowTargetMz - Constants.Proton));
                        var scorer = _ms2ScorerFactory.GetMs2Scorer(ms2ScanNum);
                        var score = seqGraph.GetFragmentScore(scorer);
                        if (score <= 2)
                        {
                            continue;
                        }

                        var precursorIon = new Ion(protCompositionWithH2O, charge);
                        var sequence = annotation.Substring(2, annotation.Length - 4);
                        var pre = annotation[0];
                        var post = annotation[annotation.Length - 1];
                        var prsm = new DatabaseSequenceSpectrumMatch(sequence, pre, post, ms2ScanNum, offset, 0, modCombinations,
                            precursorIon, score, isDecoy);

                        if (matches[ms2ScanNum] == null)
                        {
                            matches[ms2ScanNum] = new SortedSet<DatabaseSequenceSpectrumMatch> { prsm };
                        }
                        else // already exists
                        {
                            var existingMatches = matches[ms2ScanNum];
                            if (existingMatches.Count < MatchesPerSpectrumToReport)
                            {
                                existingMatches.Add(prsm);
                            }
                            else
                            {
                                var minScore = existingMatches.Min.Score;
                                if (score > minScore)
                                {
                                    existingMatches.Add(prsm);
                                    existingMatches.Remove(existingMatches.Min);
                                }
                            }
                        }
                    }
                }
            }

            return matches;
        }

        private void WriteResultsToFile(IReadOnlyList<SortedSet<DatabaseSequenceSpectrumMatch>> matches, string outputFilePath, FastaDatabase database)
        {
            using var writer = new StreamWriter(outputFilePath);

            writer.WriteLine("Scan\tPre\tSequence\tPost\tModifications\tComposition\tProteinName\tProteinDesc" +
                             "\tProteinLength\tStart\tEnd\tCharge\tMostAbundantIsotopeMz\tMass\t#MatchedFragments\tIcScore"
            );

            for (var scanNum = _run.MinLcScan; scanNum <= _run.MaxLcScan; scanNum++)
            {
                if (matches[scanNum] == null)
                {
                    continue;
                }

                foreach (var match in matches[scanNum].Reverse())
                {
                    var sequence = match.Sequence;
                    var offset = match.Offset;
                    var start = database.GetOneBasedPositionInProtein(offset) + 1 + match.NumNTermCleavages;
                    var end = start + sequence.Length - 1;
                    var proteinName = database.GetProteinName(match.Offset);
                    var protLength = database.GetProteinLength(proteinName);
                    var ion = match.Ion;

                    var scores = _bottomUpScorer.GetScores(match, ion.Composition, ion.Charge, scanNum);

                    // ReSharper disable once ConditionIsAlwaysTrueOrFalse
                    // ReSharper disable HeuristicUnreachableCode
                    if (ion == null)
                    {
                        ConsoleMsgUtils.ShowWarning("Null ion in WriteResultsToFile");
                        continue;
                    }
                    // ReSharper restore HeuristicUnreachableCode

                    if (scores == null)
                    {
                        ConsoleMsgUtils.ShowWarning("Null scores in WriteResultsToFile");
                        continue;
                    }

                    // Note for DblToString(value, 9, true), by having "9" and "true",
                    // values between 100 and 999 Da will have 7 digits after the decimal place, and
                    // values between 1000 and 9999 will have 6 digits after the decimal place

                    writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}",
                        scanNum,
                        match.Pre,
                        sequence, // Sequence
                        match.Post,
                        scores.Modifications,                                                 // Modifications
                        ion.Composition,                                                      // Composition
                        proteinName,                                                          // ProteinName
                        database.GetProteinDescription(match.Offset),                         // ProteinDescription
                        protLength,                                                           // ProteinLength
                        start,                                                                // Start
                        end,                                                                  // End
                        ion.Charge,                                                           // precursorCharge
                        StringUtilities.DblToString(ion.GetMostAbundantIsotopeMz(), 9, true), // MostAbundantIsotopeMz
                        StringUtilities.DblToString(ion.Composition.Mass, 9, true),           // Mass
                        match.Score,
                        scores.Score // Score (re-scored)
                    );
                }
            }
        }
    }
}
