using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;
using InformedProteomics.TopDown.Scoring;
using InformedProteomics.TopDown.TagBasedSearch;
using PNNLOmics.Utilities;

namespace InformedProteomics.TopDown.Execution
{
    public class IcTopDownLauncher
    {
        //public const int NumMatchesPerSpectrum = 1;
        public const string TargetFileNameEnding = "_IcTarget.tsv";
        public const string DecoyFileNameEnding = "_IcDecoy.tsv";
        public const string TdaFileNameEnding = "_IcTda.tsv";

        public IcTopDownLauncher(
            string specFilePath,
            string dbFilePath,
            string outputDir,
            AminoAcidSet aaSet,
            int minSequenceLength = 21,
            int maxSequenceLength = 500,
            int maxNumNTermCleavages = 1,
            int maxNumCTermCleavages = 0,
            int minPrecursorIonCharge = 2,
            int maxPrecursorIonCharge = 60,
            int minProductIonCharge = 1,
            int maxProductIonCharge = 20,
            double minSequenceMass = 2000.0,
            double maxSequenceMass = 50000.0,
            double precursorIonTolerancePpm = 10,
            double productIonTolerancePpm = 10,
            bool? runTargetDecoyAnalysis = true,
            int searchMode = 1,
            string featureFilePath = null,
            int maxThreads = 0,
            IEnumerable<int> scanNumbers = null,
            int numMatchesPerSpectrum = 3
            )
        {
            ErrorMessage = string.Empty;

            SpecFilePath = specFilePath;
            DatabaseFilePath = dbFilePath;
            AminoAcidSet = aaSet;
            FeatureFilePath = featureFilePath;

            OutputDir = outputDir;
            MinSequenceLength = minSequenceLength;
            MaxSequenceLength = maxSequenceLength;
            MaxNumNTermCleavages = maxNumNTermCleavages;
            MaxNumCTermCleavages = maxNumCTermCleavages;
            MinPrecursorIonCharge = minPrecursorIonCharge;
            MaxPrecursorIonCharge = maxPrecursorIonCharge;
            MinProductIonCharge = minProductIonCharge;
            MaxProductIonCharge = maxProductIonCharge;
            MinSequenceMass = minSequenceMass;
            MaxSequenceMass = maxSequenceMass;
            PrecursorIonTolerance = new Tolerance(precursorIonTolerancePpm);
            ProductIonTolerance = new Tolerance(productIonTolerancePpm);
            RunTargetDecoyAnalysis = runTargetDecoyAnalysis;
            SearchMode = searchMode;
            ScanNumbers = scanNumbers;
            NumMatchesPerSpectrum = numMatchesPerSpectrum;
            MaxNumThreads = maxThreads;
        }

        public string ErrorMessage { get; private set; }
        public string SpecFilePath { get; private set; }
        public string DatabaseFilePath { get; private set; }
        public string OutputDir { get; private set; }
        public AminoAcidSet AminoAcidSet { get; private set; }
        public string FeatureFilePath { get; private set; }
        public int MinSequenceLength { get; private set; }
        public int MaxSequenceLength { get; private set; }
        public int MaxNumNTermCleavages { get; private set; }
        public int MaxNumCTermCleavages { get; private set; }
        public int MinPrecursorIonCharge { get; private set; }
        public int MaxPrecursorIonCharge { get; private set; }
        public double MinSequenceMass { get; private set; }
        public double MaxSequenceMass { get; private set; }
        public int MinProductIonCharge { get; private set; }
        public int MaxProductIonCharge { get; private set; }
        public Tolerance PrecursorIonTolerance { get; private set; }
        public Tolerance ProductIonTolerance { get; private set; }
        public bool? RunTargetDecoyAnalysis { get; private set; }    // true: target and decoy, false: target only, null: decoy only
        public IEnumerable<int> ScanNumbers { get; private set; }
        public int NumMatchesPerSpectrum { get; private set; }
        public int MaxNumThreads { get; private set; }

        // 0: all internal sequences, 
        // 1: #NCleavages <= Max OR Cleavages <= Max (Default)
        // 2: 1: #NCleavages <= Max AND Cleavages <= Max
        public int SearchMode { get; private set; }

        private LcMsRun _run;
        //private ProductScorerBasedOnDeconvolutedSpectra _ms2ScorerFactory;
        private DeconvolutedSpectrumScorer _ms2ScorerFactory2;
        private IMassBinning _massBinComparer;
        private InformedTopDownScorer _topDownScorer;
        private ScanBasedTagSearchEngine _tagSearchEngine;

        public bool RunSearch(double corrThreshold = 0.7, CancellationToken? cancellationToken = null, IProgress<ProgressData> progress = null)
        {
            // Get the Normalized spec file/folder path
            SpecFilePath = MassSpecDataReaderFactory.NormalizeDatasetPath(SpecFilePath);

            var prog = new Progress<ProgressData>();
            var progData = new ProgressData(progress);
            if (progress != null)
            {
                prog = new Progress<ProgressData>(p =>
                {
                    progData.Status = p.Status;
                    progData.StatusInternal = p.StatusInternal;
                    progData.Report(p.Percent);
                });
            }

            var sw = new Stopwatch();
            var swAll = new Stopwatch();
            swAll.Start();
            ErrorMessage = string.Empty;

            Console.Write(@"Reading raw file...");
            progData.Status = "Reading spectra file";
            progData.StepRange(10.0);
            sw.Start();
            _run = PbfLcMsRun.GetLcMsRun(SpecFilePath, 0, 0, prog);
            _topDownScorer = new InformedTopDownScorer(_run, AminoAcidSet, MinProductIonCharge, MaxProductIonCharge, ProductIonTolerance, corrThreshold);
            sw.Stop();
            Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

            progData.StepRange(20.0);
            ISequenceFilter ms1Filter;
            if (string.IsNullOrWhiteSpace(FeatureFilePath))
            {
                // Checks whether SpecFileName.ms1ft exists
                var ms1FtFilePath = MassSpecDataReaderFactory.ChangeExtension(SpecFilePath, LcMsFeatureFinderLauncher.FileExtension);
                if (!File.Exists(ms1FtFilePath))
                {
                    Console.WriteLine(@"Running ProMex...");
                    sw.Start();
                    var param = new LcMsFeatureFinderInputParameter
                    {
                        InputPath = SpecFilePath,
                        MinSearchMass = MinSequenceMass,
                        MaxSearchMass = MaxSequenceMass,
                        MinSearchCharge = MinPrecursorIonCharge,
                        MaxSearchCharge = MaxPrecursorIonCharge,
                        CsvOutput = false,
                        ScoreReport = false,
                        LikelihoodScoreThreshold = -10
                    };
                    var featureFinder = new LcMsFeatureFinderLauncher(param);
                    featureFinder.Run();
                }
                sw.Reset();
                sw.Start();
                Console.Write(@"Reading ProMex results...");
                ms1Filter = new Ms1FtFilter(_run, PrecursorIonTolerance, ms1FtFilePath, -10);
            }
            else
            {
                sw.Reset();
                sw.Start();
                var extension = Path.GetExtension(FeatureFilePath);
                if (extension.ToLower().Equals(".csv"))
                {
                    Console.Write(@"Reading ICR2LS/Decon2LS results...");
                    ms1Filter = new IsosFilter(_run, PrecursorIonTolerance, FeatureFilePath);
                }
                else if (extension.ToLower().Equals(".ms1ft"))
                {
                    Console.Write(@"Reading ProMex results...");
                    ms1Filter = new Ms1FtFilter(_run, PrecursorIonTolerance, FeatureFilePath, -10);
                }
                else if (extension.ToLower().Equals(".msalign"))
                {
                    Console.Write(@"Reading MS-Align+ results...");
                    ms1Filter = new MsDeconvFilter(_run, PrecursorIonTolerance, FeatureFilePath);
                }
                else ms1Filter = null; //new Ms1FeatureMatrix(_run);
            }

            sw.Stop();
            Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

            // pre-generate deconvoluted spectra for scoring
            _massBinComparer = new FilteredProteinMassBinning(AminoAcidSet, MaxSequenceMass);
            _ms2ScorerFactory2 = new DeconvolutedSpectrumScorer(_run, _massBinComparer, AminoAcidSet,
                                                               MinProductIonCharge, MaxProductIonCharge, ProductIonTolerance);
            sw.Reset();
            Console.WriteLine(@"Generating deconvoluted spectra for MS/MS spectra...");
            sw.Start();
            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = MaxNumThreads,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };
            Parallel.ForEach(_run.GetScanNumbers(2), pfeOptions, ms2ScanNum =>
            {
                _ms2ScorerFactory2.DeconvonluteProductSpectrum(ms2ScanNum);
            });
            sw.Stop();
            Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

            progData.StepRange(10.0);
            progData.Status = "Reading Fasta File";

            // Target database
            var targetDb = new FastaDatabase(DatabaseFilePath);
            targetDb.Read();

            // Generate sequence tags for all MS/MS spectra
            progData.StepRange(25.0);
            progData.Status = "Generating Sequence Tags";

            sw.Reset();
            Console.WriteLine(@"Generating sequence tags for MS/MS spectra...");
            sw.Start();
            var seqTagGen = GetSequenceTagGenerator();
            _tagMs2ScanNum = seqTagGen.GetMs2ScanNumsContainingTags().ToArray();
            sw.Stop();
            Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

            _tagSearchEngine = new ScanBasedTagSearchEngine(_run, seqTagGen, new LcMsPeakMatrix(_run, ms1Filter), targetDb, ProductIonTolerance, AminoAcidSet,
                            _ms2ScorerFactory2,
                            ScanBasedTagSearchEngine.DefaultMinMatchedTagLength,
                            MaxSequenceMass, MinProductIonCharge, MaxProductIonCharge);

            var specFileName = MassSpecDataReaderFactory.RemoveExtension(Path.GetFileName(SpecFilePath));
            var targetOutputFilePath = Path.Combine(OutputDir, specFileName + TargetFileNameEnding);
            var decoyOutputFilePath = Path.Combine(OutputDir, specFileName + DecoyFileNameEnding);
            var tdaOutputFilePath = Path.Combine(OutputDir, specFileName + TdaFileNameEnding);

            progData.StepRange(60.0);
            progData.Status = "Running Target search";

            if (RunTargetDecoyAnalysis != null)
            {
                sw.Reset();
                Console.Write(@"Reading the target database...");
                sw.Start();
                targetDb.Read();
                sw.Stop();
                Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

                var targetMatches = new SortedSet<DatabaseSequenceSpectrumMatch>[_run.MaxLcScan + 1];

                progData.MaxPercentage = 42.5;

                sw.Reset();
                Console.WriteLine(@"Tag-based searching the target database");
                sw.Start();
                RunTagBasedSearch(targetMatches, targetDb, null, prog);
                Console.WriteLine(@"Target database tag-based search elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

                progData.MaxPercentage = 60.0;

                sw.Reset();
                Console.WriteLine(@"Searching the target database");
                sw.Start();
                RunSearch(targetMatches, targetDb, ms1Filter, null, prog);
                Console.WriteLine(@"Target database search elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

                // calculate spectral e-value usign generating function
                sw.Reset();
                Console.WriteLine(@"Calculating spectral E-values for target-spectrum matches");
                sw.Start();
                var bestTargetMatches = RunGeneratingFunction(targetMatches);
                WriteResultsToFile(bestTargetMatches, targetOutputFilePath, targetDb);
                sw.Stop();
                Console.WriteLine(@"Target-spectrum match E-value calculation elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);
            }

            progData.StepRange(95.0); // total to 95%
            progData.Status = "Running Decoy search";

            if (RunTargetDecoyAnalysis == true || RunTargetDecoyAnalysis == null)
            {
                // Decoy database
                sw.Reset();
                sw.Start();
                var decoyDb = targetDb.Decoy(null, true);

                Console.Write(@"Reading the decoy database...");
                decoyDb.Read();
                Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

                var decoyMatches = new SortedSet<DatabaseSequenceSpectrumMatch>[_run.MaxLcScan + 1];

                progData.MaxPercentage = 77.5;

                sw.Reset();
                Console.WriteLine(@"Tag-based searching the decoy database");
                sw.Start();
                RunTagBasedSearch(decoyMatches, decoyDb, null, prog);
                Console.WriteLine(@"Decoy database tag-based search elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

                progData.MaxPercentage = 95.0;

                sw.Reset();
                Console.WriteLine(@"Searching the decoy database");
                sw.Start();
                RunSearch(decoyMatches, decoyDb, ms1Filter, null, prog);
                Console.WriteLine(@"Decoy database search elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

                // calculate spectral e-value usign generating function
                sw.Reset();
                Console.WriteLine(@"Calculating spectral E-values for decoy-spectrum matches");
                sw.Start();
                var bestDecoyMatches = RunGeneratingFunction(decoyMatches);
                WriteResultsToFile(bestDecoyMatches, decoyOutputFilePath, decoyDb);
                sw.Stop();
                Console.WriteLine(@"Decoy-spectrum match E-value calculation elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);
            }

            progData.StepRange(100.0);
            progData.Status = "Writing combined results file";
            if (RunTargetDecoyAnalysis == true)
            {
                // Add "Qvalue" and "PepQValue"
                var fdrCalculator = new FdrCalculator(targetOutputFilePath, decoyOutputFilePath);
                if (fdrCalculator.HasError())
                {
                    ErrorMessage = fdrCalculator.ErrorMessage;
                    Console.WriteLine(@"Error computing FDR: " + fdrCalculator.ErrorMessage);
                    return false;
                }

                fdrCalculator.WriteTo(tdaOutputFilePath);
            }
            progData.Report(100.0);

            Console.WriteLine(@"Done.");
            swAll.Stop();
            Console.WriteLine(@"Total elapsed time for search: {0:f1} sec ({1:f2} min)", swAll.Elapsed.TotalSeconds, swAll.Elapsed.TotalMinutes);

            return true;
        }

        private int[] _tagMs2ScanNum;

        private IEnumerable<AnnotationAndOffset> GetAnnotationsAndOffsets(FastaDatabase database, out long estimatedProteins, CancellationToken? cancellationToken = null)
        {
            var indexedDb = new IndexedDatabase(database);
            indexedDb.Read();
            estimatedProteins = indexedDb.EstimateTotalPeptides(SearchMode, MinSequenceLength, MaxSequenceLength, MaxNumNTermCleavages, MaxNumCTermCleavages);
            IEnumerable<AnnotationAndOffset> annotationsAndOffsets;
            if (SearchMode == 0)
            {
                //annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(MinSequenceLength, MaxSequenceLength);
                annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzymeParallel(MinSequenceLength, MaxSequenceLength, MaxNumThreads, cancellationToken);
            }
            else if (SearchMode == 2)
            {
                annotationsAndOffsets = indexedDb.IntactSequenceAnnotationsAndOffsets(MinSequenceLength, MaxSequenceLength, MaxNumCTermCleavages);
            }
            else
            {
                annotationsAndOffsets = indexedDb
                    .SequenceAnnotationsAndOffsetsWithNtermOrCtermCleavageNoLargerThan(
                        MinSequenceLength, MaxSequenceLength, MaxNumNTermCleavages, MaxNumCTermCleavages);
            }

            return annotationsAndOffsets;
        }


        private void RunTagBasedSearch(SortedSet<DatabaseSequenceSpectrumMatch>[] matches, FastaDatabase db,
                                        CancellationToken? cancellationToken = null, IProgress<ProgressData> progress = null)
        {
            _tagSearchEngine.SetDatabase(db);
            //var ms2ScanNums = _run.GetScanNumbers(2);
            var progData = new ProgressData(progress)
            {
                Status = "Tag-based Searching for matches"
            };

            var sw = new Stopwatch();

            long estimatedProteins = _tagMs2ScanNum.Length;
            Console.WriteLine(@"Number of spectra containing sequence tags: " + estimatedProteins);
            var numProteins = 0;
            var lastUpdate = DateTime.MinValue; // Force original update of 0%

            sw.Reset();
            sw.Start();

            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = MaxNumThreads,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            Parallel.ForEach(_tagMs2ScanNum, pfeOptions, ms2ScanNum =>
            {
                var tagSeqMatches = _tagSearchEngine.RunSearch(ms2ScanNum);
                //var prsmList = new List<DatabaseSequenceSpectrumMatch>();
                foreach (var tagSequenceMatch in tagSeqMatches)
                {
                    var offset = _tagSearchEngine.FastaDatabase.GetOffset(tagSequenceMatch.ProteinName);
                    if (offset == null) continue;

                    var sequence = tagSequenceMatch.Sequence;
                    var numNTermCleavages = tagSequenceMatch.TagMatch.StartIndex;

                    var seqObj = Sequence.CreateSequence(sequence, tagSequenceMatch.TagMatch.Modifications, AminoAcidSet);
                    var precursorIon = new Ion(seqObj.Composition + Composition.H2O, tagSequenceMatch.TagMatch.Charge);

                    var prsm = new DatabaseSequenceSpectrumMatch(sequence, tagSequenceMatch.Pre, tagSequenceMatch.Post,
                                                                 ms2ScanNum, (long)offset, numNTermCleavages,
                                                                 null,
                                                                 precursorIon, tagSequenceMatch.TagMatch.Score)
                    {
                        ModificationText = tagSequenceMatch.TagMatch.Modifications,
                    };

                    AddMatch(matches, ms2ScanNum, prsm);    
                }

                SearchProgressReport(ref numProteins, ref lastUpdate, estimatedProteins, sw, progData, "spectra");
            });

            Console.WriteLine(@"Collected candidate matches: {0}", GetNumberOfMatches(matches));

            progData.StatusInternal = string.Empty;
            progData.Report(100.0);
        }

        private void RunSearch(SortedSet<DatabaseSequenceSpectrumMatch>[] matches, FastaDatabase db, ISequenceFilter sequenceFilter, CancellationToken? cancellationToken = null, IProgress<ProgressData> progress = null)
        {
            var progData = new ProgressData(progress)
            {
                Status = "Searching for matches"
            };

            var sw = new Stopwatch();
            long estimatedProteins;
            var annotationsAndOffsets = GetAnnotationsAndOffsets(db, out estimatedProteins, cancellationToken);
            Console.WriteLine(@"Estimated proteins: " + estimatedProteins);

            var numProteins = 0;
            var lastUpdate = DateTime.MinValue; // Force original update of 0%

            sw.Reset();
            sw.Start();
            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = MaxNumThreads,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            var maxNumNTermCleavages = SearchMode == 2 ? MaxNumNTermCleavages : 0;
            //foreach (var annotationAndOffset in annotationsAndOffsets)
            Parallel.ForEach(annotationsAndOffsets, pfeOptions, annotationAndOffset =>
            {
                if (cancellationToken != null && cancellationToken.Value.IsCancellationRequested)
                {
                    //return matches;
                    return;
                }
                SearchProgressReport(ref numProteins, ref lastUpdate, estimatedProteins, sw, progData);
                SearchForMatches(annotationAndOffset, sequenceFilter, matches, maxNumNTermCleavages, cancellationToken);
            });

            Console.WriteLine(@"Collected candidate matches: {0}", GetNumberOfMatches(matches));

            progData.StatusInternal = string.Empty;
            progData.Report(100.0);
        }

        private void SearchProgressReport(ref int numProteins, ref DateTime lastUpdate, long estimatedProteins, Stopwatch sw, ProgressData progData, string itemName = "proteins")
        {
            var tempNumProteins = Interlocked.Increment(ref numProteins) - 1;

            if (estimatedProteins < 1)
                estimatedProteins = 1;

            progData.StatusInternal = string.Format(@"Processing, {0} {1} done, {2:#0.0}% complete, {3:f1} sec elapsed",
                    tempNumProteins,
                    itemName,
                    tempNumProteins / (double)estimatedProteins * 100.0,
                    sw.Elapsed.TotalSeconds);
            progData.Report(tempNumProteins, estimatedProteins);

            int secondsThreshold;

            if (sw.Elapsed.TotalMinutes < 2)
                secondsThreshold = 15;      // Every 15 seconds
            else if (sw.Elapsed.TotalMinutes < 5)
                secondsThreshold = 30;      // Every 30 seconds
            else if (sw.Elapsed.TotalMinutes < 20)
                secondsThreshold = 60;      // Every 1 minute
            else
                secondsThreshold = 300;     // Every 5 minutes

            if (DateTime.UtcNow.Subtract(lastUpdate).TotalSeconds >= secondsThreshold)
            {
                lastUpdate = DateTime.UtcNow;

                Console.WriteLine(@"Processing, {0} {1} done, {2:#0.0}% complete, {3:f1} sec elapsed",
                    tempNumProteins,
                    itemName,
                    tempNumProteins / (double)estimatedProteins * 100.0,
                    sw.Elapsed.TotalSeconds);
            }
        }

        private const int ScoreLowerBound = 3;
        private void SearchForMatches(AnnotationAndOffset annotationAndOffset,
            ISequenceFilter sequenceFilter, SortedSet<DatabaseSequenceSpectrumMatch>[] matches, int maxNumNTermCleavages, CancellationToken? cancellationToken = null)
        {
            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = MaxNumThreads,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            var annotation = annotationAndOffset.Annotation;
            var offset = annotationAndOffset.Offset;
            //var protein = db.GetProteinName(offset);
            var protSequence = annotation.Substring(2, annotation.Length - 4);
            var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, AminoAcid.ProteinNTerm, protSequence,
                AminoAcid.ProteinCTerm);

            if (seqGraph == null) return; // No matches will be found without a sequence graph.

            for (var numNTermCleavages = 0; numNTermCleavages <= maxNumNTermCleavages; numNTermCleavages++)
            {
                if (numNTermCleavages > 0) seqGraph.CleaveNTerm();
                var numProteoforms = seqGraph.GetNumProteoforms();
                var modCombs = seqGraph.GetModificationCombinations();
                for (var modIndex = 0; modIndex < numProteoforms; modIndex++)
                {
                    seqGraph.SetSink(modIndex);
                    var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();
                    var sequenceMass = protCompositionWithH2O.Mass;

                    if (sequenceMass < MinSequenceMass || sequenceMass > MaxSequenceMass) continue;

                    var modCombinations = modCombs[modIndex];
                    var ms2ScanNums = this.ScanNumbers ?? sequenceFilter.GetMatchingMs2ScanNums(sequenceMass);

                    Parallel.ForEach(ms2ScanNums, pfeOptions, ms2ScanNum =>
                    {
                        if (ms2ScanNum > _run.MaxLcScan) return;

                        var scorer = _ms2ScorerFactory2.GetMs2Scorer(ms2ScanNum);
                        var score = seqGraph.GetFragmentScore(scorer);
                        if (score <= ScoreLowerBound) return;

                        var spec = _run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                        if (spec == null) return;
                        var charge = (int)Math.Round(sequenceMass / (spec.IsolationWindow.IsolationWindowTargetMz - Constants.Proton));

                        var precursorIon = new Ion(protCompositionWithH2O, charge);
                        var sequence = protSequence.Substring(numNTermCleavages);
                        var pre = numNTermCleavages == 0 ? annotation[0] : annotation[numNTermCleavages + 1];
                        var post = annotation[annotation.Length - 1];

                        var prsm = new DatabaseSequenceSpectrumMatch(sequence, pre, post, ms2ScanNum, offset,
                            numNTermCleavages,
                            modCombinations, precursorIon, score);

                        AddMatch(matches, ms2ScanNum, prsm);
                    });
                }
            }
        }

        private void AddMatch(SortedSet<DatabaseSequenceSpectrumMatch>[] matches, int ms2ScanNum, DatabaseSequenceSpectrumMatch prsm)
        {
            lock (matches)
            {
                if (matches[ms2ScanNum] == null)
                {
                    matches[ms2ScanNum] = new SortedSet<DatabaseSequenceSpectrumMatch> {prsm};
                }
                else // already exists
                {
                    var existingMatches = matches[ms2ScanNum];
                    var maxScore = existingMatches.Max.Score;
                    if (existingMatches.Count < NumMatchesPerSpectrum && maxScore * 0.6 < prsm.Score)
                    {
                        existingMatches.Add(prsm);
                        existingMatches.RemoveWhere(mt => mt.Score < maxScore * 0.6);
                    }
                    else
                    {
                        var minScore = existingMatches.Min.Score;
                        if (prsm.Score > minScore)
                        {
                            existingMatches.Add(prsm);
                            existingMatches.RemoveWhere(mt => mt.Score < maxScore * 0.6);
                        }
                    }                              
                }
            }
        }

        private SequenceTagGenerator GetSequenceTagGenerator(CancellationToken? cancellationToken = null, IProgress<ProgressData> progress = null)
        {
            var sequenceTagGen = new SequenceTagGenerator(_run, new Tolerance(5));
            var scanNums = _run.GetScanNumbers(2);

            var progData = new ProgressData(progress)
            {
                Status = "Generating sequence tags"
            };

            var sw = new Stopwatch();

            // Rescore and Estimate #proteins for GF calculation
            long estimatedProteins = scanNums.Count;
            Console.WriteLine(@"Number of spectra: " + estimatedProteins);
            var numProteins = 0;
            var lastUpdate = DateTime.MinValue; // Force original update of 0%
            sw.Reset();
            sw.Start();

            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = MaxNumThreads,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            Parallel.ForEach(scanNums, pfeOptions, scanNum =>
            {
                sequenceTagGen.Generate(scanNum);
                SearchProgressReport(ref numProteins, ref lastUpdate, estimatedProteins, sw, progData,
                                     "spectra");
            });

            progData.StatusInternal = string.Empty;
            progData.Report(100.0);
            Console.WriteLine(@"Generated sequence tags: " + sequenceTagGen.NumberOfGeneratedTags());
            return sequenceTagGen;
        }

        private DatabaseSequenceSpectrumMatch[] RunGeneratingFunction(SortedSet<DatabaseSequenceSpectrumMatch>[] matches, CancellationToken? cancellationToken = null, IProgress<ProgressData> progress = null)
        {
            //const double massToleranceForGf = 10000;
            var progData = new ProgressData(progress)
            {
                Status = "Calculating spectral E-values for matches"
            };

            var sw = new Stopwatch();

            // Rescore and Estimate #proteins for GF calculation
            long estimatedProteins = 0;
            for (var scanNum = _run.MinLcScan; scanNum <= _run.MaxLcScan; scanNum++)
            {
                if (matches[scanNum] == null) continue;

                var highestScore = 0d;
                foreach (var match in matches[scanNum])
                {
                    var sequence = match.Sequence;
                    var ion = match.Ion;

                    // Re-scoring
                    var scores = _topDownScorer.GetScores(AminoAcid.ProteinNTerm, sequence, AminoAcid.ProteinCTerm, ion.Composition, ion.Charge, scanNum);
                    if (scores == null) continue;
                    match.Score = scores.Ms2Score;
                    match.ModificationText = scores.Modifications;
                    highestScore = Math.Max(highestScore, scores.Ms2Score);
                }

                matches[scanNum].RemoveWhere(m => m.Score <= ScoreLowerBound || m.Score < highestScore * 0.7);
                estimatedProteins += matches[scanNum].Count;
            }

            Console.WriteLine(@"Estimated matched proteins: " + estimatedProteins);

            var numProteins = 0;
            var lastUpdate = DateTime.MinValue; // Force original update of 0%
            sw.Reset();
            sw.Start();
            var finalMatches = new DatabaseSequenceSpectrumMatch[matches.Length];

            var scanNums = new List<int>();
            for (var scanNum = _run.MinLcScan; scanNum <= _run.MaxLcScan; scanNum++)
                if (matches[scanNum] != null) scanNums.Add(scanNum);

            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = MaxNumThreads,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            Parallel.ForEach(scanNums, pfeOptions, scanNum =>
            {
                var currentTask = "?";
                try
                {
                    currentTask = "Initializing a new GeneratingFunction with NumberOfBins = " +
                                  _massBinComparer.NumberOfBins;

                    var gf = new GeneratingFunction(_massBinComparer.NumberOfBins);

                    foreach (var match in matches[scanNum])
                    {
                        var currentIteration = "for scan " + scanNum + " and mass " + match.Ion.Composition.Mass;
                        currentTask = "Calling GetMs2ScoringGraph " + currentIteration;
                        
                        var graph = _ms2ScorerFactory2.GetMs2ScoringGraph(scanNum, match.Ion.Composition.Mass);
                        if (graph == null) continue;

                        currentTask = "Calling ComputeGeneratingFunction " + currentIteration;
                        gf.ComputeGeneratingFunction(graph);

                        currentTask = "Calling GetSpectralEValue " + currentIteration + " and score " + (int)match.Score;
                        match.SpecEvalue = gf.GetSpectralEValue((int)match.Score);

                        currentTask = "Locking finalMatches " + currentIteration;
                        lock (finalMatches)
                        {
                            currentTask = "Comparing new SpecEvalue to stored SpecEvalue " + currentIteration;
                            if (finalMatches[scanNum] == null || finalMatches[scanNum].SpecEvalue > match.SpecEvalue)
                            {
                                currentTask = "Updating stored SpecEvalue " + currentIteration;
                                finalMatches[scanNum] = match;
                            }
                        }

                        currentTask = "Reporting progress " + currentIteration;
                        SearchProgressReport(ref numProteins, ref lastUpdate, estimatedProteins, sw, progData);
                    }
                }
                catch (Exception ex)
                {
                    var errMsg = string.Format("Exception while {0}: {1}", currentTask, ex.Message);
                    Console.WriteLine(errMsg);
                    throw new Exception(errMsg, ex);
                }
            });

            progData.StatusInternal = string.Empty;
            progData.Report(100.0);
            return finalMatches;
        }

        private int GetNumberOfMatches(SortedSet<DatabaseSequenceSpectrumMatch>[] matches)
        {
            var nMatches = 0;
            lock (matches)
            {
                for (var scanNum = _run.MinLcScan; scanNum <= _run.MaxLcScan; scanNum++)
                {
                    if (matches[scanNum] == null) continue;
                    nMatches += matches[scanNum].Count;
                }
            }
            return nMatches;
        }

        private double ExcelMinValue(double val)
        {
            const double minExcelValue = 9.99E-308;
            return Math.Max(val, minExcelValue);
        }

        private void WriteResultsToFile(DatabaseSequenceSpectrumMatch[] matches, string outputFilePath, FastaDatabase database)
        {
            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine("Scan\tPre\tSequence\tPost\tModifications\tComposition\tProteinName\tProteinDesc" +
                                 "\tProteinLength\tStart\tEnd\tCharge\tMostAbundantIsotopeMz\tMass\t#MatchedFragments\tSpecEValue\tEValue"
                    );
                for (var scanNum = _run.MinLcScan; scanNum <= _run.MaxLcScan; scanNum++)
                {
                    var match = matches[scanNum];
                    if (match == null)
                        continue;

                    var sequence = match.Sequence;
                    var offset = match.Offset;
                    var start = database.GetOneBasedPositionInProtein(offset) + 1 + match.NumNTermCleavages;
                    var end = start + sequence.Length - 1;
                    var proteinName = database.GetProteinName(match.Offset);
                    var protLength = database.GetProteinLength(proteinName);
                    var ion = match.Ion;

                    var proteinDescription = database.GetProteinDescription(match.Offset);

                    //var scores = _topDownScorer.GetScores(AminoAcid.ProteinNTerm, sequence, AminoAcid.ProteinCTerm, ion.Composition, ion.Charge, scanNum);

                    // Note for DblToString(value, 9, true), by having "9" and "true",
                    // values between 100 and 999 Da will have 7 digits after the decimal place, and
                    // values between 1000 and 9999 will have 6 digits after the decimal place
                    writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}",
                        scanNum,
                        match.Pre,                 // Pre
                        sequence,                  // Sequence
                        match.Post,                // Post
                        match.ModificationText,    // Modifications
                        ion.Composition,           // Composition
                        proteinName,               // ProteinName
                        proteinDescription,        // ProteinDescription
                        protLength,                // ProteinLength
                        start,                     // Start position in protein
                        end,                       // End position in protein
                        ion.Charge,                // precursorCharge
                        StringUtilities.DblToString(ion.GetMostAbundantIsotopeMz(), 9, true), // MostAbundantIsotopeMz
                        StringUtilities.DblToString(ion.Composition.Mass, 9, true),           // Mass
                        StringUtilities.DblToString(match.Score, 4),                          // Score (Number of matched fragments)
                        StringUtilities.DblToString(ExcelMinValue(match.SpecEvalue), 6, true, 0.001),                             // EValue; will be displayed using scientific notation if the value is less than 0.001
                        StringUtilities.DblToString(ExcelMinValue(match.SpecEvalue * database.GetNumEntries()), 6, true, 0.001)   // SpecEValue; will be displayed using scientific notation if the value is less than 0.001
                        );

                }
            }
        }
    }
}
