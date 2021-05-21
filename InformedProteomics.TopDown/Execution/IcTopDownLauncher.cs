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
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.SearchResults;
using InformedProteomics.FeatureFinding;
using InformedProteomics.FeatureFinding.FeatureDetection;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;
using InformedProteomics.TopDown.Scoring;
using InformedProteomics.TopDown.TagBasedSearch;
using PRISM;

namespace InformedProteomics.TopDown.Execution
{
    public class IcTopDownLauncher : EventNotifier
    {
        // Ignore Spelling: deconvolute, deconvoluted, ic, Nums, pbf, Pre

        private const bool USE_PARALLEL_FOREACH = true;

        private const bool DEBUG_MODE = false;
        private const int DEBUG_MODE_PROTEINS_TO_SEARCH = 275;

        public const string TargetFileNameEnding = "_IcTarget.tsv";
        public const string DecoyFileNameEnding = "_IcDecoy.tsv";
        public const string TdaFileNameEnding = "_IcTda.tsv";
        public const string MzidFileNameEnding = ".mzid";

        [Obsolete("Use IcTopDownLauncher(MsPfExecutionOptions)!", true)]
        public IcTopDownLauncher(
            string specFilePath,
            string dbFilePath,
            string outputDir,
            AminoAcidSet aaSet,
            string featureFilePath = null)
        {
            ErrorMessage = string.Empty;

            Options = new MsPfParameters(specFilePath, dbFilePath, outputDir, aaSet, featureFilePath);
        }

        public IcTopDownLauncher(MsPfParameters options)
        {
            ErrorMessage = string.Empty;

            Options = options;
        }

        public string ErrorMessage { get; private set; }
        public MsPfParameters Options { get; }

        private LcMsRun _run;
        private CompositeScorerFactory _ms2ScorerFactory2;

        //private ScoringGraphFactory scoringGraphFactory;
        private IMassBinning _massBinComparer;
        private ScanBasedTagSearchEngine _tagSearchEngine;
        private double[] _isolationWindowTargetMz; // spec.IsolationWindow.IsolationWindowTargetMz
        private List<int> _ms2ScanNums;
        private ProgressData searchProgressData;
        private Stopwatch searchStopwatch;
        private Timer progReportTimer;

        private void ReportOverallProgress(object searchObj)
        {
            if (!(searchObj is IcTopDownLauncher searchClass))
            {
                return;
            }

            var progData = searchClass.searchProgressData;
            var timeElapsed = searchClass.searchStopwatch.Elapsed;

            TimeSpan estimatedRemaining;
            if (progData.Percent > 5)
            {
                var estimatedTotal = TimeSpan.FromSeconds(timeElapsed.TotalSeconds / progData.Percent * 100.0);
                estimatedRemaining = estimatedTotal.Subtract(timeElapsed);
            }
            else
            {
                estimatedRemaining = TimeSpan.FromSeconds(0);
            }

            var progMsg = string.Format("Total Progress: {0:F2}%, {1:%d}d {1:%h}h {1:%m}.{2:00}m elapsed, Current Task: {3}, estimated remaining: {4:%d}d {4:%h}h {4:%m}.{5:00}m",
                progData.Percent, timeElapsed, timeElapsed.Seconds / 60.0 * 100, progData.Status, estimatedRemaining, estimatedRemaining.Seconds / 60.0 * 100);

            OnProgressUpdate(progMsg, (float)progData.Percent);
        }

        public bool RunSearch(CancellationToken? cancellationToken = null, IProgress<ProgressData> progress = null)
        {
            // Get the Normalized spec file/folder path
            Options.SpecFilePath = MassSpecDataReaderFactory.NormalizeDatasetPath(Options.SpecFilePath);

            var progData = new ProgressData(progress);
            searchProgressData = progData;
            var prog = new Progress<ProgressData>(p =>
            {
                progData.Status = p.Status;
                progData.StatusInternal = p.StatusInternal;
                progData.Report(p.Percent);
            });

            var sw = new Stopwatch();
            var swAll = new Stopwatch();
            searchStopwatch = swAll;
            swAll.Start();
            ErrorMessage = string.Empty;

            // Output a progress message every 5 minutes...
            progReportTimer = new Timer(ReportOverallProgress, this, 0, 1000 * 60 * 5);

            var specFileName = Path.GetFileName(Options.SpecFilePath);

            if (string.Equals(Path.GetExtension(Options.SpecFilePath), ".pbf", StringComparison.InvariantCultureIgnoreCase))
            {
                UpdateStatus("Reading pbf file " + specFileName, progData);
                progData.StepRange(1.0, "Reading spectra file");
            }
            else
            {
                UpdateStatus("Creating and/or reading pbf file for " + specFileName, progData);
                progData.StepRange(5.0, "Reading spectra file");
            }

            sw.Start();

            _run = PbfLcMsRun.GetLcMsRun(Options.SpecFilePath, 0, 0, prog);

            var minMs1Scan = int.MinValue;
            var maxMs1Scan = int.MaxValue;

            // Retrieve the list of MS2 scans
            _ms2ScanNums = _run.GetScanNumbers(2).ToList();

            if (Options.ScanNumbers?.Any() == true)
            {
                // Filter the MS2 scans using ScanNumbers
                _ms2ScanNums = _ms2ScanNums.Intersect(Options.ScanNumbers).ToList();

                minMs1Scan = _ms2ScanNums.Min() - 100;
                maxMs1Scan = _ms2ScanNums.Max() + 100;
            }

            _isolationWindowTargetMz = new double[_run.MaxLcScan + 1];
            foreach (var ms2Scan in _ms2ScanNums)
            {
                if (!(_run.GetSpectrum(ms2Scan) is ProductSpectrum ms2Spec))
                {
                    continue;
                }

                _isolationWindowTargetMz[ms2Scan] = ms2Spec.IsolationWindow.IsolationWindowTargetMz;
            }

            sw.Stop();
            OnStatusEvent(string.Format("Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds));

            progData.StepRange(progData.MaxPercentage + 2, "Reading Fasta File");
            OnStatusEvent(progData.Status);

            // Target database
            var targetDb = new FastaDatabase(Options.DatabaseFilePath);
            targetDb.Read();

            ISequenceFilter ms1Filter;
            if (Options.ScanNumbers?.Any() == true)
            {
                OnStatusEvent("Filtering MS1 data using the specified scan numbers");
                ms1Filter = new SelectedMsMsFilter(Options.ScanNumbers);
            }
            else if (string.IsNullOrWhiteSpace(Options.FeatureFilePath))
            {
                // Checks whether SpecFileName.ms1ft exists
                var ms1FtFilePath = MassSpecDataReaderFactory.ChangeExtension(Options.SpecFilePath, LcMsFeatureFinderLauncher.FileExtension);
                if (!File.Exists(ms1FtFilePath))
                {
                    progData.StepRange(progData.MaxPercentage + 10);
                    UpdateStatus("Running ProMex...", progData);
                    sw.Reset();
                    sw.Start();
                    var param = new LcMsFeatureFinderInputParameters
                    {
                        InputPath = Options.SpecFilePath,
                        MinSearchMass = Options.MinSequenceMass,
                        MaxSearchMass = Options.MaxSequenceMass,
                        MinSearchCharge = Options.MinPrecursorIonCharge,
                        MaxSearchCharge = Options.MaxPrecursorIonCharge,
                        CsvOutput = false,
                        ScoreReport = false,
                        LikelihoodScoreThreshold = -10
                    };
                    var featureFinder = new LcMsFeatureFinderLauncher(param);
                    featureFinder.Run();
                }
                sw.Reset();
                sw.Start();
                OnStatusEvent("Reading ProMex results from " + Path.GetFileName(ms1FtFilePath));
                ms1Filter = new Ms1FtFilter(_run, Options.PrecursorIonTolerance, ms1FtFilePath, -10);
            }
            else
            {
                var featureFileName = Path.GetFileName(Options.FeatureFilePath);

                progData.StepRange(progData.MaxPercentage + 1);
                sw.Reset();
                sw.Start();
                var extension = Path.GetExtension(Options.FeatureFilePath);
                if (extension.Equals(".csv", StringComparison.OrdinalIgnoreCase))
                {
                    OnStatusEvent("Reading ICR2LS/Decon2LS results from " + featureFileName);
                    ms1Filter = new IsosFilter(_run, Options.PrecursorIonTolerance, Options.FeatureFilePath);
                }
                else if (extension.Equals(".ms1ft", StringComparison.OrdinalIgnoreCase))
                {
                    OnStatusEvent("Reading ProMex results from " + featureFileName);
                    ms1Filter = new Ms1FtFilter(_run, Options.PrecursorIonTolerance, Options.FeatureFilePath, -10);
                }
                else if (extension.Equals(".msalign", StringComparison.OrdinalIgnoreCase))
                {
                    OnStatusEvent("Reading MS-Align+ results from " + featureFileName);
                    ms1Filter = new MsDeconvFilter(_run, Options.PrecursorIonTolerance, Options.FeatureFilePath);
                }
                else
                {
                    OnStatusEvent("Note: MS1 feature filter is not defined");
                    ms1Filter = null; //new Ms1FeatureMatrix(_run);
                }
            }

            sw.Stop();
            OnStatusEvent(string.Format("Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds));

            // Pre-generate deconvoluted spectra for scoring
            _massBinComparer = new FilteredProteinMassBinning(Options.AminoAcidSet, Options.MaxSequenceMass + 1000);

            //this.fragmentScorerFactory = new CompositionScorerFactory(_run, true);

            sw.Reset();

            UpdateStatus("Generating deconvoluted spectra for MS/MS spectra...", progData);
            sw.Start();
            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Options.MaxDegreeOfParallelism,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            // Deconvolute spectra
            if (InformedProteomics.Backend.Utils.FlipSwitch.UseFlipScoring)
            {
                var deconvoluter = new Deconvoluter(Options.MinProductIonCharge, Options.MaxProductIonCharge, 2, Options.ProductIonTolerance);
                var lcmsRunDeconvoluter = new LcmsRunDeconvoluter(_run, deconvoluter, 2, pfeOptions.MaxDegreeOfParallelism);
                if (!File.Exists(MassSpecDataReaderFactory.ChangeExtension(Options.SpecFilePath, DPbfLcMsRun.FileExtensionConst)))
                {
                    progData.StepRange(progData.MaxPercentage + 2);
                }
                var deconvolutedRun = new DPbfLcMsRun(Options.SpecFilePath, lcmsRunDeconvoluter, keepDataReaderOpen: true);

                _ms2ScorerFactory2 = new CompositeScorerFactory(deconvolutedRun, _massBinComparer, Options.AminoAcidSet,
                                                       Options.MinProductIonCharge, Options.MaxProductIonCharge, Options.ProductIonTolerance, fullRun: _run as PbfLcMsRun);
            }
            else
            {
                _ms2ScorerFactory2 = new CompositeScorerFactory(_run, _massBinComparer, Options.AminoAcidSet,
                                                                Options.MinProductIonCharge, Options.MaxProductIonCharge, Options.ProductIonTolerance);

#pragma warning disable 162
                if (USE_PARALLEL_FOREACH)
                {
                    Parallel.ForEach(_ms2ScanNums, pfeOptions, ms2ScanNum =>
                    {
                        //_ms2ScorerFactory2.GetScorer(ms2ScanNum, activationMethod: Options.ActivationMethod);
                        _ms2ScorerFactory2.DeconvoluteProductSpectrum(ms2ScanNum);
                    });
                }
                else
                {
                    foreach (var ms2ScanNum in _ms2ScanNums)
                    {
                        _ms2ScorerFactory2.DeconvoluteProductSpectrum(ms2ScanNum);
                    }
                }
#pragma warning restore 162
            }

            sw.Stop();
            OnStatusEvent(string.Format("Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds));

            // Generate sequence tags for all MS/MS spectra
            if (Options.TagBasedSearch)
            {
                progData.StepRange(progData.MaxPercentage + 2, "Generating Sequence Tags");

                UpdateStatus("Generating sequence tags for MS/MS spectra...", progData);

                sw.Reset();
                sw.Start();
                var seqTagGen = GetSequenceTagGenerator();
                _tagMs2ScanNum = seqTagGen.GetMs2ScanNumsContainingTags().ToArray();
                sw.Stop();
                OnStatusEvent(string.Format("Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds));

                var peakMatrix = new LcMsPeakMatrix(_run, ms1Filter, minMs1Scan: minMs1Scan, maxMs1Scan: maxMs1Scan);
                _tagSearchEngine = new ScanBasedTagSearchEngine(
                    _run, seqTagGen, peakMatrix, targetDb, Options.ProductIonTolerance, Options.AminoAcidSet,
                    _ms2ScorerFactory2, ScanBasedTagSearchEngine.DefaultMinMatchedTagLength,
                    Options.MaxSequenceMass, Options.MinProductIonCharge, Options.MaxProductIonCharge);
            }

            var datasetName = MassSpecDataReaderFactory.RemoveExtension(Path.GetFileName(Options.SpecFilePath));
            var targetOutputFilePath = Path.Combine(Options.OutputDir, datasetName + TargetFileNameEnding);
            var decoyOutputFilePath = Path.Combine(Options.OutputDir, datasetName + DecoyFileNameEnding);
            var tdaOutputFilePath = Path.Combine(Options.OutputDir, datasetName + TdaFileNameEnding);
            var mzidOutputFilePath = Path.Combine(Options.OutputDir, datasetName + MzidFileNameEnding);

            progData.StepRange(progData.MaxPercentage + (98 - progData.MaxPercentage) / 2.0, "Running Target search");
            List<DatabaseSearchResultData> targetSearchResults = null;

            var validTargetResults = ResultsFileHasData(targetOutputFilePath);
            if (Options.TargetDecoySearchMode.HasFlag(DatabaseSearchMode.Target) && !validTargetResults || Options.OverwriteExistingResults)
            {
                targetSearchResults = RunDatabaseSearch(targetDb, targetOutputFilePath, ms1Filter, "target", prog);
            }
            else if (validTargetResults)
            {
                OnWarningEvent(string.Format("Target results file '{0}' exists; skipping target search.", targetOutputFilePath));
            }

            progData.StepRange(98.0, "Running Decoy search"); // total to 98%
            List<DatabaseSearchResultData> decoySearchResults = null;

            var validDecoyResults = ResultsFileHasData(decoyOutputFilePath);
            if (Options.TargetDecoySearchMode.HasFlag(DatabaseSearchMode.Decoy) && !validDecoyResults || Options.OverwriteExistingResults)
            {
                var decoyDb = targetDb.Decoy(null, true);
                decoySearchResults = RunDatabaseSearch(decoyDb, decoyOutputFilePath, ms1Filter, "decoy", prog);
            }
            else if (validDecoyResults)
            {
                OnWarningEvent(string.Format("Decoy results file  '{0}' exists; skipping decoy search.", decoyOutputFilePath));
            }

            progData.StepRange(100.0, "Writing combined results file");
            if (Options.TargetDecoySearchMode.HasFlag(DatabaseSearchMode.Both))
            {
                var multiplePeptidesPerScan = Options.MatchesPerSpectrumToReport > 1;

                // Add "QValue" and "PepQValue"
                FdrCalculator fdrCalculator;
                if (targetSearchResults == null && decoySearchResults == null)
                {
                    // Objects not populated, try to read in the files.
                    fdrCalculator = new FdrCalculator(targetOutputFilePath, decoyOutputFilePath, multiplePeptidesPerScan);
                }
                else if (targetSearchResults == null)
                {
                    // Target search skipped, read in the result file
                    targetSearchResults = DatabaseSearchResultData.ReadResultsFromFile(targetOutputFilePath);
                    fdrCalculator = new FdrCalculator(targetSearchResults, decoySearchResults, multiplePeptidesPerScan);
                }
                else if (decoySearchResults == null)
                {
                    // Decoy search skipped, read in the result file
                    decoySearchResults = DatabaseSearchResultData.ReadResultsFromFile(decoyOutputFilePath);
                    fdrCalculator = new FdrCalculator(targetSearchResults, decoySearchResults, multiplePeptidesPerScan);
                }
                else
                {
                    // Just use the objects
                    fdrCalculator = new FdrCalculator(targetSearchResults, decoySearchResults, multiplePeptidesPerScan);
                }

                if (fdrCalculator.HasError())
                {
                    ErrorMessage = fdrCalculator.ErrorMessage;
                    ReportWarning("Error computing FDR: " + fdrCalculator.ErrorMessage);
                    return false;
                }

                fdrCalculator.WriteTo(tdaOutputFilePath, Options.IncludeDecoyResults);
                var mzidWriter = new MzidResultsWriter(targetDb, _run, Options);
                mzidWriter.WriteResultsToMzid(fdrCalculator.FilteredResults, mzidOutputFilePath);
            }
            else
            {
                var results = targetSearchResults;
                var filePath = targetOutputFilePath;
                if (!Options.TargetDecoySearchMode.HasFlag(DatabaseSearchMode.Target))
                {
                    results = decoySearchResults;
                    filePath = decoyOutputFilePath;
                }
                if (results == null)
                {
                    results = DatabaseSearchResultData.ReadResultsFromFile(filePath);
                }

                var mzidWriter = new MzidResultsWriter(targetDb, _run, Options);
                mzidWriter.WriteResultsToMzid(results, mzidOutputFilePath);
            }
            progData.Report(100.0);

            OnStatusEvent("Done.");

            // Stop the overall progress reports.
            progReportTimer.Dispose();
            swAll.Stop();

            var elapsed = swAll.Elapsed;
            var minutes = elapsed.TotalMinutes - ((int)elapsed.TotalHours * 60);
            OnStatusEvent(string.Format("Total elapsed time for search: {0:f1} sec ({1}d {2}h {3:f2} min)", elapsed.TotalSeconds, elapsed.Days, elapsed.Hours, minutes));

            return true;
        }

        private List<DatabaseSearchResultData> RunDatabaseSearch(
            FastaDatabase searchDb,
            string outputFilePath,
            ISequenceFilter ms1Filter,
            string searchModeString,
            IProgress<ProgressData> progress)
        {
            var progData = new ProgressData(progress);
            var sw = new Stopwatch();
            var searchModeStringCap = char.ToUpper(searchModeString[0]) + searchModeString.Substring(1);

            UpdateStatus(string.Format("Reading the {0} database...", searchModeString), progData);
            sw.Reset();
            sw.Start();

            searchDb.Read();
            sw.Stop();
            OnStatusEvent(string.Format("Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds));

            var matches = new SortedSet<DatabaseSequenceSpectrumMatch>[_run.MaxLcScan + 1];
            if (Options.TagBasedSearch)
            {
                progData.StepRange(3);
                UpdateStatus(string.Format("Tag-based searching the {0} database", searchModeString), progData);
                sw.Reset();
                sw.Start();
                var progTag = new Progress<ProgressData>(p =>
                {
                    progData.StatusInternal = p.Status;
                    progData.Report(p.Percent);
                });
                RunTagBasedSearch(matches, searchDb, null, progTag);
                OnStatusEvent(string.Format("{1} database tag-based search elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds, searchModeStringCap));
            }
            progData.StepRange(88);

            UpdateStatus(string.Format("Searching the {0} database", searchModeString), progData);
            sw.Reset();
            sw.Start();
            var prog = new Progress<ProgressData>(p =>
            {
                progData.StatusInternal = p.Status;
                progData.Report(p.Percent);
            });
            RunSearch(matches, searchDb, ms1Filter, null, prog);
            OnStatusEvent(string.Format("{1} database search elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds, searchModeStringCap));

            // calculate spectral e-value using generating function
            UpdateStatus(string.Format("Calculating spectral E-values for {0}-spectrum matches", searchModeString), progData);
            sw.Reset();
            sw.Start();
            progData.StepRange(100);
            var progGen = new Progress<ProgressData>(p =>
            {
                progData.StatusInternal = p.Status;
                progData.Report(p.Percent);
            });

            var bestMatchesByScan = RunGeneratingFunction(matches, searchDb, null, progGen);

            var results = WriteResultsToFile(bestMatchesByScan, outputFilePath, searchDb);
            sw.Stop();

            OnStatusEvent(string.Format("{1}-spectrum match E-value calculation elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds, searchModeStringCap));
            return results;
        }

        private int[] _tagMs2ScanNum;

        private IEnumerable<AnnotationAndOffset> GetAnnotationsAndOffsets(
            FastaDatabase database,
            out long estimatedProteins,
            CancellationToken? cancellationToken = null)
        {
            var indexedDb = new IndexedDatabase(database);
            indexedDb.Read();

            estimatedProteins = indexedDb.EstimateTotalPeptides(
                Options.InternalCleavageMode,
                Options.MinSequenceLength, Options.MaxSequenceLength,
                Options.MaxNumNTermCleavages, Options.MaxNumCTermCleavages);

            IEnumerable<AnnotationAndOffset> annotationsAndOffsets;
            if (Options.InternalCleavageMode == InternalCleavageType.MultipleInternalCleavages)
            {
                annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzymeParallel(
                    Options.MinSequenceLength, Options.MaxSequenceLength, Options.MaxNumThreads, cancellationToken);
            }
            else if (Options.InternalCleavageMode == InternalCleavageType.NoInternalCleavage)
            {
                annotationsAndOffsets = indexedDb.IntactSequenceAnnotationsAndOffsets(
                    Options.MinSequenceLength, Options.MaxSequenceLength, Options.MaxNumCTermCleavages);
            }
            else
            {
                annotationsAndOffsets = indexedDb.SequenceAnnotationsAndOffsetsWithNTermOrCTermCleavageNoLargerThan(
                    Options.MinSequenceLength, Options.MaxSequenceLength, Options.MaxNumNTermCleavages, Options.MaxNumCTermCleavages);
            }

            return annotationsAndOffsets;
        }

        private void RunTagBasedSearch(
            SortedSet<DatabaseSequenceSpectrumMatch>[] matches,
            FastaDatabase db,
            CancellationToken? cancellationToken = null,
            IProgress<ProgressData> progress = null)
        {
            _tagSearchEngine.SetDatabase(db);

            var progData = new ProgressData(progress)
            {
                Status = "Tag-based Searching for matches"
            };

            var sw = new Stopwatch();

            long estimatedProteins = _tagMs2ScanNum.Length;
            OnStatusEvent("Number of spectra containing sequence tags: " + estimatedProteins);
            var numProteins = 0;
            var lastUpdate = DateTime.MinValue; // Force original update of 0%

            sw.Reset();
            sw.Start();

            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Options.MaxDegreeOfParallelism,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            Parallel.ForEach(_tagMs2ScanNum, pfeOptions, ms2ScanNum =>
            {
                var tagSeqMatches = _tagSearchEngine.RunSearch(ms2ScanNum);

                foreach (var tagSequenceMatch in tagSeqMatches)
                {
                    var offset = _tagSearchEngine.FastaDatabase.GetOffset(tagSequenceMatch.ProteinName);
                    if (offset == null)
                    {
                        continue;
                    }

                    var sequence = tagSequenceMatch.Sequence;
                    var numNTermCleavages = tagSequenceMatch.TagMatch.StartIndex;

                    var seqObj = Sequence.CreateSequence(sequence, tagSequenceMatch.TagMatch.ModificationText, Options.AminoAcidSet);
                    var precursorIon = new Ion(seqObj.Composition + Composition.H2O, tagSequenceMatch.TagMatch.Charge);

                    var prsm = new DatabaseSequenceSpectrumMatch(sequence, tagSequenceMatch.Pre, tagSequenceMatch.Post,
                                                                 ms2ScanNum, (long)offset, numNTermCleavages,
                                                                 tagSequenceMatch.TagMatch.Modifications,
                                                                 precursorIon, tagSequenceMatch.TagMatch.Score, db.IsDecoy)
                    {
                        ModificationText = tagSequenceMatch.TagMatch.ModificationText,
                    };

                    AddMatch(matches, ms2ScanNum, prsm);
                }

                SearchProgressReport(ref numProteins, ref lastUpdate, estimatedProteins, sw, progData, matches, "spectra");
            });

            OnStatusEvent(string.Format("Collected candidate matches: {0}", GetNumberOfMatches(matches)));

            progData.StatusInternal = string.Empty;
            progData.Report(100.0);
        }

        private void RunSearch(
            SortedSet<DatabaseSequenceSpectrumMatch>[] matches,
            FastaDatabase db,
            ISequenceFilter sequenceFilter,
            CancellationToken? cancellationToken = null,
            IProgress<ProgressData> progress = null)
        {
            var progData = new ProgressData(progress)
            {
                Status = "Searching for matches"
            };

            var sw = new Stopwatch();
            var annotationsAndOffsets = GetAnnotationsAndOffsets(db, out var estimatedProteins, cancellationToken);
            OnStatusEvent("Estimated Sequences: " + estimatedProteins.ToString("#,##0"));

            var numProteins = 0;
            var lastUpdate = DateTime.MinValue; // Force original update of 0%

            sw.Reset();
            sw.Start();
            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Options.MaxDegreeOfParallelism,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            var maxNumNTermCleavages = Options.InternalCleavageMode == InternalCleavageType.NoInternalCleavage ? Options.MaxNumNTermCleavages : 0;

            if (USE_PARALLEL_FOREACH && !DEBUG_MODE)
            {
                Parallel.ForEach(annotationsAndOffsets, pfeOptions, annotationAndOffset =>
                {
                    if (cancellationToken?.IsCancellationRequested == true)
                    {
                        //return matches;
                        return;
                    }

                    SearchProgressReport(ref numProteins, ref lastUpdate, estimatedProteins, sw, progData, matches);
                    SearchForMatches(annotationAndOffset, sequenceFilter, matches, maxNumNTermCleavages, db.IsDecoy, cancellationToken);
                });
            }
            else
            {
                foreach(var annotationAndOffset in annotationsAndOffsets)
                {
                    if (cancellationToken?.IsCancellationRequested == true)
                    {
                        //return matches;
                        return;
                    }

                    SearchProgressReport(ref numProteins, ref lastUpdate, estimatedProteins, sw, progData, matches);
                    SearchForMatches(annotationAndOffset, sequenceFilter, matches, maxNumNTermCleavages, db.IsDecoy, cancellationToken);

                    if (numProteins > DEBUG_MODE_PROTEINS_TO_SEARCH)
                    {
                        ConsoleMsgUtils.ShowWarning("Debug mode");
                        ConsoleMsgUtils.ShowWarning(string.Format("Exiting ForEach since {0} proteins have been searched", DEBUG_MODE_PROTEINS_TO_SEARCH));
                        Console.WriteLine();
                        break;
                    }
                }
            }

            OnStatusEvent(string.Format("Collected candidate matches: {0}", GetNumberOfMatches(matches)));

            progData.StatusInternal = string.Empty;
            progData.Report(100.0);
        }

        private void SearchProgressReport(
            ref int numProteins,
            ref DateTime lastUpdate,
            long estimatedProteins,
            Stopwatch sw,
            ProgressData progData,
            SortedSet<DatabaseSequenceSpectrumMatch>[] matches,
            string itemName = "proteins")
        {
            var tempNumProteins = Interlocked.Increment(ref numProteins) - 1;

            if (estimatedProteins < 1)
            {
                estimatedProteins = 1;
            }

            progData.StatusInternal = string.Format("Processing, {0} {1} done, {2:#0.0}% complete, {3:f1} sec elapsed",
                tempNumProteins,
                itemName,
                tempNumProteins / (double)estimatedProteins * 100.0,
                sw.Elapsed.TotalSeconds);
            progData.Report(tempNumProteins, estimatedProteins);

            int secondsThreshold;

            if (sw.Elapsed.TotalMinutes < 2)
            {
                secondsThreshold = 15;      // Every 15 seconds
            }
            else if (sw.Elapsed.TotalMinutes < 5)
            {
                secondsThreshold = 30;      // Every 30 seconds
            }
            else if (sw.Elapsed.TotalMinutes < 20)
            {
                secondsThreshold = 60;      // Every 1 minute
            }
            else
            {
                secondsThreshold = 300;     // Every 5 minutes
            }

            if (DateTime.UtcNow.Subtract(lastUpdate).TotalSeconds < secondsThreshold)
            {
                return;
            }

            lastUpdate = DateTime.UtcNow;

            string matchCountStats;
            if (matches != null)
            {
                var numMatches = GetNumberOfMatches(matches);
                if (numMatches == 1)
                {
                    matchCountStats = ", 1 match";
                }
                else
                {
                    matchCountStats = string.Format(", {0} matches", numMatches);
                }
            }
            else
            {
                matchCountStats = string.Empty;
            }

            // Processing, 633 proteins done, 34.2% complete, 60.4 sec elapsed, 48 matches
            OnStatusEvent(string.Format("Processing, {0} {1} done, {2:#0.0}% complete, {3:f1} sec elapsed{4}",
                                        tempNumProteins,
                                        itemName,
                                        tempNumProteins / (double)estimatedProteins * 100.0,
                                        sw.Elapsed.TotalSeconds,
                                        matchCountStats));
        }

        private void SearchForMatches(
            AnnotationAndOffset annotationAndOffset,
            ISequenceFilter sequenceFilter,
            SortedSet<DatabaseSequenceSpectrumMatch>[] matches,
            int maxNumNTermCleavages,
            bool isDecoy,
            CancellationToken? cancellationToken = null)
        {
            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Options.MaxDegreeOfParallelism,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            var annotation = annotationAndOffset.Annotation;
            var offset = annotationAndOffset.Offset;
            //var protein = db.GetProteinName(offset);
            var proteinSequence = annotation.Substring(2, annotation.Length - 4);
            var seqGraph = SequenceGraph.CreateGraph(Options.AminoAcidSet, AminoAcid.ProteinNTerm, proteinSequence,
                AminoAcid.ProteinCTerm);

            if (seqGraph == null)
            {
                return; // No matches will be found without a sequence graph.
            }

            for (var numNTermCleavages = 0; numNTermCleavages <= maxNumNTermCleavages; numNTermCleavages++)
            {
                if (numNTermCleavages > 0)
                {
                    seqGraph.CleaveNTerm();
                }

                var numProteoforms = seqGraph.GetNumProteoformCompositions();
                var modCombs = seqGraph.GetModificationCombinations();
                for (var modIndex = 0; modIndex < numProteoforms; modIndex++)
                {
                    seqGraph.SetSink(modIndex);
                    var proteinCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();
                    var sequenceMass = proteinCompositionWithH2O.Mass;

                    if (sequenceMass < Options.MinSequenceMass || sequenceMass > Options.MaxSequenceMass)
                    {
                        continue;
                    }

                    var modCombinations = modCombs[modIndex];
                    var ms2ScanNums = Options.ScanNumbers ?? sequenceFilter.GetMatchingMs2ScanNums(sequenceMass);
                    var filter = sequenceFilter as Ms1FtFilter;
                    var featureIds = filter?.GetMatchingFeatureIds(sequenceMass).ToList() ?? new List<int>();

                    var nTermCleavages = numNTermCleavages;

                    Parallel.ForEach(ms2ScanNums, pfeOptions, ms2ScanNum =>
                    {
                        if (ms2ScanNum > _ms2ScanNums.Last() || ms2ScanNum < _ms2ScanNums.First())
                        {
                            return;
                        }

                        var isoTargetMz = _isolationWindowTargetMz[ms2ScanNum];
                        if (!(isoTargetMz > 0))
                        {
                            return;
                        }

                        var charge = (int)Math.Round(sequenceMass / (isoTargetMz - Constants.Proton));

                        //var scorer = _ms2ScorerFactory2.GetScorer(ms2ScanNum, sequenceMass, charge, Options.ActivationMethod);
                        IScorer scorer;
                        if (InformedProteomics.Backend.Utils.FlipSwitch.UseFlipScoring)
                        {
                            scorer = _ms2ScorerFactory2.GetScorer(ms2ScanNum);
                        }
                        else
                        {
                            scorer = _ms2ScorerFactory2.GetMs2Scorer(ms2ScanNum);
                        }

                        var score = seqGraph.GetFragmentScore(scorer);

                        var precursorIon = new Ion(proteinCompositionWithH2O, charge);
                        var sequence = proteinSequence.Substring(nTermCleavages);
                        var pre = nTermCleavages == 0 ? annotation[0] : annotation[nTermCleavages + 1];
                        var post = annotation[annotation.Length - 1];

                        // MS1 FeatureIds for MS2 scans
                        var ms2FeatureIds = new List<int>();
                        if (filter != null)
                        {
                            foreach (var featureId in featureIds)
                            {
                                var scanRange = filter.Ms1FtIndexToScanRange[featureId];
                                if (ms2ScanNum >= scanRange.Item1 && ms2ScanNum <= scanRange.Item2)
                                {
                                    ms2FeatureIds.Add(featureId);
                                }
                            }
                        }
                        else
                        {
                            ms2FeatureIds = featureIds;
                        }

                        var prsm = new DatabaseSequenceSpectrumMatch(sequence, pre, post, ms2ScanNum, offset, nTermCleavages,
                            modCombinations, precursorIon, score, isDecoy, featureId: Math.Max(ms2FeatureIds.FirstOrDefault(), 1));

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
                    matches[ms2ScanNum] = new SortedSet<DatabaseSequenceSpectrumMatch> { prsm };
                }
                else // already exists
                {
                    var existingMatches = matches[ms2ScanNum];
                    //var maxScore = existingMatches.Max.Score;
                    if (existingMatches.Count < Options.MatchesPerSpectrumToKeepInMemory)
                    {
                        //if (!(maxScore*0.7 < prsm.Score)) return;
                        existingMatches.Add(prsm);
                    }
                    else
                    {
                        var minScore = existingMatches.Min.Score;
                        if (!(prsm.Score > minScore))
                        {
                            return;
                        }

                        existingMatches.Add(prsm);
                        existingMatches.Remove(existingMatches.Min);
                    }
                    //if (MatchesPerSpectrumToKeepInMemory > 1) existingMatches.RemoveWhere(mt => mt.Score < maxScore * 0.7);
                }
            }
        }

        private SequenceTagGenerator GetSequenceTagGenerator(CancellationToken? cancellationToken = null, IProgress<ProgressData> progress = null)
        {
            var sequenceTagGen = new SequenceTagGenerator(_run, new Tolerance(5));
            var scanNums = _ms2ScanNums;

            var progData = new ProgressData(progress)
            {
                Status = "Generating sequence tags"
            };

            var sw = new Stopwatch();

            // Re-score and Estimate #proteins for GF calculation
            long estimatedProteins = scanNums.Count;
            OnStatusEvent("Number of spectra: " + estimatedProteins);
            var numProteins = 0;
            var lastUpdate = DateTime.MinValue; // Force original update of 0%
            sw.Reset();
            sw.Start();

            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Options.MaxDegreeOfParallelism,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            Parallel.ForEach(scanNums, pfeOptions, scanNum =>
            {
                sequenceTagGen.Generate(scanNum);
                SearchProgressReport(ref numProteins, ref lastUpdate, estimatedProteins, sw, progData, null, "spectra");
            });

            progData.StatusInternal = string.Empty;
            progData.Report(100.0);
            OnStatusEvent("Generated sequence tags: " + sequenceTagGen.NumberOfGeneratedTags());
            return sequenceTagGen;
        }

        private LinkedList<Tuple<double, ScoreDistribution>>[] _cachedScoreDistributions;

        private Dictionary<int, List<DatabaseSequenceSpectrumMatch>> RunGeneratingFunction(
            IReadOnlyList<SortedSet<DatabaseSequenceSpectrumMatch>> sortedMatches,
            FastaDatabase searchDb,
            CancellationToken? cancellationToken = null,
            IProgress<ProgressData> progress = null)
        {
            var progData = new ProgressData(progress)
            {
                Status = "Calculating spectral E-values for matches"
            };

            var sw = new Stopwatch();

            long estimatedSequences = 0;
            var matches = new LinkedList<DatabaseSequenceSpectrumMatch>[sortedMatches.Count];

            var currentTask = "?";

            try
            {
                currentTask = "Validate _cachedScoreDistributions";

                if (_cachedScoreDistributions == null)
                {
                    if (_run == null)
                    {
                        ReportError("_run is null in RunGeneratingFunction; cannot initialize _cachedScoreDistributions");
                    }

                    if (_ms2ScanNums == null)
                    {
                        ReportError("_ms2ScanNums is null in RunGeneratingFunction; cannot initialize _cachedScoreDistributions");
                    }

                    _cachedScoreDistributions = new LinkedList<Tuple<double, ScoreDistribution>>[_run.MaxLcScan + 1];
                    foreach (var scanNum in _ms2ScanNums)
                    {
                        _cachedScoreDistributions[scanNum] = new LinkedList<Tuple<double, ScoreDistribution>>();
                    }
                }

                currentTask = "Instantiate InformedTopDownScorer";

                var topDownScorer = new InformedTopDownScorer(_run, Options.AminoAcidSet, Options.MinProductIonCharge, Options.MaxProductIonCharge, Options.ProductIonTolerance, activationMethod: Options.ActivationMethod);

                currentTask = "Re-score and Estimate #proteins for GF calculation";

                foreach (var scanNum in _ms2ScanNums)
                {
                    var prsms = sortedMatches[scanNum];
                    if (prsms == null)
                    {
                        continue;
                    }

                    if (!(_run.GetSpectrum(scanNum) is ProductSpectrum spec))
                    {
                        continue;
                    }

                    if (spec.Peaks.Length == 0)
                    {
                        continue;
                    }

                    currentTask = "Looping over PRSMs for scan " + scanNum;

                    var matchIndex = 0;
                    foreach (var match in prsms)
                    {
                        var sequence = match.Sequence;

                        var ion = match.Ion;

                        if (ion == null)
                        {
                            ReportError("Ion is null for index " + matchIndex + " in scan " + scanNum);
                            continue;
                        }

                        // Re-scoring
                        var scores = topDownScorer.GetScores(spec, sequence, ion.Composition, ion.Charge, scanNum);
                        if (scores == null)
                        {
                            ReportWarning("scores is null for index " + matchIndex + " in scan " + scanNum);
                            continue;
                        }

                        match.Score = scores.Score;
                        match.ModificationText = scores.Modifications;
                        match.NumMatchedFragments = scores.NumMatchedFrags;
                        if (match.Score > CompositeScorer.ScoreParam.Cutoff)
                        {
                            if (matches[scanNum] == null)
                            {
                                matches[scanNum] = new LinkedList<DatabaseSequenceSpectrumMatch>();
                            }

                            matches[scanNum].AddLast(match);
                        }

                        matchIndex++;
                    }

                    if (matches[scanNum] != null)
                    {
                        estimatedSequences += matches[scanNum].Count;
                    }
                }

                currentTask = "Parallel.ForEach";

                OnStatusEvent("Estimated matched sequences: " + estimatedSequences.ToString("#,##0"));
            }
            catch (Exception ex)
            {
                ReportError(string.Format("Exception while {0} in RunGeneratingFunction: {1}", currentTask, ex.Message), ex);
            }

            var numProteins = 0;
            var lastUpdate = DateTime.MinValue; // Force original update of 0%
            sw.Reset();
            sw.Start();

            var scanNums = _ms2ScanNums.Where(scanNum => matches[scanNum] != null).ToArray();

            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Options.MaxDegreeOfParallelism,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            Parallel.ForEach(scanNums, pfeOptions, scanNum =>
            {
                currentTask = "Inside Parallel.ForEach";
                try
                {
                    var scoreDistributions = _cachedScoreDistributions[scanNum];
                    foreach (var match in matches[scanNum])
                    {
                        var currentIteration = "for scan " + scanNum + " and mass " + match.Ion.Composition.Mass;
                        currentTask = "Calling GetMs2ScoringGraph " + currentIteration;

                        var graph = _ms2ScorerFactory2.GetMs2ScoringGraph(scanNum, match.Ion.Composition.Mass);
                        if (graph == null)
                        {
                            continue;
                        }

                        currentTask = "Calling ComputeGeneratingFunction " + currentIteration;

                        var scoreDist = (from distribution in scoreDistributions
                                         where Math.Abs(distribution.Item1 - match.Ion.Composition.Mass) < Options.PrecursorIonTolerance.GetToleranceAsMz(match.Ion.Composition.Mass)
                                         select distribution.Item2).FirstOrDefault();
                        if (scoreDist == null)
                        {
                            var gf = new GeneratingFunction(graph);
                            gf.ComputeGeneratingFunction();
                            scoreDist = gf.GetScoreDistribution();
                            scoreDistributions.AddLast(new Tuple<double, ScoreDistribution>(match.Ion.Composition.Mass, scoreDist));
                        }

                        currentTask = "Calling GetSpectralEValue " + currentIteration + " and score " + (int)match.Score;
                        match.SpecEvalue = scoreDist.GetSpectralEValue(match.Score);

                        currentTask = "Reporting progress " + currentIteration;
                        SearchProgressReport(ref numProteins, ref lastUpdate, estimatedSequences, sw, progData, null);

                        if (DEBUG_MODE && numProteins > DEBUG_MODE_PROTEINS_TO_SEARCH)
                        {
                            ConsoleMsgUtils.ShowWarning("Debug mode");
                            ConsoleMsgUtils.ShowWarning(string.Format(
                                "Exiting Parallel.ForEach since {0} proteins have been searched", DEBUG_MODE_PROTEINS_TO_SEARCH));
                            Console.WriteLine();
                            break;
                        }
                    }
                }
                catch (Exception ex)
                {
                    ReportError(string.Format("Exception while {0} in RunGeneratingFunction: {1}", currentTask, ex.Message), ex);
                }
            });

            // Values in this SortedSet are Sequence_Modification_ProteinName
            var storedSequences = new SortedSet<string>();

            // Keys in this dictionary are scan number
            // Values are the list of matches for each scan
            var finalMatches = new Dictionary<int, List<DatabaseSequenceSpectrumMatch>>();

            foreach (var scanNum in scanNums)
            {
                var matchesForScan = matches[scanNum].OrderBy(m => m.SpecEvalue).ToList();

                var matchesToStore = new List<DatabaseSequenceSpectrumMatch>();

                storedSequences.Clear();

                double comparisonSpecEValue = -1;
                var ms1FeatureId = 0;

                for (var i = 0; i < matchesForScan.Count; i++)
                {
                    // Multiple matches can have the same sequence, same protein, and same score
                    // Typically they will have differing Offsets
                    // Check for this and avoid duplicates

                    var proteinName = searchDb.GetProteinName(matchesForScan[i].Offset);
                    var sequenceKey = matchesForScan[i].Sequence + "_" + matchesForScan[i].ModificationText ?? string.Empty + "_" + proteinName;

                    if (!storedSequences.Contains(sequenceKey))
                    {
                        matchesToStore.Add(matchesForScan[i]);
                        storedSequences.Add(sequenceKey);
                    }

                    if (matchesToStore.Count < Options.MatchesPerSpectrumToReport)
                    {
                        continue;
                    }

                    if (ms1FeatureId == 0 && matchesForScan[i].FeatureId > 0)
                    {
                        ms1FeatureId = matchesForScan[i].FeatureId;
                    }

                    if (i == matchesForScan.Count - 1)
                    {
                        // No more matches
                        break;
                    }

                    if (i == 0)
                    {
                        comparisonSpecEValue = matchesForScan[0].SpecEvalue;
                    }

                    // Compare the next lowest scoring match's score to the comparison Spec EValue
                    if (Math.Abs(comparisonSpecEValue - matchesForScan[i + 1].SpecEvalue) > double.Epsilon)
                    {
                        // Scores are sufficiently different
                        break;
                    }
                }

                // If any of the matches have a FeatureId of 0, change it to ms1FeatureId
                foreach (var match in matchesToStore)
                {
                    if (match.FeatureId == 0)
                    {
                        match.UpdateFeatureId(ms1FeatureId);
                    }
                }

                finalMatches.Add(scanNum, matchesToStore);
            }

            progData.StatusInternal = string.Empty;
            progData.Report(100.0);
            return finalMatches;
        }

        private int GetNumberOfMatches(SortedSet<DatabaseSequenceSpectrumMatch>[] matches)
        {
            var nMatches = 0;
            lock (matches)
            {
                // ReSharper disable once UselessBinaryOperation
                nMatches += _ms2ScanNums.Where(scanNum => matches[scanNum] != null).Sum(scanNum => matches[scanNum].Count);
            }
            return nMatches;
        }

        private bool ResultsFileHasData(string targetOutputFilePath)
        {
            try
            {
                var outputFile = new FileInfo(targetOutputFilePath);
                if (!outputFile.Exists || outputFile.Length == 0)
                {
                    return false;
                }

                using (var reader = new StreamReader(new FileStream(outputFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    var resultCount = 0;
                    while (!reader.EndOfStream)
                    {
                        var dataLine = reader.ReadLine();
                        if (string.IsNullOrWhiteSpace(dataLine))
                        {
                            continue;
                        }

                        resultCount++;
                    }

                    if (resultCount > 1)
                    {
                        return true;
                    }

                    if (resultCount == 1)
                    {
                        ReportWarning("Results file only has a header line; will re-generate " + targetOutputFilePath);
                    }

                    ReportWarning("Results file is empty: " + targetOutputFilePath);
                }

                return false;
            }
            catch (Exception ex)
            {
                ReportWarning(string.Format(
                                  "Error validating existing results in file {0}: {1}",
                                  targetOutputFilePath, ex.Message));
                return false;
            }
        }

        private List<DatabaseSearchResultData> WriteResultsToFile(
            IReadOnlyDictionary<int, List<DatabaseSequenceSpectrumMatch>> matchesByScan,
            string outputFilePath,
            FastaDatabase database)
        {
            var results = CreateResults(matchesByScan, database).ToList();
            DatabaseSearchResultData.WriteResultsToFile(outputFilePath, results);
            return results;
        }

        private IEnumerable<DatabaseSearchResultData> CreateResults(
            IReadOnlyDictionary<int, List<DatabaseSequenceSpectrumMatch>> matchesByScan,
            FastaDatabase database)
        {
            foreach (var scanNum in _ms2ScanNums)
            {
                if (!matchesByScan.TryGetValue(scanNum, out var matches))
                {
                    continue;
                }

                foreach (var match in matches)
                {
                    var start = database.GetOneBasedPositionInProtein(match.Offset) + 1 + match.NumNTermCleavages;
                    var proteinName = database.GetProteinName(match.Offset);
                    var result = new DatabaseSearchResultData
                    {
                        ScanNum = scanNum,
                        Pre = match.Pre.ToString(),
                        Sequence = match.Sequence,
                        Post = match.Post.ToString(),
                        Modifications = match.ModificationText,
                        Composition = match.Ion.Composition.ToString(),
                        ProteinName = proteinName,
                        ProteinLength = database.GetProteinLength(proteinName),
                        ProteinDescription = database.GetProteinDescription(match.Offset),
                        Start = start,
                        End = start + match.Sequence.Length - 1,
                        Charge = match.Ion.Charge,
                        MostAbundantIsotopeMz = match.Ion.GetMostAbundantIsotopeMz(),
                        Mass = match.Ion.Composition.Mass,
                        NumMatchedFragments = match.NumMatchedFragments,
                        Probability = CompositeScorer.GetProbability(match.Score),
                        SpecEValue = match.SpecEvalue,
                        EValue = match.SpecEvalue * database.GetNumEntries(),
                        Ms1Feature = match.FeatureId
                    };

                    yield return result;
                }
            }
        }

        /// <summary>
        /// Show an error message at the console then throw an exception
        /// </summary>
        /// <param name="errMsg">Error message</param>
        /// <param name="ex">Exception (can be null)</param>
        /// <param name="throwException">True to throw an exception</param>
        private void ReportError(string errMsg, Exception ex = null, bool throwException = true)
        {
            OnErrorEvent(errMsg, ex);

            if (!throwException)
            {
                return;
            }

            if (ex == null)
            {
                throw new Exception(errMsg);
            }

            throw new Exception(errMsg, ex);
        }

        private void ReportWarning(string msg)
        {
            OnWarningEvent(msg);
        }

        private void UpdateStatus(string message, ProgressData progData)
        {
            OnStatusEvent(message);
            progData.Status = message;
        }
    }
}
