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
using InformedProteomics.Backend.Utils;
using InformedProteomics.FeatureFinding;
using InformedProteomics.FeatureFinding.FeatureDetection;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;
using InformedProteomics.TopDown.Scoring;
using InformedProteomics.TopDown.TagBasedSearch;

using InformedProteomics.Scoring.Interfaces;

namespace InformedProteomics.TopDown.Execution
{
    public class IcTopDownLauncher
    {
        //public const int NumMatchesPerSpectrum = 1;
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

            //Options.MinSequenceLength = 21;
            //Options.MaxSequenceLength = 300;
            //Options.MaxNumNTermCleavages = 1;
            //Options.MaxNumCTermCleavages = 0;
            //Options.MinPrecursorIonCharge = 2;
            //Options.MaxPrecursorIonCharge = 60;
            //Options.MinProductIonCharge = 1;
            //Options.MaxProductIonCharge = 20;
            //Options.MinSequenceMass = 2000.0;
            //Options.MaxSequenceMass = 50000.0;
            //Options.PrecursorIonTolerance = new Tolerance(10);
            //Options.ProductIonTolerance = new Tolerance(10);
            //Options.RunTargetDecoyAnalysis = DatabaseSearchMode.Both;
            //Options.SearchMode = InternalCleavageType.SingleInternalCleavage;
            //Options.MaxNumThreads = 4;
            //Options.ScanNumbers = null;
            //Options.NumMatchesPerSpectrum = 3;
            //Options.TagBasedSearch = true;
        }

        public IcTopDownLauncher(MsPfParameters options)
        {
            ErrorMessage = string.Empty;

            Options = options;
        }

        public string ErrorMessage { get; private set; }
        public MsPfParameters Options { get; }

        [Obsolete("Use options!", true)]
        public string SpecFilePath => Options.SpecFilePath;

        [Obsolete("Use options!", true)]
        public string DatabaseFilePath => Options.DatabaseFilePath;

        [Obsolete("Use options!", true)]
        public string OutputDir => Options.OutputDir;

        [Obsolete("Use options!", true)]
        public AminoAcidSet AminoAcidSet => Options.AminoAcidSet;

        [Obsolete("Use options!", true)]
        public string FeatureFilePath { get => Options.FeatureFilePath;
            set => Options.FeatureFilePath = value;
        }

        /// <remarks>default 21</remarks>
        [Obsolete("Use options!", true)]
        public int MinSequenceLength { get => Options.MinSequenceLength;
            set => Options.MinSequenceLength = value;
        }

        /// <remarks>default 300</remarks>
        [Obsolete("Use options!", true)]
        public int MaxSequenceLength { get => Options.MaxSequenceLength;
            set => Options.MaxSequenceLength = value;
        }

        /// <remarks>default 1</remarks>
        [Obsolete("Use options!", true)]
        public int MaxNumNTermCleavages { get => Options.MaxNumNTermCleavages;
            set => Options.MaxNumNTermCleavages = value;
        }

        /// <remarks>default 0</remarks>
        [Obsolete("Use options!", true)]
        public int MaxNumCTermCleavages { get => Options.MaxNumCTermCleavages;
            set => Options.MaxNumCTermCleavages = value;
        }

        /// <remarks>default 2</remarks>
        [Obsolete("Use options!", true)]
        public int MinPrecursorIonCharge { get => Options.MinPrecursorIonCharge;
            set => Options.MinPrecursorIonCharge = value;
        }

        /// <remarks>default 60</remarks>
        [Obsolete("Use options!", true)]
        public int MaxPrecursorIonCharge { get => Options.MaxPrecursorIonCharge;
            set => Options.MaxPrecursorIonCharge = value;
        }

        /// <remarks>default 2000</remarks>
        [Obsolete("Use options!", true)]
        public double MinSequenceMass { get => Options.MinSequenceMass;
            set => Options.MinSequenceMass = value;
        }

        /// <remarks>default 50000</remarks>
        [Obsolete("Use options!", true)]
        public double MaxSequenceMass { get => Options.MaxSequenceMass;
            set => Options.MaxSequenceMass = value;
        }

        /// <remarks>default 1</remarks>
        [Obsolete("Use options!", true)]
        public int MinProductIonCharge { get => Options.MinProductIonCharge;
            set => Options.MinProductIonCharge = value;
        }

        /// <remarks>default 20</remarks>
        [Obsolete("Use options!", true)]
        public int MaxProductIonCharge { get => Options.MaxProductIonCharge;
            set => Options.MaxProductIonCharge = value;
        }

        /// <remarks>default 10 ppm</remarks>
        [Obsolete("Use options!", true)]
        public Tolerance PrecursorIonTolerance { get => Options.PrecursorIonTolerance;
            set => Options.PrecursorIonTolerance = value;
        }

        /// <remarks>default 10 ppm</remarks>
        [Obsolete("Use options!", true)]
        public Tolerance ProductIonTolerance { get => Options.ProductIonTolerance;
            set => Options.ProductIonTolerance = value;
        }

        /// <remarks>default true
        /// true: target and decoy, false: target only, null: decoy only</remarks>
        [Obsolete("Use options!", true)]
        public bool? RunTargetDecoyAnalysisBool { get => Options.RunTargetDecoyAnalysisBool;
            set => Options.RunTargetDecoyAnalysisBool = value;
        }

        /// <remarks>default Both</remarks>
        [Obsolete("Use options!", true)]
        public DatabaseSearchMode RunTargetDecoyAnalysis { get => Options.TargetDecoySearchMode;
            set => Options.TargetDecoySearchMode = value;
        }

        [Obsolete("Use options!", true)]
        public bool TagBasedSearch { get => Options.TagBasedSearch;
            set => Options.TagBasedSearch = value;
        }

        /// <summary>
        /// Specific MS2 scan numbers to process
        /// </summary>
        [Obsolete("Use options!", true)]
        public IEnumerable<int> ScanNumbers { get => Options.ScanNumbers;
            set => Options.ScanNumbers = value;
        }

        /// <remarks>default 3</remarks>
        [Obsolete("Use options!", true)]
        public int NumMatchesPerSpectrum { get => Options.NumMatchesPerSpectrum;
            set => Options.NumMatchesPerSpectrum = value;
        }

        /// <remarks>default 4</remarks>
        [Obsolete("Use options!", true)]
        public int MaxNumThreads { get => Options.MaxNumThreads;
            set => Options.MaxNumThreads = value;
        }

        /// <remarks>default 10 ppm</remarks>
        [Obsolete("Use options!", true)]
        public double PrecursorIonTolerancePpm { get => Options.PrecursorIonTolerancePpm;
            set => Options.PrecursorIonTolerancePpm = value;
        }

        /// <remarks>default 10 ppm</remarks>
        [Obsolete("Use options!", true)]
        public double ProductIonTolerancePpm { get => Options.ProductIonTolerancePpm;
            set => Options.ProductIonTolerancePpm = value;
        }

        /// <summary>
        /// 0: all internal sequences,
        /// 1: #NCleavages &lt;= Max OR Cleavages &lt;= Max (Default)
        /// 2: 1: #NCleavages &lt;= Max AND Cleavages &lt;= Max
        /// </summary>
        /// <remarks>default 1</remarks>
        [Obsolete("Use options!", true)]
        public int SearchModeInt { get => Options.SearchModeInt;
            set => Options.SearchModeInt = value;
        }

        /// <remarks>default SingleInternalCleavage</remarks>
        [Obsolete("Use options!", true)]
        public InternalCleavageType SearchMode { get => Options.InternalCleavageMode;
            set => Options.InternalCleavageMode = value;
        }

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

        public void PrintOverallProgressReport(object searchObj)
        {
            var searchClass = (IcTopDownLauncher)searchObj;
            var progData = searchClass.searchProgressData;
            var timeElapsed = searchClass.searchStopwatch.Elapsed;
            var minutes = timeElapsed.TotalMinutes - ((int) timeElapsed.TotalHours * 60);

            Console.WriteLine(@"Total Progress: {0:F2}%, {1}d {2}h {3:F2}m elapsed, Current Task: {4}", progData.Percent, timeElapsed.Days, timeElapsed.Hours, minutes, progData.Status);
        }

        public bool RunSearch(double corrThreshold = 0.7, CancellationToken? cancellationToken = null, IProgress<ProgressData> progress = null)
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
            // Output a progress report every 5 minutes...
            progReportTimer = new Timer(PrintOverallProgressReport, this, 0, 1000 * 60 * 5);

            if (string.Equals(Path.GetExtension(Options.SpecFilePath), ".pbf", StringComparison.InvariantCultureIgnoreCase))
                Console.WriteLine(@"Reading pbf file...");
            else
                Console.WriteLine(@"Creating and loading pbf file...");

            progData.StepRange(5.0, "Reading spectra file");
            sw.Start();

            _run = PbfLcMsRun.GetLcMsRun(Options.SpecFilePath, 0, 0, prog);

            var minMs1Scan = int.MinValue;
            var maxMs1Scan = int.MaxValue;

            // Retrieve the list of MS2 scans
            _ms2ScanNums = _run.GetScanNumbers(2).ToList();

            if (Options.ScanNumbers != null && Options.ScanNumbers.Any())
            {
                // Filter the MS2 scans using ScanNumbers
                _ms2ScanNums = _ms2ScanNums.Intersect(Options.ScanNumbers).ToList();

                minMs1Scan = _ms2ScanNums.Min() - 100;
                maxMs1Scan = _ms2ScanNums.Max() + 100;
            }

            _isolationWindowTargetMz = new double[_run.MaxLcScan + 1];
            foreach (var ms2Scan in _ms2ScanNums)
            {
                var ms2Spec = _run.GetSpectrum(ms2Scan) as ProductSpectrum;
                if (ms2Spec == null) continue;
                _isolationWindowTargetMz[ms2Scan] = ms2Spec.IsolationWindow.IsolationWindowTargetMz;
            }

            sw.Stop();
            Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

            progData.StepRange(10, "Reading Fasta File");
            Console.WriteLine(progData.Status);

            // Target database
            var targetDb = new FastaDatabase(Options.DatabaseFilePath);
            targetDb.Read();

            progData.StepRange(20.0);
            ISequenceFilter ms1Filter;
            if (Options.ScanNumbers != null && Options.ScanNumbers.Any())
            {
                ms1Filter = new SelectedMsMsFilter(Options.ScanNumbers);
            }
            else if (string.IsNullOrWhiteSpace(Options.FeatureFilePath))
            {
                // Checks whether SpecFileName.ms1ft exists
                var ms1FtFilePath = MassSpecDataReaderFactory.ChangeExtension(Options.SpecFilePath, LcMsFeatureFinderLauncher.FileExtension);
                if (!File.Exists(ms1FtFilePath))
                {
                    Console.WriteLine(@"Running ProMex...");
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
                Console.Write(@"Reading ProMex results...");
                ms1Filter = new Ms1FtFilter(_run, Options.PrecursorIonTolerance, ms1FtFilePath, -10);
            }
            else
            {
                sw.Reset();
                sw.Start();
                var extension = Path.GetExtension(Options.FeatureFilePath);
                if (extension.ToLower().Equals(".csv"))
                {
                    Console.Write(@"Reading ICR2LS/Decon2LS results...");
                    ms1Filter = new IsosFilter(_run, Options.PrecursorIonTolerance, Options.FeatureFilePath);
                }
                else if (extension.ToLower().Equals(".ms1ft"))
                {
                    Console.Write(@"Reading ProMex results...");
                    ms1Filter = new Ms1FtFilter(_run, Options.PrecursorIonTolerance, Options.FeatureFilePath, -10);
                }
                else if (extension.ToLower().Equals(".msalign"))
                {
                    Console.Write(@"Reading MS-Align+ results...");
                    ms1Filter = new MsDeconvFilter(_run, Options.PrecursorIonTolerance, Options.FeatureFilePath);
                }
                else ms1Filter = null; //new Ms1FeatureMatrix(_run);
            }

            sw.Stop();
            Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

            // pre-generate deconvoluted spectra for scoring
            _massBinComparer = new FilteredProteinMassBinning(Options.AminoAcidSet, Options.MaxSequenceMass + 1000);

            //this.fragmentScorerFactory = new CompositionScorerFactory(_run, true);

            sw.Reset();

            Console.WriteLine(@"Generating deconvoluted spectra for MS/MS spectra...");
            sw.Start();
            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Options.MaxNumThreads,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            // Deconvolute spectra
            var deconvoluter = new Deconvoluter(Options.MinProductIonCharge, Options.MaxProductIonCharge, 2, 0.1, Options.ProductIonTolerance);
            var lcmsRunDeconvoluter = new LcmsRunDeconvoluter(_run, deconvoluter, 2, pfeOptions.MaxDegreeOfParallelism);
            var deconvolutedRun = new DPbfLcMsRun(Options.SpecFilePath, lcmsRunDeconvoluter, keepDataReaderOpen: true);

            _ms2ScorerFactory2 = new CompositeScorerFactory(deconvolutedRun, _massBinComparer, Options.AminoAcidSet,
                                                   Options.MinProductIonCharge, Options.MaxProductIonCharge, Options.ProductIonTolerance);
            //Parallel.ForEach(_ms2ScanNums, pfeOptions, ms2ScanNum =>
            //{
            //    _ms2ScorerFactory2.GetScorer(ms2ScanNum, activationMethod: Options.ActivationMethod);
            //});

            sw.Stop();
            Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

            // Generate sequence tags for all MS/MS spectra
            if (Options.TagBasedSearch)
            {
                progData.StepRange(25.0, "Generating Sequence Tags");

                sw.Reset();
                Console.WriteLine(@"Generating sequence tags for MS/MS spectra...");
                sw.Start();
                var seqTagGen = GetSequenceTagGenerator();
                _tagMs2ScanNum = seqTagGen.GetMs2ScanNumsContainingTags().ToArray();
                sw.Stop();
                Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

                var peakMatrix = new LcMsPeakMatrix(_run, ms1Filter, minMs1Scan: minMs1Scan, maxMs1Scan: maxMs1Scan);
                _tagSearchEngine = new ScanBasedTagSearchEngine(
                    _run, seqTagGen, peakMatrix, targetDb, Options.ProductIonTolerance, Options.AminoAcidSet,
                    _ms2ScorerFactory2, ScanBasedTagSearchEngine.DefaultMinMatchedTagLength,
                    Options.MaxSequenceMass, Options.MinProductIonCharge, Options.MaxProductIonCharge);
            }

            var specFileName = MassSpecDataReaderFactory.RemoveExtension(Path.GetFileName(Options.SpecFilePath));
            var targetOutputFilePath = Path.Combine(Options.OutputDir, specFileName + TargetFileNameEnding);
            var decoyOutputFilePath = Path.Combine(Options.OutputDir, specFileName + DecoyFileNameEnding);
            var tdaOutputFilePath = Path.Combine(Options.OutputDir, specFileName + TdaFileNameEnding);
            var mzidOutputFilePath = Path.Combine(Options.OutputDir, specFileName + MzidFileNameEnding);

            progData.StepRange(60.0, "Running Target search");
            List<DatabaseSearchResultData> targetSearchResults = null;

            if (Options.TargetDecoySearchMode.HasFlag(DatabaseSearchMode.Target) && !File.Exists(targetOutputFilePath))
            {
                targetSearchResults = RunDatabaseSearch(targetDb, targetOutputFilePath, ms1Filter, "target", prog);
            }
            else if (File.Exists(targetOutputFilePath))
            {
                Console.WriteLine(@"Target results file '{0}' exists; skipping target search.", targetOutputFilePath);
            }

            progData.StepRange(95.0, "Running Decoy search"); // total to 95%
            List<DatabaseSearchResultData> decoySearchResults = null;

            if (Options.TargetDecoySearchMode.HasFlag(DatabaseSearchMode.Decoy) && !File.Exists(decoyOutputFilePath))
            {
                var decoyDb = targetDb.Decoy(null, true);
                decoySearchResults = RunDatabaseSearch(decoyDb, decoyOutputFilePath, ms1Filter, "decoy", prog);
            }
            else if (File.Exists(decoyOutputFilePath))
            {
                Console.WriteLine(@"Decoy results file '{0}' exists; skipping decoy search.", decoyOutputFilePath);
            }

            progData.StepRange(100.0, "Writing combined results file");
            if (Options.TargetDecoySearchMode.HasFlag(DatabaseSearchMode.Both))
            {
                // Add "Qvalue" and "PepQValue"
                FdrCalculator fdrCalculator;
                if (targetSearchResults == null && decoySearchResults == null)
                {
                    // Objects not populated, try to read in the files.
                    fdrCalculator = new FdrCalculator(targetOutputFilePath, decoyOutputFilePath);
                }
                else if (targetSearchResults == null)
                {
                    // Target search skipped, read in the result file
                    targetSearchResults = DatabaseSearchResultData.ReadResultsFromFile(targetOutputFilePath);
                    fdrCalculator = new FdrCalculator(targetSearchResults, decoySearchResults);
                }
                else if (decoySearchResults == null)
                {
                    // Decoy search skipped, read in the result file
                    decoySearchResults = DatabaseSearchResultData.ReadResultsFromFile(decoyOutputFilePath);
                    fdrCalculator = new FdrCalculator(targetSearchResults, decoySearchResults);
                }
                else
                {
                    // Just use the objects
                    fdrCalculator = new FdrCalculator(targetSearchResults, decoySearchResults);
                }

                if (fdrCalculator.HasError())
                {
                    ErrorMessage = fdrCalculator.ErrorMessage;
                    Console.WriteLine(@"Error computing FDR: " + fdrCalculator.ErrorMessage);
                    return false;
                }

                fdrCalculator.WriteTo(tdaOutputFilePath);
                var mzidWriter = new MzidResultsWriter(targetDb, _run, Options);
                mzidWriter.WriteResultsToMzid(fdrCalculator.FilteredResults, mzidOutputFilePath);
            }
            progData.Report(100.0);

            Console.WriteLine(@"Done.");
            // Stop the overall progress reports.
            progReportTimer.Dispose();
            swAll.Stop();

            var elapsed = swAll.Elapsed;
            var minutes = elapsed.TotalMinutes - ((int)elapsed.TotalHours * 60);
            Console.WriteLine(@"Total elapsed time for search: {0:f1} sec ({1}d {2}h {3:f2} min)", elapsed.TotalSeconds, elapsed.Days, elapsed.Hours, minutes);

            return true;
        }

        private List<DatabaseSearchResultData> RunDatabaseSearch(FastaDatabase searchDb, string outputFilePath, ISequenceFilter ms1Filter, string searchModeString, IProgress<ProgressData> progress)
        {
            var progressData = new ProgressData(progress);
            var sw = new Stopwatch();
            var searchModeStringCap = char.ToUpper(searchModeString[0]) + searchModeString.Substring(1);

            sw.Reset();
            sw.Start();
            Console.WriteLine(@"Reading the {0} database...", searchModeString);
            searchDb.Read();
            sw.Stop();
            Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

            var matches = new SortedSet<DatabaseSequenceSpectrumMatch>[_run.MaxLcScan + 1];
            progressData.StepRange(50);
            if (Options.TagBasedSearch)
            {
                sw.Reset();
                Console.WriteLine(@"Tag-based searching the {0} database", searchModeString);
                sw.Start();
                var progTag = new Progress<ProgressData>(p =>
                {
                    progressData.StatusInternal = p.Status;
                    progressData.Report(p.Percent);
                });
                RunTagBasedSearch(matches, searchDb, null, progTag);
                Console.WriteLine(@"{1} database tag-based search elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds, searchModeStringCap);
            }
            progressData.StepRange(100);

            sw.Reset();
            Console.WriteLine(@"Searching the {0} database", searchModeString);
            sw.Start();
            var prog = new Progress<ProgressData>(p =>
            {
                progressData.StatusInternal = p.Status;
                progressData.Report(p.Percent);
            });
            RunSearch(matches, searchDb, ms1Filter, null, prog);
            Console.WriteLine(@"{1} database search elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds, searchModeStringCap);

            // calculate spectral e-value using generating function
            sw.Reset();
            Console.WriteLine(@"Calculating spectral E-values for {0}-spectrum matches", searchModeString);
            sw.Start();
            var bestMatches = RunGeneratingFunction(matches);
            var results = WriteResultsToFile(bestMatches, outputFilePath, searchDb);
            sw.Stop();
            Console.WriteLine(@"{1}-spectrum match E-value calculation elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds, searchModeStringCap);
            return results;
        }

        private int[] _tagMs2ScanNum;

        private IEnumerable<AnnotationAndOffset> GetAnnotationsAndOffsets(FastaDatabase database, out long estimatedProteins, CancellationToken? cancellationToken = null)
        {
            var indexedDb = new IndexedDatabase(database);
            indexedDb.Read();
            estimatedProteins = indexedDb.EstimateTotalPeptides(Options.InternalCleavageMode, Options.MinSequenceLength, Options.MaxSequenceLength, Options.MaxNumNTermCleavages, Options.MaxNumCTermCleavages);
            IEnumerable<AnnotationAndOffset> annotationsAndOffsets;
            if (Options.InternalCleavageMode == InternalCleavageType.MultipleInternalCleavages)
            {
                //annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(MinSequenceLength, MaxSequenceLength);
                annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzymeParallel(Options.MinSequenceLength, Options.MaxSequenceLength, Options.MaxNumThreads, cancellationToken);
            }
            else if (Options.InternalCleavageMode == InternalCleavageType.NoInternalCleavage)
            {
                annotationsAndOffsets = indexedDb.IntactSequenceAnnotationsAndOffsets(Options.MinSequenceLength, Options.MaxSequenceLength, Options.MaxNumCTermCleavages);
            }
            else
            {
                annotationsAndOffsets = indexedDb
                    .SequenceAnnotationsAndOffsetsWithNtermOrCtermCleavageNoLargerThan(
                        Options.MinSequenceLength, Options.MaxSequenceLength, Options.MaxNumNTermCleavages, Options.MaxNumCTermCleavages);
            }

            return annotationsAndOffsets;
        }

        private void RunTagBasedSearch(SortedSet<DatabaseSequenceSpectrumMatch>[] matches, FastaDatabase db,
                                        CancellationToken? cancellationToken = null, IProgress<ProgressData> progress = null)
        {
            _tagSearchEngine.SetDatabase(db);

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
                MaxDegreeOfParallelism = Options.MaxNumThreads,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            Parallel.ForEach(_tagMs2ScanNum, pfeOptions, ms2ScanNum =>
            {
                var tagSeqMatches = _tagSearchEngine.RunSearch(ms2ScanNum);

                foreach (var tagSequenceMatch in tagSeqMatches)
                {
                    var offset = _tagSearchEngine.FastaDatabase.GetOffset(tagSequenceMatch.ProteinName);
                    if (offset == null) continue;

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
            var annotationsAndOffsets = GetAnnotationsAndOffsets(db, out var estimatedProteins, cancellationToken);
            Console.WriteLine(@"Estimated Sequences: " + estimatedProteins.ToString("#,##0"));

            var numProteins = 0;
            var lastUpdate = DateTime.MinValue; // Force original update of 0%

            sw.Reset();
            sw.Start();
            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Options.MaxNumThreads,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            var maxNumNTermCleavages = Options.InternalCleavageMode == InternalCleavageType.NoInternalCleavage ? Options.MaxNumNTermCleavages : 0;
            //foreach (var annotationAndOffset in annotationsAndOffsets)
            Parallel.ForEach(annotationsAndOffsets, pfeOptions, annotationAndOffset =>
            {
                if (cancellationToken != null && cancellationToken.Value.IsCancellationRequested)
                {
                    //return matches;
                    return;
                }

                SearchProgressReport(ref numProteins, ref lastUpdate, estimatedProteins, sw, progData);
                SearchForMatches(annotationAndOffset, sequenceFilter, matches, maxNumNTermCleavages, db.IsDecoy, cancellationToken);
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

        private void SearchForMatches(AnnotationAndOffset annotationAndOffset,
            ISequenceFilter sequenceFilter, SortedSet<DatabaseSequenceSpectrumMatch>[] matches, int maxNumNTermCleavages, bool isDecoy, CancellationToken? cancellationToken = null)
        {
            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Options.MaxNumThreads,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            var annotation = annotationAndOffset.Annotation;
            var offset = annotationAndOffset.Offset;
            //var protein = db.GetProteinName(offset);
            var protSequence = annotation.Substring(2, annotation.Length - 4);
            var seqGraph = SequenceGraph.CreateGraph(Options.AminoAcidSet, AminoAcid.ProteinNTerm, protSequence,
                AminoAcid.ProteinCTerm);

            if (seqGraph == null) return; // No matches will be found without a sequence graph.

            for (var numNTermCleavages = 0; numNTermCleavages <= maxNumNTermCleavages; numNTermCleavages++)
            {
                if (numNTermCleavages > 0) seqGraph.CleaveNTerm();
                var numProteoforms = seqGraph.GetNumProteoformCompositions();
                var modCombs = seqGraph.GetModificationCombinations();
                for (var modIndex = 0; modIndex < numProteoforms; modIndex++)
                {
                    seqGraph.SetSink(modIndex);
                    var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();
                    var sequenceMass = protCompositionWithH2O.Mass;

                    if (sequenceMass < Options.MinSequenceMass || sequenceMass > Options.MaxSequenceMass) continue;

                    var modCombinations = modCombs[modIndex];
                    var ms2ScanNums = Options.ScanNumbers ?? sequenceFilter.GetMatchingMs2ScanNums(sequenceMass);
                    var filter = sequenceFilter as Ms1FtFilter;
                    var featureIds = filter?.GetMatchingFeatureIds(sequenceMass).ToList() ?? new List<int>();

                    Parallel.ForEach(ms2ScanNums, pfeOptions, ms2ScanNum =>
                    {
                        if (ms2ScanNum > _ms2ScanNums.Last() || ms2ScanNum < _ms2ScanNums.First()) return;

                        var isoTargetMz = _isolationWindowTargetMz[ms2ScanNum];
                        if (!(isoTargetMz > 0)) return;
                        var charge = (int)Math.Round(sequenceMass / (isoTargetMz - Constants.Proton));

                        //var scorer = _ms2ScorerFactory2.GetScorer(ms2ScanNum, sequenceMass, charge, Options.ActivationMethod);
                        var scorer = _ms2ScorerFactory2.GetScorer(ms2ScanNum);
                        var score = seqGraph.GetFragmentScore(scorer);

                        var precursorIon = new Ion(protCompositionWithH2O, charge);
                        var sequence = protSequence.Substring(numNTermCleavages);
                        var pre = numNTermCleavages == 0 ? annotation[0] : annotation[numNTermCleavages + 1];
                        var post = annotation[annotation.Length - 1];

                        var ms2FeatureIds = new List<int>();
                        if (filter != null)
                        {
                            foreach (var featureId in featureIds)
                            {
                                var scanRange = filter.Ms1FtIndexToScanRange[featureId];
                                if (ms2ScanNum >= scanRange.Item1 && ms2ScanNum <= scanRange.Item2) ms2FeatureIds.Add(featureId);
                            }
                        }
                        else ms2FeatureIds = featureIds;
                        var prsm = new DatabaseSequenceSpectrumMatch(sequence, pre, post, ms2ScanNum, offset, numNTermCleavages,
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
                    matches[ms2ScanNum] = new SortedSet<DatabaseSequenceSpectrumMatch> {prsm};
                }
                else // already exists
                {
                    var existingMatches = matches[ms2ScanNum];
                    //var maxScore = existingMatches.Max.Score;
                    if (existingMatches.Count < Options.NumMatchesPerSpectrum)
                    {
                        //if (!(maxScore*0.7 < prsm.Score)) return;
                        existingMatches.Add(prsm);
                    }
                    else
                    {
                        var minScore = existingMatches.Min.Score;
                        if (!(prsm.Score > minScore)) return;
                        existingMatches.Add(prsm);
                        existingMatches.Remove(existingMatches.Min);
                    }
                    //if (NumMatchesPerSpectrum > 1) existingMatches.RemoveWhere(mt => mt.Score < maxScore * 0.7);
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

            // Rescore and Estimate #proteins for GF calculation
            long estimatedProteins = scanNums.Count;
            Console.WriteLine(@"Number of spectra: " + estimatedProteins);
            var numProteins = 0;
            var lastUpdate = DateTime.MinValue; // Force original update of 0%
            sw.Reset();
            sw.Start();

            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Options.MaxNumThreads,
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

        private LinkedList<Tuple<double, ScoreDistribution>>[] _cachedScoreDistributions;
        private DatabaseSequenceSpectrumMatch[] RunGeneratingFunction(IReadOnlyList<SortedSet<DatabaseSequenceSpectrumMatch>> sortedMatches, CancellationToken? cancellationToken = null, IProgress<ProgressData> progress = null)
        {
            var progData = new ProgressData(progress)
            {
                Status = "Calculating spectral E-values for matches"
            };

            if (_cachedScoreDistributions == null)
            {
                _cachedScoreDistributions = new LinkedList<Tuple<double, ScoreDistribution>>[_run.MaxLcScan + 1];
                foreach (var scanNum in _ms2ScanNums)
                    _cachedScoreDistributions[scanNum] = new LinkedList<Tuple<double, ScoreDistribution>>();
            }

            var sw = new Stopwatch();

            var topDownScorer = new InformedTopDownScorer(_run, Options.AminoAcidSet, Options.MinProductIonCharge, Options.MaxProductIonCharge, Options.ProductIonTolerance, activationMethod: Options.ActivationMethod);

            // Rescore and Estimate #proteins for GF calculation
            var matches = new LinkedList<DatabaseSequenceSpectrumMatch>[sortedMatches.Length];
            long estimatedSequences = 0;
            foreach(var scanNum in _ms2ScanNums)
            {
                var prsms = sortedMatches[scanNum];
                if (prsms == null) continue;
                var spec = _run.GetSpectrum(scanNum) as ProductSpectrum;
                if (spec == null || spec.Peaks.Length == 0)
                    continue;

                foreach (var match in prsms)
                {
                    //if (CompositeScorer.GetProbability(match.Score) < CompositeScorer.ProbabilityCutoff)
                    //{
                    //    continue;
                    //}

                    var ion = match.Ion;

                    // Re-scoring
                    //var scores = topDownScorer.GetScores(spec, match.Sequence, ion.Composition, ion.Charge, scanNum);
                    //if (scores == null) continue;
                    var scorer = this._ms2ScorerFactory2.GetScorer(this._run.GetSpectrum(match.ScanNum) as ProductSpectrum, ion.Composition.Mass, ion.Charge, this.Options.ActivationMethod);
                    var informedScorer = scorer as IInformedScorer;
                    var scores = topDownScorer.GetIcScores(informedScorer, scorer, match.Sequence, ion.Composition);

                    match.ModificationText = scores.Modifications;

                    //double s1 = scores.Score, s2;
                    //int nf1;
                    //topDownScorer.GetCompositeScores(sequence, ion.Charge, scanNum, out s1, out nf1);
                    //topDownScorer.GetCompositeScores(ipSequence, this._ms2ScorerFactory2.GetMs2Scorer(scanNum) as CompositeScorerBasedOnDeconvolutedSpectrum, out s2);

                    match.Score = scores.Score;
                    match.NumMatchedFragments = scores.NumMatchedFrags;
                    //match.ModificationText = match.Modifications.ToString();
                    //match.NumMatchedFragments = scores.NumMatchedFrags;
                    if (match.Score > 0.25*informedScorer.ScoreCutOff)
                    //if (match.Score > CompositeScorer.ScoreParam.Cutoff)
                    {
                        if (matches[scanNum] == null) matches[scanNum] = new LinkedList<DatabaseSequenceSpectrumMatch>();
                        matches[scanNum].AddLast(match);
                    }
                }

                if (matches[scanNum] != null) estimatedSequences += matches[scanNum].Count;
            }

            Console.WriteLine(@"Estimated matched sequences: " + estimatedSequences.ToString("#,##0"));

            var numProteins = 0;
            var lastUpdate = DateTime.MinValue; // Force original update of 0%
            sw.Reset();
            sw.Start();

            var scanNums = _ms2ScanNums.Where(scanNum => matches[scanNum] != null).ToArray();

            var pfeOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Options.MaxNumThreads,
                CancellationToken = cancellationToken ?? CancellationToken.None
            };

            Parallel.ForEach(scanNums, pfeOptions, scanNum =>
            //foreach (var scanNum in scanNums)
            {
                var currentTask = "?";
                try
                {
                    var scoreDistributions = _cachedScoreDistributions[scanNum];
                    foreach (var match in matches[scanNum])
                    {
                        var currentIteration = "for scan " + scanNum + " and mass " + match.Ion.Composition.Mass;
                        currentTask = "Calling GetMs2ScoringGraph " + currentIteration;

                        //var scorer = _ms2ScorerFactory2.GetScorer(match.ScanNum);
                        //var graph = this.scoringGraphFactory.GetScoringGraph(scorer, match.Ion.Composition.Mass);
                        var graph = _ms2ScorerFactory2.GetMs2ScoringGraph(scanNum, match.Ion.Composition.Mass);
                        if (graph == null) continue;

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
                        //var ipSequence = Sequence.CreateSequence(match.Sequence, match.ModificationText, topDownScorer.AminoAcidSet);
                        //var graphScore = graph.ScoreSequence(ipSequence);

                        match.SpecEvalue = scoreDist.GetSpectralEValue(match.Score); //* 1e-34;

                        currentTask = "Reporting progress " + currentIteration;
                        SearchProgressReport(ref numProteins, ref lastUpdate, estimatedSequences, sw, progData);
                    }
            }
                catch (Exception ex)
            {
                var errMsg = string.Format("Exception while {0}: {1}", currentTask, ex.Message);
                Console.WriteLine(errMsg);
                throw new Exception(errMsg, ex);
            }
        });

            var finalMatches = new DatabaseSequenceSpectrumMatch[matches.Length];
            foreach (var scanNum in scanNums)
            {
                finalMatches[scanNum] = matches[scanNum].OrderBy(m => m.SpecEvalue).First();
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
                nMatches += _ms2ScanNums.Where(scanNum => matches[scanNum] != null).Sum(scanNum => matches[scanNum].Count);
            }
            return nMatches;
        }

        private List<DatabaseSearchResultData> WriteResultsToFile(IReadOnlyList<DatabaseSequenceSpectrumMatch> matches, string outputFilePath, FastaDatabase database)
        {
            var results = CreateResults(matches, database).ToList();
            DatabaseSearchResultData.WriteResultsToFile(outputFilePath, results);
            return results;
        }

        private IEnumerable<DatabaseSearchResultData> CreateResults(IReadOnlyList<DatabaseSequenceSpectrumMatch> matches, FastaDatabase database)
        {
            foreach (var scanNum in _ms2ScanNums)
            {
                var match = matches[scanNum];
                if (match == null)
                    continue;

                var start = database.GetOneBasedPositionInProtein(match.Offset) + 1 + match.NumNTermCleavages;
                var proteinName = database.GetProteinName(match.Offset);
                var result = new DatabaseSearchResultData()
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
}
