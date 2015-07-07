using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Scoring;

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
            int maxSequenceLength = 300,
            int maxNumNTermCleavages = 1,
            int maxNumCTermCleavages = 0,
            int minPrecursorIonCharge = 2,
            int maxPrecursorIonCharge = 30,
            int minProductIonCharge = 1,
            int maxProductIonCharge = 15,
            double minSequenceMass = 3000.0,
            double maxSequenceMass = 50000.0,
            double precursorIonTolerancePpm = 10,
            double productIonTolerancePpm = 10,
            bool? runTargetDecoyAnalysis = true,
            int searchMode = 1,
            string featureFilePath = null,
            double minFeatureProbability = 0.15,
            IEnumerable<int> scanNumbers = null,
            int numMatchesPerSpectrum = 1
            )
        {
            ErrorMessage = string.Empty;

            SpecFilePath = specFilePath;
            DatabaseFilePath = dbFilePath;
            AminoAcidSet = aaSet;
            FeatureFilePath = featureFilePath;

            MinFeatureProbability = minFeatureProbability;
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
            MaxNumThreads = 0;
            ForceParallel = false;
        }

        public string ErrorMessage { get; private set; }
        public string SpecFilePath { get; private set; }
        public string DatabaseFilePath { get; private set; }
        public string OutputDir { get; private set; }
        public AminoAcidSet AminoAcidSet { get; private set; }
        public string FeatureFilePath { get; private set; }
        public double MinFeatureProbability { get; private set; }
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
        public bool? RunTargetDecoyAnalysis { get; private set; } // true: target and decoy, false: target only, null: decoy only
        public IEnumerable<int> ScanNumbers { get; private set; }
        public int NumMatchesPerSpectrum { get; private set; }
        public int MaxNumThreads { get; set; }
        public bool ForceParallel { get; set; }

        // 0: all internal sequences, 
        // 1: #NCleavages <= Max OR Cleavages <= Max (Default)
        // 2: 1: #NCleavages <= Max AND Cleavages <= Max
        public int SearchMode { get; private set; } 

        private LcMsRun _run;
        private ProductScorerBasedOnDeconvolutedSpectra _ms2ScorerFactory;
        private InformedTopDownScorer _topDownScorer;

        public bool RunSearch(double corrThreshold = 0.7, CancellationToken? cancellationToken=null, IProgress<ProgressData> progress = null)
        {
            Progress<ProgressData> prog = new Progress<ProgressData>();
            var progData = new ProgressData();
            if (progress != null)
            {
                prog = new Progress<ProgressData>(p =>
                {
                    progData.Status = p.Status;
                    progress.Report(progData.UpdatePercent(p.Percent));
                });
            }
            else
            {
                progress = new Progress<ProgressData>();
            }

            var sw = new Stopwatch();
            var swAll = new Stopwatch();
            swAll.Start();
            ErrorMessage = string.Empty;

            Console.Write(@"Reading raw file...");
            progData.Status = "Reading spectra file";
            progData.IsPartialRange = true;
            progData.MaxPercentage = 10.0;
            sw.Start();
            //_run = InMemoryLcMsRun.GetLcMsRun(SpecFilePath, MassSpecDataType.XCaliburRun, 0, 0);   // 1.4826
            _run = PbfLcMsRun.GetLcMsRun(SpecFilePath, 0, 0, prog);
            _topDownScorer = new InformedTopDownScorer(_run, AminoAcidSet, MinProductIonCharge, MaxProductIonCharge, ProductIonTolerance, corrThreshold);
            sw.Stop();
            Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);


            //var sequenceFilter = new Ms1IsotopeAndChargeCorrFilter(_run, PrecursorIonTolerance, MinPrecursorIonCharge,
            //    MaxPrecursorIonCharge,
            //    MinSequenceMass, MaxSequenceMass, corrThreshold, 0.2, 0.2); //corrThreshold, corrThreshold);

            progData.StepRange(20.0);
            ISequenceFilter ms1Filter;
            if (string.IsNullOrWhiteSpace(FeatureFilePath))
            {
                // Checks whether SpecFileName.ms1ft exists
                var ms1FtFilePath = MassSpecDataReaderFactory.ChangeExtension(SpecFilePath, Ms1FeatureFinderLauncher.FileExtension);
                if (!File.Exists(ms1FtFilePath))
                {
                    Console.Write(@"Running ProMex...");
                    sw.Start();
                    var param = new Ms1FeatureFinderInputParameter
                    {
                        InputPath = SpecFilePath,
                        MinSearchCharge = MinPrecursorIonCharge,
                        MaxSearchCharge = MaxPrecursorIonCharge
                    };
                    var featureFinder = new Ms1FeatureFinderLauncher(param);
                    featureFinder.Run();
                    //var extractor = new Ms1FeatureMatrix(_run, MinPrecursorIonCharge, MaxPrecursorIonCharge);
                    //ms1FtFilePath = extractor.GetFeatureFile(SpecFilePath, MinSequenceMass, MaxSequenceMass);
                }
                sw.Reset();
                sw.Start();
                Console.Write(@"Reading ProMex results...");
                ms1Filter = new Ms1FtFilter(_run, PrecursorIonTolerance, ms1FtFilePath, MinFeatureProbability);
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
                    ms1Filter = new Ms1FtFilter(_run, PrecursorIonTolerance, FeatureFilePath, MinFeatureProbability);
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

            _ms2ScorerFactory = new ProductScorerBasedOnDeconvolutedSpectra(
                _run,
                MinProductIonCharge, MaxProductIonCharge,
                ProductIonTolerance
                );

            progData.StepRange(25.0);
            progData.Status = "Reading Fasta File";
            progress.Report(progData.UpdatePercent(100.0)); // Output 25.0%

            // Target database
            var targetDb = new FastaDatabase(DatabaseFilePath);
            targetDb.Read();

//            string dirName = OutputDir ?? Path.GetDirectoryName(SpecFilePath);
            var specFileName = MassSpecDataReaderFactory.RemoveExtension(Path.GetFileName(SpecFilePath));
            var targetOutputFilePath = Path.Combine(OutputDir, specFileName + TargetFileNameEnding);
            var decoyOutputFilePath = Path.Combine(OutputDir, specFileName + DecoyFileNameEnding);
            var tdaOutputFilePath = Path.Combine(OutputDir, specFileName + TdaFileNameEnding);

            progData.StepRange(60.0);
            progData.Status = "Running Target search";
            progress.Report(progData.UpdatePercent(0.0));
            if (RunTargetDecoyAnalysis != null)
            {
                sw.Reset();
                Console.Write(@"Reading the target database...");
                sw.Start();
                targetDb.Read();
                sw.Stop();
                Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

                sw.Reset();
                Console.WriteLine(@"Searching the target database");
                sw.Start();
                SortedSet<DatabaseSequenceSpectrumMatch>[] targetMatches;
                if (ForceParallel || (SearchMode == 0 && MaxNumThreads != 1))
                {
                    targetMatches = RunSearchParallel(targetDb, ms1Filter, null, prog);
                }
                else
                {
                    targetMatches = RunSearch(targetDb, ms1Filter, null, prog);
                }
                WriteResultsToFile(targetMatches, targetOutputFilePath, targetDb);
                sw.Stop();

                Console.WriteLine(@"Target database search elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);
            }

            progData.StepRange(95.0);
            progData.Status = "Running Decoy search";
            progress.Report(progData.UpdatePercent(0.0));
            if (RunTargetDecoyAnalysis == true || RunTargetDecoyAnalysis == null)
            {
                // Decoy database
                sw.Reset();
                Console.Write(@"Reading the decoy database...");
                sw.Start();
                var decoyDb = targetDb.Decoy(null, true);
                decoyDb.Read();

                Console.WriteLine(@"Elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);

                sw.Reset();
                Console.WriteLine(@"Searching the decoy database");
                sw.Start();
                SortedSet<DatabaseSequenceSpectrumMatch>[] decoyMatches;
                if (ForceParallel || (SearchMode == 0 && MaxNumThreads != 1))
                {
                    decoyMatches = RunSearchParallel(decoyDb, ms1Filter, null, prog);
                }
                else
                {
                    decoyMatches = RunSearch(decoyDb, ms1Filter, null, prog);
                }
                WriteResultsToFile(decoyMatches, decoyOutputFilePath, decoyDb);
                sw.Stop();

                Console.WriteLine(@"Decoy database search elapsed Time: {0:f1} sec", sw.Elapsed.TotalSeconds);
            }

            progData.StepRange(100.0);
            progData.Status = "Writing combined results file";
            progress.Report(progData.UpdatePercent(0.0));
            if (RunTargetDecoyAnalysis == true)
            {
                var fdrCalculator = new FdrCalculator(targetOutputFilePath, decoyOutputFilePath);
                if (fdrCalculator.HasError())
                {
                    ErrorMessage = fdrCalculator.ErrorMessage;
                    Console.WriteLine(@"Error computing FDR: " + fdrCalculator.ErrorMessage);
                    return false;
                }
                
                fdrCalculator.WriteTo(tdaOutputFilePath);
            }
            progress.Report(progData.UpdatePercent(100.0));

            Console.WriteLine(@"Done.");
            swAll.Stop();
            Console.WriteLine(@"Total elapsed time for search: {0:f1} sec ({1:f2} min)", swAll.Elapsed.TotalSeconds, swAll.Elapsed.TotalMinutes);

            return true;
        }

        private IEnumerable<AnnotationAndOffset> GetAnnotationsAndOffsets(FastaDatabase database, out long estimatedProteins, CancellationToken? cancellationToken = null)
        {
            var indexedDb = new IndexedDatabase(database);
            indexedDb.Read();
            estimatedProteins = indexedDb.EstimateTotalPeptides(SearchMode, MinSequenceLength, MaxSequenceLength, MaxNumNTermCleavages, MaxNumCTermCleavages);
            IEnumerable<AnnotationAndOffset> annotationsAndOffsets;
            if (SearchMode == 0)
            {
                if (ForceParallel || (SearchMode == 0 && MaxNumThreads != 1))
                {
                    annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzymeParallel(MinSequenceLength, MaxSequenceLength, MaxNumThreads, cancellationToken);
                }
                else
                {
                    annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(MinSequenceLength, MaxSequenceLength);
                }
            }
            else if (SearchMode == 2)
            {
                annotationsAndOffsets = indexedDb.IntactSequenceAnnotationsAndOffsets(MinSequenceLength,
                    MaxSequenceLength, MaxNumCTermCleavages);
            }
            else
            {
                annotationsAndOffsets = indexedDb
                    .SequenceAnnotationsAndOffsetsWithNtermOrCtermCleavageNoLargerThan(
                        MinSequenceLength, MaxSequenceLength, MaxNumNTermCleavages, MaxNumCTermCleavages);
            }

            return annotationsAndOffsets;
        }

        private SortedSet<DatabaseSequenceSpectrumMatch>[] RunSearch(FastaDatabase db, ISequenceFilter sequenceFilter, CancellationToken? cancellationToken=null, IProgress<ProgressData> progress = null)
        {
            if (progress == null)
            {
                progress = new Progress<ProgressData>();
            }
            var progData = new ProgressData();
            progData.Status = "Searching for matches";
            var sw = new Stopwatch();
            long estimatedProteins;
            var annotationsAndOffsets = GetAnnotationsAndOffsets(db, out estimatedProteins);
            Console.WriteLine(@"Estimated proteins: " + estimatedProteins);
            
            var numProteins = 0;
            var lastUpdate = DateTime.UtcNow;

            sw.Reset();
            sw.Start();

            var matches = new SortedSet<DatabaseSequenceSpectrumMatch>[_run.MaxLcScan+1];

            var maxNumNTermCleavages = SearchMode == 2 ? MaxNumNTermCleavages : 0;

            foreach (var annotationAndOffset in annotationsAndOffsets)
            {
                if (cancellationToken != null && cancellationToken.Value.IsCancellationRequested)
                {
                    return matches;
                }

                var annotation = annotationAndOffset.Annotation;
                var offset = annotationAndOffset.Offset;

                var protein = db.GetProteinName(offset);

                progress.Report(progData.UpdatePercent((double)numProteins / (double)estimatedProteins * 100.0));

                //++numProteins;
                Interlocked.Increment(ref numProteins);

                if (DateTime.UtcNow.Subtract(lastUpdate).TotalSeconds >= 15)
                {
                    lastUpdate = DateTime.UtcNow;

                    Console.WriteLine(@"Processing, {0} proteins done, {1:#0.0}% complete, {2:f1} sec elapsed",
                        numProteins,
                        numProteins / (double)estimatedProteins * 100.0,
                        sw.Elapsed.TotalSeconds);           
                }

                var protSequence = annotation.Substring(2, annotation.Length - 4);

                var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, AminoAcid.ProteinNTerm, protSequence,
                    AminoAcid.ProteinCTerm);
                if (seqGraph == null) continue;

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
                        var modCombinations = modCombs[modIndex];

                        var ms2ScanNums = this.ScanNumbers ?? sequenceFilter.GetMatchingMs2ScanNums(sequenceMass);

                        foreach (var ms2ScanNum in ms2ScanNums)
                        {
                            if (ms2ScanNum > _run.MaxLcScan) continue;

                            var spec = _run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                            if (spec == null) continue;
                            var charge =
                                (int)
                                    Math.Round(sequenceMass/
                                               (spec.IsolationWindow.IsolationWindowTargetMz - Constants.Proton));
                            var scorer = _ms2ScorerFactory.GetMs2Scorer(ms2ScanNum);
                            var score = seqGraph.GetFragmentScore(scorer);
                            if (score <= 3) continue;

                            var precursorIon = new Ion(protCompositionWithH2O, charge);
                            var sequence = protSequence.Substring(numNTermCleavages);
                            var pre = numNTermCleavages == 0 ? annotation[0] : annotation[numNTermCleavages + 1];
                            var post = annotation[annotation.Length - 1];

                            var prsm = new DatabaseSequenceSpectrumMatch(sequence, pre, post, ms2ScanNum, offset,
                                numNTermCleavages,
                                modCombinations, precursorIon, score);

                            if (matches[ms2ScanNum] == null)
                            {
                                matches[ms2ScanNum] = new SortedSet<DatabaseSequenceSpectrumMatch> {prsm};
                            }
                            else // already exists
                            {
                                var existingMatches = matches[ms2ScanNum];
                                if (existingMatches.Count < NumMatchesPerSpectrum) existingMatches.Add(prsm);
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
            }
            progress.Report(progData.UpdatePercent(100.0));
            return matches;
        }

        private SortedSet<DatabaseSequenceSpectrumMatch>[] RunSearchParallel(FastaDatabase db, ISequenceFilter sequenceFilter, CancellationToken? cancellationToken = null, IProgress<ProgressData> progress = null)
        {
            if (progress == null)
            {
                progress = new Progress<ProgressData>();
            }
            var progData = new ProgressData();
            progData.Status = "Searching for matches";

            var threads = MaxNumThreads;
            var sw = new Stopwatch();

            // Try to get the number of physical cores in the system - requires System.Management.dll and a WMI query, but the performance penalty for 
            // using the number of logical processors in a hyperthreaded system is significant, and worse than the penalty for using fewer than all physical cores.
            int coreCount = 0;
            try
            {
                foreach (var item in new System.Management.ManagementObjectSearcher("Select NumberOfCores from Win32_Processor").Get())
                {
                    coreCount += int.Parse(item["NumberOfCores"].ToString());
                }
                //Console.WriteLine(@"Number Of Cores: {0}", coreCount);
            }
            catch (Exception)
            {
                // Use the logical processor count, divided by 2 to avoid the greater performance penalty of over-threading.
                coreCount = (int)(Math.Ceiling(System.Environment.ProcessorCount / 2.0));
            }

            if (threads <= 0 || threads > coreCount)
            {
                threads = coreCount;
            }

            long estimatedProteins;
            var annotationsAndOffsets = GetAnnotationsAndOffsets(db, out estimatedProteins, cancellationToken);
            Console.WriteLine(@"Estimated proteins: " + estimatedProteins);
            
            var numProteins = 0;
            var lastUpdate = DateTime.UtcNow;

            sw.Reset();
            sw.Start();

            var matches = new SortedSet<DatabaseSequenceSpectrumMatch>[_run.MaxLcScan + 1];

            var maxNumNTermCleavages = SearchMode == 2 ? MaxNumNTermCleavages : 0;

            var pfeOptions = new ParallelOptions();
            pfeOptions.MaxDegreeOfParallelism = threads;
            pfeOptions.CancellationToken = cancellationToken != null ? (CancellationToken)cancellationToken : CancellationToken.None;

            //foreach (var annotationAndOffset in annotationsAndOffsets)
            Parallel.ForEach(annotationsAndOffsets, pfeOptions, annotationAndOffset =>
            {
                var annotation = annotationAndOffset.Annotation;
                var offset = annotationAndOffset.Offset;

                var protein = db.GetProteinName(offset);

                var tempNumProteins = Interlocked.Increment(ref numProteins);
                //lock (progress)
                //{
                    progress.Report(progData.UpdatePercent((double)(tempNumProteins - 1) / (double)estimatedProteins * 100.0));
                //}

                if (DateTime.UtcNow.Subtract(lastUpdate).TotalSeconds >= 15)
                {
                    lastUpdate = DateTime.UtcNow;

                    Console.WriteLine(@"Processing, {0} proteins done, {1:#0.0}% complete, {2:f1} sec elapsed",
                        tempNumProteins,
                        tempNumProteins / (double)estimatedProteins * 100.0,
                        sw.Elapsed.TotalSeconds);                    
                }

                var protSequence = annotation.Substring(2, annotation.Length - 4);

                var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, AminoAcid.ProteinNTerm, protSequence,
                    AminoAcid.ProteinCTerm);
                if (seqGraph == null) return; // Early exit from this iteration, not for the function.

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
                        var modCombinations = modCombs[modIndex];

                        var ms2ScanNums = this.ScanNumbers ?? sequenceFilter.GetMatchingMs2ScanNums(sequenceMass);

                        foreach (var ms2ScanNum in ms2ScanNums)
                        {
                            if (ms2ScanNum > _run.MaxLcScan) continue;

                            var spec = _run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                            if (spec == null) continue;
                            var charge =
                                (int)
                                    Math.Round(sequenceMass /
                                               (spec.IsolationWindow.IsolationWindowTargetMz - Constants.Proton));
                            var scorer = _ms2ScorerFactory.GetMs2Scorer(ms2ScanNum);
                            var score = seqGraph.GetFragmentScore(scorer);
                            if (score <= 3) continue;

                            var precursorIon = new Ion(protCompositionWithH2O, charge);
                            var sequence = protSequence.Substring(numNTermCleavages);
                            var pre = numNTermCleavages == 0 ? annotation[0] : annotation[numNTermCleavages + 1];
                            var post = annotation[annotation.Length - 1];

                            var prsm = new DatabaseSequenceSpectrumMatch(sequence, pre, post, ms2ScanNum, offset,
                                numNTermCleavages,
                                modCombinations, precursorIon, score);

                            lock (matches)
                            {
                                if (matches[ms2ScanNum] == null)
                                {
                                    matches[ms2ScanNum] = new SortedSet<DatabaseSequenceSpectrumMatch> {prsm};
                                }
                                else // already exists
                                {
                                    var existingMatches = matches[ms2ScanNum];
                                    if (existingMatches.Count < NumMatchesPerSpectrum) existingMatches.Add(prsm);
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
                }
            });
            progress.Report(progData.UpdatePercent(100.0));
            return matches;
        }

        private void WriteResultsToFile(SortedSet<DatabaseSequenceSpectrumMatch>[] matches, string outputFilePath, FastaDatabase database)
        {
            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine("Scan\tPre\tSequence\tPost\tModifications\tComposition\tProteinName\tProteinDesc" +
                             "\tProteinLength\tStart\tEnd\tCharge\tMostAbundantIsotopeMz\tMass\t#MatchedFragments"
                             );
                for (var scanNum = _run.MinLcScan; scanNum <= _run.MaxLcScan; scanNum++)
                {
                    if (matches[scanNum] == null) continue;
                    foreach (var match in matches[scanNum].Reverse())
                    {
                        var sequence = match.Sequence;
                        var offset = match.Offset;
                        var start = database.GetOneBasedPositionInProtein(offset) + 1 + match.NumNTermCleavages;
                        var end = start + sequence.Length - 1;
                        var proteinName = database.GetProteinName(match.Offset);
                        var protLength = database.GetProteinLength(proteinName);
                        var ion = match.Ion;

                        // Re-scoring
                        var scores = _topDownScorer.GetScores(AminoAcid.ProteinNTerm, sequence, AminoAcid.ProteinCTerm, ion.Composition, ion.Charge, scanNum);

                        writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}",
                            scanNum,
                            match.Pre,  // Pre
                            sequence, // Sequence
                            match.Post, // Post
                            scores.Modifications, // Modifications
                            ion.Composition, // Composition
                            proteinName, // ProteinName
                            database.GetProteinDescription(match.Offset), // ProteinDescription
                            protLength, // ProteinLength
                            start, // Start
                            end, // End
                            ion.Charge, // precursorCharge
                            ion.GetMostAbundantIsotopeMz(), // MostAbundantIsotopeMz
                            ion.Composition.Mass,   // Mass
                            scores.Ms2Score    // Score (re-scored)
                            );
                    }
                }
            }
        }
    }
}
