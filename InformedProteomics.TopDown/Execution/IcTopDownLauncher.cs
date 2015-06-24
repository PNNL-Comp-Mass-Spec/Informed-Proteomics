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
        public const string TargetFileExtension = "_IcTarget.tsv";
        public const string DecoyFileExtension = "_IcDecoy.tsv";
        public const string TdaFileExtension = "_IcTda.tsv";

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
        }

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

        // 0: all internal sequences, 
        // 1: #NCleavages <= Max OR Cleavages <= Max (Default)
        // 2: 1: #NCleavages <= Max AND Cleavages <= Max
        public int SearchMode { get; private set; } 

        private LcMsRun _run;
        private ProductScorerBasedOnDeconvolutedSpectra _ms2ScorerFactory;
        private InformedTopDownScorer _topDownScorer;

        public void RunSearch(double corrThreshold = 0.7, CancellationToken? cancellationToken=null)
        {
            var sw = new Stopwatch();

            Console.Write("Reading raw file...");
            sw.Start();
            //_run = InMemoryLcMsRun.GetLcMsRun(SpecFilePath, MassSpecDataType.XCaliburRun, 0, 0);   // 1.4826
            _run = PbfLcMsRun.GetLcMsRun(SpecFilePath, MassSpecDataType.XCaliburRun, 0, 0);
            _topDownScorer = new InformedTopDownScorer(_run, AminoAcidSet, MinProductIonCharge, MaxProductIonCharge, ProductIonTolerance, corrThreshold);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);


            //var sequenceFilter = new Ms1IsotopeAndChargeCorrFilter(_run, PrecursorIonTolerance, MinPrecursorIonCharge,
            //    MaxPrecursorIonCharge,
            //    MinSequenceMass, MaxSequenceMass, corrThreshold, 0.2, 0.2); //corrThreshold, corrThreshold);

            ISequenceFilter ms1Filter;
            if (FeatureFilePath == null)
            {
                // Checks whether SpecFileName.ms1ft exists
                var ms1FtFilePath = Path.ChangeExtension(SpecFilePath, Ms1FeatureFinderLauncher.FileExtension);
                if (!File.Exists(ms1FtFilePath))
                {
                    Console.Write("Running ProMex...");
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
                Console.Write("Reading ProMex results...");
                ms1Filter = new Ms1FtFilter(_run, PrecursorIonTolerance, ms1FtFilePath, MinFeatureProbability);
            }
            else
            {
                sw.Reset();
                sw.Start();
                var extension = Path.GetExtension(FeatureFilePath);
                if (extension.ToLower().Equals(".csv"))
                {
                    Console.Write("Reading ICR2LS/Decon2LS results...");
                    ms1Filter = new IsosFilter(_run, PrecursorIonTolerance, FeatureFilePath);
                }
                else if (extension.ToLower().Equals(".ms1ft"))
                {
                    Console.Write("Reading ProMex results...");
                    ms1Filter = new Ms1FtFilter(_run, PrecursorIonTolerance, FeatureFilePath, MinFeatureProbability);
                }
                else if (extension.ToLower().Equals(".msalign"))
                {
                    Console.Write("Reading MS-Align+ results...");
                    ms1Filter = new MsDeconvFilter(_run, PrecursorIonTolerance, FeatureFilePath);
                }
                else ms1Filter = null; //new Ms1FeatureMatrix(_run);
            }

            sw.Stop();
            sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

            _ms2ScorerFactory = new ProductScorerBasedOnDeconvolutedSpectra(
                _run,
                MinProductIonCharge, MaxProductIonCharge,
                ProductIonTolerance
                );

            // Target database
            var targetDb = new FastaDatabase(DatabaseFilePath);
            targetDb.Read();

//            string dirName = OutputDir ?? Path.GetDirectoryName(SpecFilePath);

            var targetOutputFilePath = OutputDir + Path.DirectorySeparatorChar +
                                       Path.GetFileNameWithoutExtension(SpecFilePath) + TargetFileExtension;
            var decoyOutputFilePath = OutputDir + Path.DirectorySeparatorChar +
                                       Path.GetFileNameWithoutExtension(SpecFilePath) + DecoyFileExtension;
            var tdaOutputFilePath = OutputDir + Path.DirectorySeparatorChar +
                                    Path.GetFileNameWithoutExtension(SpecFilePath) + TdaFileExtension;

            if (RunTargetDecoyAnalysis != null)
            {
                sw.Reset();
                Console.Write("Reading the target database...");
                sw.Start();
                targetDb.Read();
                sw.Stop();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

                sw.Reset();
                Console.WriteLine("Searching the target database");
                sw.Start();
                var targetMatches = RunSearch(targetDb, ms1Filter);
                WriteResultsToFile(targetMatches, targetOutputFilePath, targetDb);
                sw.Stop();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine(@"Target database search elapsed Time: {0:f4} sec", sec);                
            }

            if (RunTargetDecoyAnalysis == true || RunTargetDecoyAnalysis == null)
            {
                // Decoy database
                sw.Reset();
                Console.Write("Reading the decoy database...");
                sw.Start();
                var decoyDb = targetDb.Decoy(null, true);
                decoyDb.Read();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

                sw.Reset();
                Console.WriteLine("Searching the decoy database");
                sw.Start();
                var decoyMatches = RunSearch(decoyDb, ms1Filter);
                WriteResultsToFile(decoyMatches, decoyOutputFilePath, decoyDb);
                sw.Stop();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine(@"Decoy database search elapsed Time: {0:f4} sec", sec);
            }

            if (RunTargetDecoyAnalysis == true)
            {
                var fdrCalculator = new FdrCalculator(targetOutputFilePath, decoyOutputFilePath);
                fdrCalculator.WriteTo(tdaOutputFilePath);
            }

            Console.WriteLine("Done.");
        }

        public void RunSearchParallel(double corrThreshold = 0.7, CancellationToken? cancellationToken = null)
        {
            var sw = new Stopwatch();

            Console.Write("Reading raw file...");
            sw.Start();
            //_run = InMemoryLcMsRun.GetLcMsRun(SpecFilePath, MassSpecDataType.XCaliburRun, 0, 0);   // 1.4826
            _run = PbfLcMsRun.GetLcMsRun(SpecFilePath, MassSpecDataType.XCaliburRun, 0, 0);
            _topDownScorer = new InformedTopDownScorer(_run, AminoAcidSet, MinProductIonCharge, MaxProductIonCharge, ProductIonTolerance, corrThreshold);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);


            //var sequenceFilter = new Ms1IsotopeAndChargeCorrFilter(_run, PrecursorIonTolerance, MinPrecursorIonCharge,
            //    MaxPrecursorIonCharge,
            //    MinSequenceMass, MaxSequenceMass, corrThreshold, 0.2, 0.2); //corrThreshold, corrThreshold);

            ISequenceFilter ms1Filter;
            if (FeatureFilePath == null)
            {
                // Checks whether SpecFileName.ms1ft exists
                var ms1FtFilePath = Path.ChangeExtension(SpecFilePath, Ms1FeatureFinderLauncher.FileExtension);
                if (!File.Exists(ms1FtFilePath))
                {
                    Console.Write("Running ProMex...");
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
                Console.Write("Reading ProMex results...");
                ms1Filter = new Ms1FtFilter(_run, PrecursorIonTolerance, ms1FtFilePath, MinFeatureProbability);
            }
            else
            {
                sw.Reset();
                sw.Start();
                var extension = Path.GetExtension(FeatureFilePath);
                if (extension.ToLower().Equals(".csv"))
                {
                    Console.Write("Reading ICR2LS/Decon2LS results...");
                    ms1Filter = new IsosFilter(_run, PrecursorIonTolerance, FeatureFilePath);
                }
                else if (extension.ToLower().Equals(".ms1ft"))
                {
                    Console.Write("Reading ProMex results...");
                    ms1Filter = new Ms1FtFilter(_run, PrecursorIonTolerance, FeatureFilePath, MinFeatureProbability);
                }
                else if (extension.ToLower().Equals(".msalign"))
                {
                    Console.Write("Reading MS-Align+ results...");
                    ms1Filter = new MsDeconvFilter(_run, PrecursorIonTolerance, FeatureFilePath);
                }
                else ms1Filter = null; //new Ms1FeatureMatrix(_run);
            }

            sw.Stop();
            sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

            _ms2ScorerFactory = new ProductScorerBasedOnDeconvolutedSpectra(
                _run,
                MinProductIonCharge, MaxProductIonCharge,
                ProductIonTolerance
                );

            // Target database
            var targetDb = new FastaDatabase(DatabaseFilePath);
            targetDb.Read();

            //            string dirName = OutputDir ?? Path.GetDirectoryName(SpecFilePath);

            var targetOutputFilePath = OutputDir + Path.DirectorySeparatorChar +
                                       Path.GetFileNameWithoutExtension(SpecFilePath) + TargetFileExtension;
            var decoyOutputFilePath = OutputDir + Path.DirectorySeparatorChar +
                                       Path.GetFileNameWithoutExtension(SpecFilePath) + DecoyFileExtension;
            var tdaOutputFilePath = OutputDir + Path.DirectorySeparatorChar +
                                    Path.GetFileNameWithoutExtension(SpecFilePath) + TdaFileExtension;

            if (RunTargetDecoyAnalysis != null)
            {
                sw.Reset();
                Console.Write("Reading the target database...");
                sw.Start();
                targetDb.Read();
                sw.Stop();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

                sw.Reset();
                Console.WriteLine("Searching the target database");
                sw.Start();
                var targetMatches = RunSearchParallel(targetDb, ms1Filter);
                WriteResultsToFile(targetMatches, targetOutputFilePath, targetDb);
                sw.Stop();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine(@"Target database search elapsed Time: {0:f4} sec", sec);
            }

            if (RunTargetDecoyAnalysis == true || RunTargetDecoyAnalysis == null)
            {
                // Decoy database
                sw.Reset();
                Console.Write("Reading the decoy database...");
                sw.Start();
                var decoyDb = targetDb.Decoy(null, true);
                decoyDb.Read();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

                sw.Reset();
                Console.WriteLine("Searching the decoy database");
                sw.Start();
                var decoyMatches = RunSearchParallel(decoyDb, ms1Filter);
                WriteResultsToFile(decoyMatches, decoyOutputFilePath, decoyDb);
                sw.Stop();
                sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
                Console.WriteLine(@"Decoy database search elapsed Time: {0:f4} sec", sec);
            }

            if (RunTargetDecoyAnalysis == true)
            {
                var fdrCalculator = new FdrCalculator(targetOutputFilePath, decoyOutputFilePath);
                fdrCalculator.WriteTo(tdaOutputFilePath);
            }

            Console.WriteLine("Done.");
        }

        private IEnumerable<AnnotationAndOffset> GetAnnotationsAndOffsets(FastaDatabase database, out long estimatedProteins)
        {
            var indexedDb = new IndexedDatabase(database);
            indexedDb.Read();
            estimatedProteins = indexedDb.EstimateTotalPeptides(MinSequenceLength, MaxSequenceLength, MaxNumNTermCleavages, MaxNumCTermCleavages);
            IEnumerable<AnnotationAndOffset> annotationsAndOffsets;
            if (SearchMode == 0)
            {
                annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(MinSequenceLength, MaxSequenceLength);
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

        private SortedSet<DatabaseSequenceSpectrumMatch>[] RunSearch(FastaDatabase db, ISequenceFilter sequenceFilter, CancellationToken? cancellationToken=null)
        {
            var sw = new Stopwatch();
            long estimatedProteins;
            var annotationsAndOffsets = GetAnnotationsAndOffsets(db, out estimatedProteins);
            Console.WriteLine("Estimated proteins: " + estimatedProteins);
            var numProteins = 0;
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

                //++numProteins;
                Interlocked.Increment(ref numProteins);

                if (numProteins%100000 == 0)
                //if(numProteins % 10 == 0)
                {
                    Console.Write("Processing {0}{1} proteins..., {3:00.0}%", numProteins,
                        numProteins == 1 ? "st" : numProteins == 2 ? "nd" : numProteins == 3 ? "rd" : "th", (double)numProteins / (double)estimatedProteins);
                    if (numProteins != 0)
                    {
                        sw.Stop();
                        var sec = sw.ElapsedTicks/(double) Stopwatch.Frequency;
                        Console.WriteLine("Elapsed Time: {0:f4} sec", sec);
                        sw.Reset();
                        sw.Start();
                    }
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
            return matches;
        }

        private IEnumerable<AnnotationAndOffset> GetAnnotationsAndOffsetsParallel(FastaDatabase database, out long estimatedProteins, int threads = 0, CancellationToken? cancellationToken = null)
        {
            var indexedDb = new IndexedDatabase(database);
            indexedDb.Read();
            estimatedProteins = indexedDb.EstimateTotalPeptides(MinSequenceLength, MaxSequenceLength, MaxNumNTermCleavages, MaxNumCTermCleavages);
            IEnumerable<AnnotationAndOffset> annotationsAndOffsets;
            if (SearchMode == 0)
            {
                annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzymeParallel(MinSequenceLength, MaxSequenceLength, threads, cancellationToken);
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

        private SortedSet<DatabaseSequenceSpectrumMatch>[] RunSearchParallel(FastaDatabase db, ISequenceFilter sequenceFilter, int threads = 0, CancellationToken? cancellationToken = null)
        {
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
                //Console.WriteLine("Number Of Cores: {0}", coreCount);
            }
            catch (Exception)
            {
                // Use the logical processor count, divided by 2 to avoid the greater performance penalty of over-threading.
                coreCount = (int)(Math.Ceiling(System.Environment.ProcessorCount / 2.0));
            }

            if (threads == 0 || threads > coreCount)
            {
                threads = coreCount;
            }

            long estimatedProteins;
            var annotationsAndOffsets = GetAnnotationsAndOffsetsParallel(db, out estimatedProteins, threads, cancellationToken);
            Console.WriteLine("Estimated proteins: " + estimatedProteins);
            var numProteins = 0;
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

                //lock (numProteins)
                //{
                    //++numProteins;
                    Interlocked.Increment(ref numProteins);

                    if (numProteins % 100000 == 0)
                        //if(numProteins % 10 == 0)
                    {
                        Console.Write("Processing {0}{1} proteins..., {3:00.0}%", numProteins,
                            numProteins == 1 ? "st" : numProteins == 2 ? "nd" : numProteins == 3 ? "rd" : "th", (double)numProteins / (double)estimatedProteins);
                        if (numProteins != 0)
                        {
                            lock (sw)
                            {
                                sw.Stop();
                                var sec = sw.ElapsedTicks / (double) Stopwatch.Frequency;
                                Console.WriteLine("Elapsed Time: {0:f4} sec", sec);
                                sw.Reset();
                                sw.Start();
                            }
                        }
                    }
                //}

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
