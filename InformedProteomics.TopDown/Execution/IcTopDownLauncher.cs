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
    public class IcTopDownLauncher
    {
        public const int NumMatchesPerSpectrum = 1;
        public const string TargetFileExtension = ".ictresult";
        public const string DecoyFileExtension = ".icdresult";

        // Consider all subsequenes of lengths [minSequenceLength,maxSequenceLength]
        public IcTopDownLauncher(
            string specFilePath,
            string dbFilePath,
            AminoAcidSet aaSet,
            int minSequenceLength = 30,
            int maxSequenceLength = 250,
            int minPrecursorIonCharge = 3,
            int maxPrecursorIonCharge = 30,
            int minProductIonCharge = 1,
            int maxProductIonCharge = 10,
            double precursorIonTolerancePpm = 10,
            double productIonTolerancePpm = 15,
            bool runTargetDecoyAnalysis = true)
        {
            SpecFilePath = specFilePath;
            DatabaseFilePath = dbFilePath;
            AminoAcidSet = aaSet;
            MinSequenceLength = minSequenceLength;
            MaxSequenceLength = maxSequenceLength;
            MaxNumNTermCleavages = null;
            MaxNumCTermCleavages = null;
            MinPrecursorIonCharge = minPrecursorIonCharge;
            MaxPrecursorIonCharge = maxPrecursorIonCharge;
            MinProductIonCharge = minProductIonCharge;
            MaxProductIonCharge = maxProductIonCharge;
            PrecursorIonTolerance = new Tolerance(precursorIonTolerancePpm);
            ProductIonTolerance = new Tolerance(productIonTolerancePpm);
            RunTargetDecoyAnalysis = runTargetDecoyAnalysis;
        }

        // Consider intact sequences with N- and C-terminal cleavages
        public IcTopDownLauncher(
            string specFilePath,
            string dbFilePath, 
            AminoAcidSet aaSet,
            int minSequenceLength = 30, 
            int maxSequenceLength = 250, 
            int maxNumNTermCleavages = 30,
            int maxNumCTermCleavages = 0,
            int minPrecursorIonCharge = 3, 
            int maxPrecursorIonCharge = 30,
            int minProductIonCharge = 1, 
            int maxProductIonCharge = 10,
            double precursorIonTolerancePpm = 10,
            double productIonTolerancePpm = 15,
            bool runTargetDecoyAnalysis = true): 
            this(
            specFilePath, 
            dbFilePath, 
            aaSet, 
            minSequenceLength, 
            maxSequenceLength, 
            minPrecursorIonCharge, 
            maxPrecursorIonCharge, 
            minProductIonCharge, 
            maxProductIonCharge, 
            precursorIonTolerancePpm, 
            productIonTolerancePpm, 
            runTargetDecoyAnalysis)
        {
            MaxNumNTermCleavages = maxNumNTermCleavages;
            MaxNumCTermCleavages = maxNumCTermCleavages;
        }

        public string SpecFilePath { get; private set; }
        public string DatabaseFilePath { get; private set; }
        public AminoAcidSet AminoAcidSet { get; private set; }
        public int MinSequenceLength { get; private set; }
        public int MaxSequenceLength { get; private set; }
        public int? MaxNumNTermCleavages { get; private set; }
        public int? MaxNumCTermCleavages { get; private set; }
        public int MinPrecursorIonCharge { get; private set; }
        public int MaxPrecursorIonCharge { get; private set; }
        public int MinProductIonCharge { get; private set; }
        public int MaxProductIonCharge { get; private set; }
        public Tolerance PrecursorIonTolerance { get; private set; }
        public Tolerance ProductIonTolerance { get; private set; }
        public bool RunTargetDecoyAnalysis { get; private set; }

        private LcMsRun _run;
        private ISequenceFilter _ms1Filter;
        private ProductScorerBasedOnDeconvolutedSpectra _ms2ScorerFactory;

        public void RunSearch()
        {
            var sw = new System.Diagnostics.Stopwatch();

            Console.Write("Reading raw file...");
            sw.Start();
            _run = LcMsRun.GetLcMsRun(SpecFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

            sw.Reset();
            Console.Write("Determining precursor masses...");
            sw.Start();
            _ms1Filter = new Ms1IsotopeAndChargeCorrFilter(_run, MinPrecursorIonCharge, MaxPrecursorIonCharge, 10, 3000, 50000, 0.7, 0.7);
            sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

            sw.Reset();
            Console.Write("Deconvoluting MS2 spectra...");
            sw.Start();
            _ms2ScorerFactory = new ProductScorerBasedOnDeconvolutedSpectra(
                _run,
                MinProductIonCharge, MaxProductIonCharge,
                ProductIonTolerance
                );
            _ms2ScorerFactory.DeconvoluteProductSpectra();
            sw.Stop();
            sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

            // Target database
            sw.Reset();
            Console.Write("Reading the target database...");
            sw.Start();
            var targetDb = new FastaDatabase(DatabaseFilePath);
            targetDb.Read();
            var indexedDbTarget = new IndexedDatabase(targetDb);
            var annotationsAndOffsets = MaxNumNTermCleavages != null && MaxNumCTermCleavages != null
                ? indexedDbTarget.IntactSequenceAnnotationsAndOffsets(MinSequenceLength, MaxSequenceLength, (int)MaxNumCTermCleavages) 
                : indexedDbTarget.AnnotationsAndOffsetsNoEnzyme(MinSequenceLength, MaxSequenceLength);
            sw.Stop();
            sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

            sw.Reset();
            Console.WriteLine("Searching the target database");
            sw.Start();
            var targetMatches = RunSearch(annotationsAndOffsets);
            var targetOutputFilePath = Path.ChangeExtension(SpecFilePath, TargetFileExtension);
            WriteResultsToFile(targetMatches, targetOutputFilePath, targetDb);
            sw.Stop();
            sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Target database search elapsed Time: {0:f4} sec", sec);

            // Decoy database
            if (RunTargetDecoyAnalysis)
            {
                sw.Reset();
                Console.Write("Reading the decoy database...");
                sw.Start();
                var decoyDb = targetDb.Decoy(null, true);
                decoyDb.Read();
                var indexedDbDecoy = new IndexedDatabase(decoyDb);
                var annotationsAndOffsetsDecoy = MaxNumNTermCleavages != null
                    ? indexedDbDecoy.IntactSequenceAnnotationsAndOffsets(MinSequenceLength, MaxSequenceLength, (int)MaxNumNTermCleavages)
                    : indexedDbDecoy.AnnotationsAndOffsetsNoEnzyme(MinSequenceLength, MaxSequenceLength);
                sw.Stop();
                sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
                Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

                sw.Reset();
                Console.WriteLine("Searching the decoy database");
                sw.Start();
                var decoyMatches = RunSearch(annotationsAndOffsetsDecoy);
                var decoyOutputFilePath = Path.ChangeExtension(SpecFilePath, DecoyFileExtension);
                WriteResultsToFile(decoyMatches, decoyOutputFilePath, decoyDb);
                sw.Stop();
                sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
                Console.WriteLine(@"Decoy database search elapsed Time: {0:f4} sec", sec);
            }
            Console.WriteLine("Done.");
        }

        private void WriteResultsToFile(SortedSet<ProteinSpectrumMatch>[] matches, string outputFilePath, FastaDatabase database)
        {
            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine("ScanNum\tSequence\tModifications\tComposition\tProteinName\tProteinDesc\tProteinLength\tStart\tEnd\tCharge\tMostAbundantIsotopeMz\tMass\tScore");
                for (var scanNum = _run.MinLcScan; scanNum <= _run.MaxLcScan; scanNum++)
                {
                    if (matches[scanNum] == null) continue;
                    foreach (var match in matches[scanNum].Reverse())
                    {
                        var sequence = match.Sequence;
                        var offset = match.Offset;
                        var start = database.GetZeroBasedPositionInProtein(offset) + 1 + match.NumNTermCleavages;
                        var end = start + sequence.Length - 1;
                        var proteinName = database.GetProteinName(match.Offset);
                        var protLength = database.GetProteinLength(proteinName);
                        var ion = match.Ion;

                        writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}",
                            scanNum,
                            sequence, // Sequence
                            match.Modifications, // Modifications
                            ion.Composition, // Composition
                            proteinName, // ProteinName
                            database.GetProteinDescription(match.Offset), // ProteinDescription
                            protLength, // ProteinLength
                            start, // Start
                            end, // End
                            ion.Charge, // precursorCharge
                            ion.GetMostAbundantIsotopeMz(), // MostAbundantIsotopeMz
                            ion.Composition.Mass,   // Mass
                            match.Score);   // Score
                    }
                }
            }
        }

        private SortedSet<ProteinSpectrumMatch>[] RunSearch(IEnumerable<AnnotationAndOffset> annotationsAndOffsets)
        {
            var sw = new System.Diagnostics.Stopwatch();
            var numProteins = 0;
            sw.Reset();
            sw.Start();

            var matches = new SortedSet<ProteinSpectrumMatch>[_run.MaxLcScan+1];

            foreach (var annotationAndOffset in annotationsAndOffsets)
            {
                ++numProteins;

                var annotation = annotationAndOffset.Annotation;
                var offset = annotationAndOffset.Offset;

                if (numProteins % 100000 == 0)
                {
                    Console.Write("Processing {0}{1} proteins...", numProteins,
                        numProteins == 1 ? "st" : numProteins == 2 ? "nd" : numProteins == 3 ? "rd" : "th");
                    if (numProteins != 0)
                    {
                        sw.Stop();
                        var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
                        Console.WriteLine("Elapsed Time: {0:f4} sec", sec);
                        sw.Reset();
                        sw.Start();
                    }
                }

                var seqGraph = SequenceGraph.CreateGraph(AminoAcidSet, annotation);
                if (seqGraph == null)
                {
                    Console.WriteLine("Ignoring illegal protein: {0}", annotation);
                    continue;
                }

                for (var numNTermCleavage = 0; numNTermCleavage <= (MaxNumNTermCleavages ?? 0); numNTermCleavage++)
                {
                    var protCompositions = seqGraph.GetSequenceCompositionsWithNTermCleavage(numNTermCleavage);
                    for (var modIndex = 0; modIndex < protCompositions.Length; modIndex++)
                    {
                        seqGraph.SetSink(modIndex, numNTermCleavage);
                        var protCompositionWithH2O = seqGraph.GetSinkSequenceCompositionWithH2O();
                        var sequenceMass = protCompositionWithH2O.Mass;
                        var modCombinations = seqGraph.ModificationParams.GetModificationCombination(modIndex);

                        foreach (var ms2ScanNum in _ms1Filter.GetMatchingMs2ScanNums(sequenceMass))
                        {
                            var spec = _run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                            if (spec == null) continue;
                            var charge =
                                (int)Math.Round(sequenceMass / (spec.IsolationWindow.IsolationWindowTargetMz - Constants.Proton));
                            var scorer = _ms2ScorerFactory.GetMs2Scorer(ms2ScanNum);
                            var score = seqGraph.GetScore(charge, scorer);
                            if (score <= 3) continue;

                            var precursorIon = new Ion(protCompositionWithH2O, charge);
                            var sequence = annotation.Substring(numNTermCleavage + 2,
                                annotation.Length - 4 - numNTermCleavage);
                            var prsm = new ProteinSpectrumMatch(sequence, ms2ScanNum, offset, numNTermCleavage, modCombinations,
                                precursorIon, score);

                            
                            if (matches[ms2ScanNum] == null)
                            {
                                matches[ms2ScanNum] = new SortedSet<ProteinSpectrumMatch> { prsm };
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
    }
}
