using System;
using System.Collections;
using System.Diagnostics;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Execution
{
    public class QuickProteinId
    {
        public LcMsRun Run { get; private set; }

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
      
        public void QuickId()
        {
            const string rawFilePath = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            const string fastaFilePath = @"H:\Research\QCShew_TopDown\Production\ID_002216_235ACCEA.fasta";
            const string modFilePath = @"H:\Research\QCShew_TopDown\Production\Mods.txt";
            const int numBits = 29; // max error: 4ppm
            const int minCharge = 1;
            const int maxCharge = 20;
            var tolerance = new Tolerance(10);
            const double corrThreshold = 0.7;

            var comparer = new MzComparerWithBinning(numBits);
            const double minFragmentMass = 200.0;
            const double maxFragmentMass = 50000.0;
            var minFragMassBin = comparer.GetBinNumber(minFragmentMass);
            var maxFragMassBin = comparer.GetBinNumber(maxFragmentMass);

            var aminoAcidSet = new AminoAcidSet(modFilePath);

            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);
            var ms2ScanNumArr = run.GetScanNumbers(2).ToArray();

            var sw = new Stopwatch();

            sw.Start();
            Console.Write("Building Spectrum Arrays...");
            var massVectors = new BitArray[maxFragMassBin - minFragMassBin + 1];
            for (var i = minFragMassBin; i <= maxFragMassBin; i++)
            {
                massVectors[i - minFragMassBin] = new BitArray(run.MaxLcScan + 1);
            }

            foreach (var ms2ScanNum in ms2ScanNumArr)
            {
                var productSpec = run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
                if (productSpec == null) continue;

                var deconvolutedPeaks = Deconvoluter.GetDeconvolutedPeaks(productSpec.Peaks, minCharge, maxCharge, 2, 1.1, tolerance, corrThreshold);

                if (deconvolutedPeaks == null) continue;

                foreach (var p in deconvolutedPeaks)
                {
                    var mass = p.Mass;
                    var deltaMass = tolerance.GetToleranceAsDa(mass, 1);
                    var minMass = mass - deltaMass;
                    var maxMass = mass + deltaMass;

                    var minBinNum = comparer.GetBinNumber(minMass);
                    var maxBinNum = comparer.GetBinNumber(maxMass);
                    for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
                    {
                        if (binNum >= minFragMassBin && binNum <= maxFragMassBin) massVectors[binNum - minFragMassBin][ms2ScanNum] = true;
                    }
                }
            }
            sw.Stop();
            Console.WriteLine(@"{0:f1} sec.", sw.Elapsed.TotalSeconds);

            sw.Reset();
            sw.Start();
            var fastaDb = new FastaDatabase(fastaFilePath);
            fastaDb.Read();
            var indexedDb = new IndexedDatabase(fastaDb);
            var numProteins = 0;
            var intactProteinAnnotationAndOffsets =
                indexedDb.IntactSequenceAnnotationsAndOffsets(0, int.MaxValue);

            var bestProtein = new string[run.MaxLcScan + 1];
            var bestScore = new int[run.MaxLcScan + 1];
            foreach (var annotationAndOffset in intactProteinAnnotationAndOffsets)
            {
                if (++numProteins % 10 == 0)
                {
                    Console.WriteLine(@"Processing, {0} proteins done, {1:f1} sec elapsed",
                        numProteins,
                        sw.Elapsed.TotalSeconds);           
                }
                var annotation = annotationAndOffset.Annotation;
                var offset = annotationAndOffset.Offset;

                var protSequence = annotation.Substring(2, annotation.Length - 4);

                // suffix
                var seqGraph = SequenceGraph.CreateGraph(aminoAcidSet, AminoAcid.ProteinNTerm, protSequence,
                    AminoAcid.ProteinCTerm);
                if (seqGraph == null) continue;

                for (var numNTermCleavage = 0; numNTermCleavage <= 1; numNTermCleavage++)
                {
                    if (numNTermCleavage > 0) seqGraph.CleaveNTerm();
                    var allCompositions = seqGraph.GetAllFragmentNodeCompositions();

                    var scoreArr = new int[run.MaxLcScan + 1];
                    foreach (var fragComp in allCompositions)
                    {
                        var suffixMass = fragComp.Mass + BaseIonType.Y.OffsetComposition.Mass;
                        var binNum = comparer.GetBinNumber(suffixMass);
                        if (binNum < minFragMassBin || binNum > maxFragMassBin) continue;

                        var vector = massVectors[binNum - minFragMassBin];
                        foreach (var ms2ScanNum in ms2ScanNumArr)
                        {
                            if (vector[ms2ScanNum]) ++scoreArr[ms2ScanNum];
                        }
                    }
                    foreach (var ms2ScanNum in ms2ScanNumArr)
                    {
                        if (scoreArr[ms2ScanNum] > bestScore[ms2ScanNum])
                        {
                            bestScore[ms2ScanNum] = scoreArr[ms2ScanNum];
                            var proteinName = fastaDb.GetProteinName(offset);
                            bestProtein[ms2ScanNum] = proteinName + (numNTermCleavage == 1 ? "'" : "");
                        }
                    }
                }
                // prefix
            }

            Console.WriteLine("ScanNum\tBestProtein\tScore");
            foreach (var ms2ScanNum in ms2ScanNumArr)
            {
                Console.WriteLine("{0}\t{1}\t{2}", ms2ScanNum, bestScore[ms2ScanNum], bestProtein[ms2ScanNum] ?? "");
            }            
        }
    }
}
