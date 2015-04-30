using System.Collections.Generic;
//using LibSVMsharp;
//using LibSVMsharp.Extensions;
using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Execution
{
    public class Ms1FeatureFinderLauncher
    {
        public Ms1FeatureFinderLauncher(Ms1FeatureFinderInputParameter parameters = null)
        {
            Parameters = parameters ?? new Ms1FeatureFinderInputParameter();
        }

        public void Run()
        {
            var attr = File.GetAttributes(Parameters.InputPath);

            if ((attr & FileAttributes.Directory) == FileAttributes.Directory)
            {
                ProcessDirectory(Parameters.InputPath);
            }
            else
            {
                if (!MsRawFile(Parameters.InputPath) && !MsPbfFile(Parameters.InputPath))
                {
                    Console.WriteLine("Not supported file extension");
                }
                else
                {
                    ProcessFile(Parameters.InputPath);
                }
            }
        }

        private void ProcessDirectory(string targetDirectory)
        {
            var fileEntries = Directory.GetFiles(targetDirectory);
            foreach (string fileName in fileEntries)
            {
                if (MsRawFile(fileName))
                {
                    var pbfFilePath = Path.ChangeExtension(fileName, "pbf");
                    if (!File.Exists(pbfFilePath)) ProcessFile(fileName);
                }
                else if (MsPbfFile(fileName)) ProcessFile(fileName);
            }
        }

        private bool MsRawFile(string path)
        {
            return (path.EndsWith(".raw") || path.EndsWith(".mzML"));
        }

        private bool MsPbfFile(string path)
        {
            return path.EndsWith(".pbf");
        }

        private void ProcessFile(string rawFile)
        {
            var outTsvFilePath = Path.ChangeExtension(rawFile, FileExtension);
            var j = rawFile.LastIndexOf('.');
            var outPath = rawFile.Substring(0, j);
            var outCsvFilePath = string.Format("{0}_{1}.csv", outPath, FileExtension);


            if (File.Exists(outTsvFilePath))
            {
                Console.WriteLine("ProMex output already exists: {0}", outTsvFilePath);
                return;
            }

            if (!File.Exists(rawFile))
            {
                Console.WriteLine("Cannot find input file: {0}", rawFile);
                return;
            }

            var stopwatch = Stopwatch.StartNew();
            Console.WriteLine("Start loading MS1 data from {0}", rawFile);
            var run = PbfLcMsRun.GetLcMsRun(rawFile, rawFile.EndsWith(".mzML") ? MassSpecDataType.MzMLFile : MassSpecDataType.XCaliburRun);

            var extractor = new Ms1FeatureMatrix(run, Parameters.MinSearchCharge, Parameters.MaxSearchCharge, Parameters.MaxThreads);
            Console.WriteLine("Complete loading MS1 data. Elapsed Time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            //extractor.GetFeatureFile(rawFile, _minSearchMass, _maxSearchMass, _scoreReport, _csvOutput, _tmpOutput);
            //var outTsvFilePath = Path.ChangeExtension(rawFile, FileExtension);
            //if (File.Exists(outTsvFilePath)) return outTsvFilePath;


            var container = new Ms1FeatureContainer(extractor.Spectrums);
            var minSearchMassBin = extractor.Comparer.GetBinNumber(Parameters.MinSearchMass);
            var maxSearchMassBin = extractor.Comparer.GetBinNumber(Parameters.MaxSearchMass);
            double totalMassBin = maxSearchMassBin - minSearchMassBin + 1;

            Console.WriteLine("Start MS1 feature extraction.");
            //if (Parameters.CsvOutput) Console.WriteLine("Csv Output File\t{0}", outCsvFilePath);

            stopwatch.Restart();
            for (var binNum = minSearchMassBin; binNum <= maxSearchMassBin; binNum++)
            {
                var clusters = extractor.GetProbableClusters(binNum);
                container.Add(clusters);

                if (binNum > minSearchMassBin && (binNum - minSearchMassBin) % 3000 == 0)
                {
                    var elapsed = (stopwatch.ElapsedMilliseconds) / 1000.0d;
                    var processedBins = binNum - minSearchMassBin;
                    var processedPercentage = ((double)processedBins / totalMassBin) * 100;
                    Console.WriteLine("Processing {0:0.0}% of mass bins ({1:0.0} Da); elapsed time = {2:0.000} sec; # of features = {3}", processedPercentage, extractor.Comparer.GetMzEnd(binNum), elapsed, container.NumberOfFeatures);
                }
            }

            Console.WriteLine("Complete MS1 feature extraction.");
            Console.WriteLine(" - Elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            Console.WriteLine(" - Number of extracted features = {0}", container.NumberOfFeatures);

            // write result files
            StreamWriter csvWriter = null;
            if (Parameters.CsvOutput)
            {
                csvWriter = new StreamWriter(outCsvFilePath);
                csvWriter.WriteLine("scan_num,charge,abundance,mz,fit,monoisotopic_mw,FeatureID");
            }

            Console.WriteLine("Start selecting mutually independent features from feature network graph");
            var connectedFeatures = container.GetAllConnectedFeatures();
            if (Parameters.TmpOutput)
            {
                var tmpTsvFilePath = string.Format("{0}.tmp", outTsvFilePath);
                var tmpTsvWriter = new StreamWriter(tmpTsvFilePath);
                tmpTsvWriter.WriteLine(GetHeaderString());
                var clusterId = 0;
                foreach (var featureSet in connectedFeatures)
                {
                    clusterId++;
                    foreach (var feature in featureSet)
                    {
                        tmpTsvWriter.WriteLine("{0}\t{1}", clusterId, GetString(feature));
                    }
                }
                tmpTsvWriter.Close();
            }

            var filteredFeatures = container.GetFilteredFeatures(connectedFeatures);
            stopwatch.Stop();

            var featureId = 0;
            var ms1ScanNums = run.GetMs1ScanVector();
            var tsvWriter = new StreamWriter(outTsvFilePath);
            tsvWriter.WriteLine(GetHeaderString());
            foreach (var cluster in filteredFeatures)
            {
                featureId++;
                tsvWriter.WriteLine("{0}\t{1}", featureId, GetString(cluster));

                if (csvWriter != null)
                {
                    foreach (var envelope in cluster.Envelopes)
                    {
                        //var mostAbuIsotopeInternalIndex = cluster.IsotopeList.SortedIndexByIntensity[0];
                        var refInternalIndex = envelope.RefIsotopeInternalIndex;
                        var mostAbuPeak = envelope.Peaks[refInternalIndex];
                        var charge = envelope.Row + extractor.ChargeRange.Min;
                        var mass = Ion.GetMonoIsotopicMass(mostAbuPeak.Mz, charge, cluster.TheoreticalEnvelope[refInternalIndex].Index);
                        var abundance = envelope.Abundance;
                        csvWriter.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6}", ms1ScanNums[envelope.Col], charge, abundance, mostAbuPeak.Mz, 1.0 - cluster.Probability, mass, featureId));
                    }
                }
            }
            tsvWriter.Close();

            Console.WriteLine("Complete feature filteration");
            Console.WriteLine(" - Elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            Console.WriteLine(" - Number of filtered features = {0}", filteredFeatures.Count);
            Console.WriteLine(" - ProMex output: {0}", outTsvFilePath);

            if (Parameters.CsvOutput)
            {
                csvWriter.Close();
                Console.WriteLine(" - ProMex output in ICR2LS format: {0}", outCsvFilePath);
            }

            //return outTsvFilePath;
        }

        public const string FileExtension = "ms1ft";
        public readonly Ms1FeatureFinderInputParameter Parameters;

        public static readonly string[] TsvHeaderWithScore = new string[]
        {
            "FeatureID", "MinScan", "MaxScan", "MinCharge", "MaxCharge", "MonoMass", "RepScan", "RepCharge", "RepMz", "Abundance",
            "ObservedEnvelopes", "ScanLength",
            "BestCorr", "SummedCorr", 
            "RanksumScore", "PoissonScore",
            "BcDist", "SummedBcDist", 
            "SummedBcDistOverCharges", "SummedBcDistOverTimes", 
            "XicDist", "MinXicDist",
            "MzDiffPpm", "TotalMzDiffPpm",
            "bcDistEvenCharge", "bcDistanceOddCharge", "AbundanceChange",
            "Envelope", 
            "Probability", "GoodEnough"
        };

        private static readonly string[] TsvHeader = new string[]
        {
            "FeatureID", "MinScan", "MaxScan", "MinCharge", "MaxCharge", "MonoMass", "RepScan", "RepCharge", "RepMz", "Abundance",
            "BestCorr", "SummedCorr", 
            "XicDist", 
            "MzDiffPpm", "Envelope",
            "Probability", "GoodEnough"
        };

        public string GetHeaderString()
        {
            return (Parameters.ScoreReport) ? ArrayUtil.ToString(TsvHeaderWithScore) : ArrayUtil.ToString(TsvHeader);
        }

        public string GetString(Ms1FeatureCluster feature)
        {
            // should be called after calling UpdateScore & UpdateAbundance
            var sb = new StringBuilder(string.Format("{0}\t{1}\t{2}\t{3}\t{4:0.0000}\t{5}\t{6}\t{7:0.0000}\t{8:0.00}",
                                        feature.MinScanNum, feature.MaxScanNum,
                                        feature.MinCharge, feature.MaxCharge,
                                        feature.RepresentativeMass,
                                        feature.RepresentativeScanNum,
                                        feature.RepresentativeCharge,
                                        feature.RepresentativeMz,
                                        feature.Abundance));

            if (Parameters.ScoreReport)
            {
                sb.AppendFormat("\t{0}", feature.GoodEnvelopeCount);
                sb.AppendFormat("\t{0}", feature.ScanLength);
                //foreach (var score in _scores) sb.AppendFormat("\t{0}", score);
                for (var i = 0; i < Ms1FeatureScore.Count; i++) sb.AppendFormat("\t{0}", feature.GetScore((byte)i));
            }
            else
            {
                sb.AppendFormat("\t{0:0.00000}", feature.GetScore(Ms1FeatureScore.EnvelopeCorrelation));
                sb.AppendFormat("\t{0:0.00000}", feature.GetScore(Ms1FeatureScore.EnvelopeCorrelationSummed));
                sb.AppendFormat("\t{0:0.00000}", feature.GetScore(Ms1FeatureScore.TotalMzError));
                sb.AppendFormat("\t{0:0.00000}", feature.GetScore(Ms1FeatureScore.XicCorrMean));
            }

            sb.Append("\t");

            var intensity = feature.SummedEnvelope;
            var maxIntensity = intensity.Max();
            for (var i = 0; i < intensity.Length; i++)
            {
                if (i != 0) sb.Append(";");
                sb.AppendFormat("{0},{1:0.000}", feature.TheoreticalEnvelope[i].Index, intensity[i] / maxIntensity);
            }

            sb.AppendFormat("\t{0:0.0000}\t{1}", feature.Probability, feature.GoodEnough ? 1 : 0);
            return sb.ToString();
        }

    }
}
