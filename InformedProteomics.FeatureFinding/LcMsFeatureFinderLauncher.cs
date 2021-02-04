using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.FeatureFindingResults;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.FeatureFinding.Clustering;
using InformedProteomics.FeatureFinding.Data;
using InformedProteomics.FeatureFinding.FeatureDetection;
using InformedProteomics.FeatureFinding.Graphics;
using InformedProteomics.FeatureFinding.Scoring;
using InformedProteomics.FeatureFinding.Util;
using PRISM;

namespace InformedProteomics.FeatureFinding
{
    public class LcMsFeatureFinderLauncher
    {
        private readonly LcMsFeatureLikelihood _likelihoodScorer;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="parameters"></param>
        public LcMsFeatureFinderLauncher(LcMsFeatureFinderInputParameters parameters = null)
        {
            Parameters = parameters ?? new LcMsFeatureFinderInputParameters();
            try
            {
                _likelihoodScorer = new LcMsFeatureLikelihood(Parameters.LikelihoodScoreThreshold);
            }
            catch (FileNotFoundException fe)
            {
                ShowErrorMessage(fe.Message);
            }
        }

        /// <summary>
        /// Process the given file or folder
        /// </summary>
        /// <returns>0 if success, otherwise an error code</returns>
        /// <remarks>Folder could either be a supported mass spec data folder, or a normal directory with several supported data files</remarks>
        public int Run()
        {
            // Normalize the input path. Only affects paths to a file/folder in a folder-type dataset
            Parameters.InputPath = MassSpecDataReaderFactory.NormalizeDatasetPath(Parameters.InputPath);

            var attr = File.GetAttributes(Parameters.InputPath);
            int errorCode;

            if ((attr & FileAttributes.Directory) == FileAttributes.Directory &&
                !MassSpecDataReaderFactory.IsADirectoryDataset(Parameters.InputPath))
            {
                errorCode = ProcessDirectory(Parameters.InputPath);
            }
            else
            {
                if (!MsRawFile(Parameters.InputPath) && !MsPbfFile(Parameters.InputPath))
                {
                    ShowErrorMessage(@"File extension not supported, " + Parameters.InputPath);
                    return -1;
                }

                errorCode = ProcessFile(Parameters.InputPath);
            }

            return errorCode;
        }

        public static string GetHeaderString(bool scoreReport = false)
        {
            var header = string.Join("\t", TsvHeader);
            if (scoreReport)
            {
                header = header + "\t" + string.Join("\t", TsvExtraScoreHeader);
            }

            return header;
        }

        /// <summary>
        /// Process each file in the target folder
        /// </summary>
        /// <param name="targetDirectory"></param>
        /// <returns>0 if success, number of errors if a problem</returns>
        private int ProcessDirectory(string targetDirectory)
        {
            var fileEntries = Directory.GetFiles(targetDirectory);
            var errorCount = 0;

            foreach (var fileName in fileEntries)
            {
                if ((MsRawFile(fileName) && !HasExistingPbfFile(fileName)) || MsPbfFile(fileName))
                {
                    var returnCode = ProcessFile(fileName);
                    if (returnCode != 0)
                    {
                        ConsoleMsgUtils.ShowWarning(string.Format("  error code {0} for {1}", returnCode, Path.GetFileName(fileName)));
                        errorCount++;
                    }
                }
            }

            return errorCount;
        }

        private bool MsRawFile(string specFilePath)
        {
            var types = MassSpecDataReaderFactory.MassSpecDataTypeFilterList;
            types.Remove(".pbf");

            // Only supposed to affect execution when running on a directory; however having this test here will affect single file execution
            // i.e., if run with a raw file as input when a .pbf file exists for the dataset, this will return false, and kill the run erroneously.
            //var pbfFilePath = MassSpecDataReaderFactory.ChangeExtension(specFilePath, "pbf");
            //if (File.Exists(pbfFilePath)) return false;

            return types.Any(ext => specFilePath.ToLower().EndsWith(ext));
            /*
            Console.WriteLine("--------------");
            Console.WriteLine(specFilePath);
            Console.WriteLine(types.Select(ext => specFilePath.ToLower().EndsWith(ext)).Any(x => x == true));
            foreach (var x in types.Select(ext => specFilePath.ToLower().EndsWith(ext)))
            {
                Console.WriteLine(x);
            }

            Console.WriteLine("--------------");
            */
        }

        private bool HasExistingPbfFile(string path)
        {
            return File.Exists(MassSpecDataReaderFactory.ChangeExtension(path, ".pbf"));
        }

        private bool MsPbfFile(string path)
        {
            return path.EndsWith(".pbf", StringComparison.OrdinalIgnoreCase);
        }

        public const string FileExtension = "ms1ft";
        public readonly LcMsFeatureFinderInputParameters Parameters;

        /// <summary>
        /// Find features in the data file
        /// </summary>
        /// <param name="rawFile">Data file (either a pbf file or a file type from which a pbf file can be auto-created)</param>
        /// <returns>0 if success; negative number on error</returns>
        private int ProcessFile(string rawFile)
        {
            var outDirectory = GetOutputDirectory(rawFile);
            if (string.IsNullOrEmpty(outDirectory))
            {
                return -1;
            }

            var baseName = Path.GetFileName(MassSpecDataReaderFactory.RemoveExtension(rawFile));
            var ms1FeaturesFilePath = Path.Combine(outDirectory, baseName + "." + FileExtension);
            var outCsvFilePath = Path.Combine(outDirectory, baseName + "_" + FileExtension + ".csv");
            var pngFilePath = Path.Combine(outDirectory, baseName + "_" + FileExtension + ".png");

            if (File.Exists(ms1FeaturesFilePath))
            {
                ShowErrorMessage("ProMex output already exists: " + ms1FeaturesFilePath);
                return -2;
            }

            if (!File.Exists(rawFile))
            {
                ShowErrorMessage("Cannot find input file: " + rawFile);
                return -3;
            }

            var stopwatch = Stopwatch.StartNew();
            Console.WriteLine("Start loading MS1 data from {0}", rawFile);
            var run = PbfLcMsRun.GetLcMsRun(rawFile);

            var featureFinder = new LcMsPeakMatrix(run, _likelihoodScorer, 1, 60, Parameters.MaxThreads);
            Console.WriteLine("Complete loading MS1 data. Elapsed Time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);

            if (run.GetMs1ScanVector().Length == 0)
            {
                ShowErrorMessage(@"Data file has no MS1 spectra: " + Path.GetFileName(rawFile));
                return -4;
            }

            if (featureFinder.Ms1PeakCount == 0)
            {
                ShowErrorMessage(@"Data file has no MS1 peaks: " + Path.GetFileName(rawFile));
                return -5;
            }

            var comparer = featureFinder.Comparer;
            var container = new LcMsFeatureContainer(featureFinder.Ms1Spectra, _likelihoodScorer, new LcMsFeatureMergeComparer(new Tolerance(10)));
            var minSearchMassBin = comparer.GetBinNumber(Parameters.MinSearchMass);
            var maxSearchMassBin = comparer.GetBinNumber(Parameters.MaxSearchMass);
            double totalMassBin = maxSearchMassBin - minSearchMassBin + 1;

            Console.WriteLine("Start MS1 feature extraction.");
            stopwatch.Restart();
            for (var binNum = minSearchMassBin; binNum <= maxSearchMassBin; binNum++)
            {
                var clusters = featureFinder.FindFeatures(binNum);
                container.Add(clusters);

                if (binNum > minSearchMassBin && (binNum - minSearchMassBin) % 1000 == 0)
                {
                    var elapsed = (stopwatch.ElapsedMilliseconds) / 1000.0d;
                    var processedBins = binNum - minSearchMassBin;
                    var processedPercentage = processedBins / totalMassBin * 100;
                    Console.WriteLine("Processing {0:0.0}% of mass bins ({1:0.0} Da); elapsed time = {2:0.000} sec; # of features = {3}",
                        processedPercentage, featureFinder.Comparer.GetMzEnd(binNum), elapsed,
                        container.NumberOfFeatures);
                }
            }

            Console.WriteLine("Complete MS1 feature extraction.");
            Console.WriteLine(" - Elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            Console.WriteLine(" - Number of extracted features = {0}", container.NumberOfFeatures);
            Console.WriteLine("Start selecting mutually independent features from feature network graph");
            stopwatch.Restart();

            var featureId = FilterAndOutputFeatures(container, featureFinder, outCsvFilePath, ms1FeaturesFilePath);

            Console.WriteLine("Complete feature filtration");
            Console.WriteLine(" - Elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            Console.WriteLine(" - Number of filtered features = {0}", featureId);
            Console.WriteLine(" - ProMex output: {0}", ms1FeaturesFilePath);

            if (Parameters.CsvOutput)
            {
                Console.WriteLine(" - ProMex output in ICR2LS format: {0}", outCsvFilePath);
            }

            if (Parameters.FeatureMapImage)
            {
                CreateFeatureMapImage(run, ms1FeaturesFilePath, pngFilePath);
            }

            return 0;
        }

        [Obsolete("Unused")]
        private int FilterAndOutputFeaturesOld(LcMsFeatureContainer container, LcMsPeakMatrix featureFinder, string outCsvFilePath, string ms1FeaturesFilePath)
        {
            var featureId = 0;

            Stream csvStream = new MemoryStream();
            if (Parameters.CsvOutput)
            {
                csvStream = new FileStream(outCsvFilePath, FileMode.Create, FileAccess.Write, FileShare.ReadWrite);
            }
            // write result files
            using (var tsvWriter = new StreamWriter(ms1FeaturesFilePath))
            using (var csvWriter = new StreamWriter(csvStream))
            {
                tsvWriter.WriteLine(GetHeaderString(Parameters.ScoreReport));

                if (Parameters.CsvOutput)
                {
                    csvWriter.WriteLine("scan_num,charge,abundance,mz,fit,monoisotopic_mw,FeatureID");
                }

                var filteredFeatures = container.GetFilteredFeatures(featureFinder);
                foreach (var feature in filteredFeatures)
                {
                    featureId++;
                    tsvWriter.WriteLine("{0}\t{1}", featureId, GetString(feature, Parameters.ScoreReport));

                    var mostAbuIdx = feature.TheoreticalEnvelope.IndexOrderByRanking[0];

                    if (Parameters.CsvOutput)
                    {
                        foreach (var envelope in feature.EnumerateEnvelopes())
                        {
                            //var mostAbuIsotopeInternalIndex = cluster.IsotopeList.SortedIndexByIntensity[0];
                            var mostAbuPeak = envelope.Peaks[mostAbuIdx];
                            if (mostAbuPeak == null || !mostAbuPeak.Active)
                            {
                                continue;
                            }

                            var fitScore = 1.0 - feature.BestCorrelationScore;
                            csvWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6}", envelope.ScanNum, envelope.Charge, envelope.Abundance,
                                mostAbuPeak.Mz, fitScore, envelope.MonoMass, featureId);
                        }
                    }
                }
            }

            return featureId;
        }

        private int FilterAndOutputFeatures(LcMsFeatureContainer container, LcMsPeakMatrix featureFinder, string outCsvFilePath, string ms1FeaturesFilePath)
        {
            var featureCounter = new int[1];
            Ms1FtEntry.WriteToFile(ms1FeaturesFilePath, FilterFeaturesWithOutput(container, featureFinder, outCsvFilePath, featureCounter), Parameters.ScoreReport);

            return featureCounter[0];
        }

        /// <summary>
        /// Create <see cref="Ms1FtEntry"/> objects for the features, and output to csv (if desired).
        /// </summary>
        /// <param name="container"></param>
        /// <param name="featureFinder"></param>
        /// <param name="outCsvFilePath"></param>
        /// <param name="featureCounter"></param>
        /// <returns></returns>
        private IEnumerable<Ms1FtEntry> FilterFeaturesWithOutput(
            LcMsFeatureContainer container,
            LcMsPeakMatrix featureFinder,
            string outCsvFilePath,
            IList<int> featureCounter)
        {
            // Using an array, since we can't use ref or out parameters
            featureCounter[0] = 0;

            Stream csvStream = new MemoryStream();
            if (Parameters.CsvOutput)
            {
                csvStream = new FileStream(outCsvFilePath, FileMode.Create, FileAccess.Write, FileShare.ReadWrite);
            }
            using (var csvWriter = new StreamWriter(csvStream))
            {
                if (Parameters.CsvOutput)
                {
                    csvWriter.WriteLine("scan_num,charge,abundance,mz,fit,monoisotopic_mw,FeatureID");
                }

                var filteredFeatures = container.GetFilteredFeatures(featureFinder);
                foreach (var feature in filteredFeatures)
                {
                    featureCounter[0]++;

                    if (Parameters.CsvOutput)
                    {
                        var mostAbuIdx = feature.TheoreticalEnvelope.IndexOrderByRanking[0];

                        foreach (var envelope in feature.EnumerateEnvelopes())
                        {
                            //var mostAbuIsotopeInternalIndex = cluster.IsotopeList.SortedIndexByIntensity[0];
                            var mostAbuPeak = envelope.Peaks[mostAbuIdx];
                            if (mostAbuPeak == null || !mostAbuPeak.Active)
                            {
                                continue;
                            }

                            var fitScore = 1.0 - feature.BestCorrelationScore;
                            csvWriter.WriteLine("{0},{1},{2},{3},{4},{5},{6}", envelope.ScanNum, envelope.Charge, envelope.Abundance,
                                mostAbuPeak.Mz, fitScore, envelope.MonoMass, featureCounter);
                        }
                    }

                    yield return feature.ToMs1FtEntry(featureCounter[0]);
                }
            }
        }

        /// <summary>
        /// Create a PNG image of previously found MS1 features
        /// </summary>
        /// <param name="pbfFilePath">.pbf file path</param>
        /// <param name="ms1FeaturesFilePath">.ms1ft file path</param>
        /// <returns>0 if success, otherwise an error code</returns>
        /// <remarks>
        /// If ms1FeaturesFilePath is an empty string, it is auto-determined based on the .pbf file name
        /// Mass range is determined using Parameters.MinSearchMass and Parameters.MaxSearchMass
        /// </remarks>
        public int CreateFeatureMapImage(string pbfFilePath, string ms1FeaturesFilePath)
        {
            if (!File.Exists(pbfFilePath))
            {
                ShowErrorMessage("Data file not found: " + pbfFilePath);
                return -1;
            }

            var outDirectory = GetOutputDirectory(pbfFilePath);
            if (string.IsNullOrEmpty(outDirectory))
            {
                return -2;
            }

            var baseName = Path.GetFileName(MassSpecDataReaderFactory.RemoveExtension(pbfFilePath));

            if (string.IsNullOrEmpty(ms1FeaturesFilePath) || string.Equals(ms1FeaturesFilePath, "."))
            {
                ms1FeaturesFilePath = Path.Combine(outDirectory, baseName + "." + FileExtension);
            }

            var pngFilePath = Path.Combine(outDirectory, baseName + "_" + FileExtension + ".png");

            if (!File.Exists(ms1FeaturesFilePath))
            {
                ShowErrorMessage("MS1 features file not found: " + ms1FeaturesFilePath);
                return -3;
            }

            Console.WriteLine("Start loading MS1 data from {0}", pbfFilePath);
            var run = PbfLcMsRun.GetLcMsRun(pbfFilePath);

            CreateFeatureMapImage(run, ms1FeaturesFilePath, pngFilePath);

            return 0;
        }

        private void CreateFeatureMapImage(LcMsRun run, string featuresFilePath, string imgFilePath)
        {
            var map = new LcMsFeatureMap(run, featuresFilePath, Math.Max(0, Parameters.MinSearchMass - 500), Parameters.MaxSearchMass);
            map.SaveImage(imgFilePath);
            Console.WriteLine(" - Feature map image output: {0}", imgFilePath);
        }

        private string GetOutputDirectory(string dataFilePath)
        {
            var outDirectory = Parameters.OutputPath ?? Path.GetDirectoryName(Path.GetFullPath(dataFilePath));

            if (string.IsNullOrEmpty(outDirectory))
            {
                Console.WriteLine("Unable to determine the folder path for the data file: {0}", dataFilePath);
                return string.Empty;
            }

            return outDirectory;
        }

        public static string GetString(LcMsFeature feature)
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

            sb.AppendFormat("\t{0:0}", 0);
            sb.AppendFormat("\t{0:0.00}", 0);

            sb.AppendFormat("\t{0:0.0}", feature.MinElutionTime);
            sb.AppendFormat("\t{0:0.0}", feature.MaxElutionTime);
            sb.AppendFormat("\t{0:0.0}", feature.ElutionLength);

            sb.Append("\t");
            sb.Append("");
            sb.Append(string.Format("\t{0:0.0}", feature.Score));
            return sb.ToString();
        }

        public static string GetString(LcMsPeakCluster feature, bool scoreReport = false)
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

            sb.AppendFormat("\t{0:0}", feature.ApexScanNum);
            sb.AppendFormat("\t{0:0.00}", feature.ApexIntensity);

            sb.AppendFormat("\t{0:0.000}", feature.MinElutionTime);
            sb.AppendFormat("\t{0:0.000}", feature.MaxElutionTime);
            sb.AppendFormat("\t{0:0.000}", feature.ElutionLength);

            sb.Append("\t");
            var intensity = feature.RepresentativeSummedEnvelop;
            var maxIntensity = intensity.Max();
            for (var i = 0; i < intensity.Length; i++)
            {
                if (i != 0)
                {
                    sb.Append(";");
                }

                sb.AppendFormat("{0},{1:0.000}", feature.TheoreticalEnvelope.Isotopes[i].Index, intensity[i] / maxIntensity);
            }

            // LikelihoodRatio
            sb.Append(string.Format("\t{0:0.0000}", feature.Score));

            if (scoreReport)
            {
                sb.AppendFormat("\t{0}", feature.BestCharge[LcMsPeakCluster.EvenCharge]);
                sb.AppendFormat("\t{0}", feature.BestCharge[LcMsPeakCluster.OddCharge]);

                sb.AppendFormat("\t{0:0.000}", feature.BestCorrelationScoreAcrossCharge[LcMsPeakCluster.EvenCharge]);
                sb.AppendFormat("\t{0:0.000}", feature.BestCorrelationScoreAcrossCharge[LcMsPeakCluster.OddCharge]);

                sb.AppendFormat("\t{0:0.000}", feature.BestIntensityScoreAcrossCharge[LcMsPeakCluster.EvenCharge]);
                sb.AppendFormat("\t{0:0.000}", feature.BestIntensityScoreAcrossCharge[LcMsPeakCluster.OddCharge]);

                sb.AppendFormat("\t{0:0.000}", feature.EnvelopeCorrelationScoreAcrossCharge[LcMsPeakCluster.EvenCharge]);
                sb.AppendFormat("\t{0:0.000}", feature.EnvelopeCorrelationScoreAcrossCharge[LcMsPeakCluster.OddCharge]);

                sb.AppendFormat("\t{0:0.000}", feature.EnvelopeIntensityScoreAcrossCharge[LcMsPeakCluster.EvenCharge]);
                sb.AppendFormat("\t{0:0.000}", feature.EnvelopeIntensityScoreAcrossCharge[LcMsPeakCluster.OddCharge]);

                sb.AppendFormat("\t{0:0.000}", feature.XicCorrelationBetweenBestCharges[LcMsPeakCluster.EvenCharge]);
                sb.AppendFormat("\t{0:0.000}", feature.XicCorrelationBetweenBestCharges[LcMsPeakCluster.OddCharge]);

                sb.AppendFormat("\t{0:0.000}", feature.AbundanceDistributionAcrossCharge[LcMsPeakCluster.EvenCharge]);
                sb.AppendFormat("\t{0:0.000}", feature.AbundanceDistributionAcrossCharge[LcMsPeakCluster.OddCharge]);
            }

            return sb.ToString();
        }

        private void ShowErrorMessage(string errorMessage)
        {
            ConsoleMsgUtils.ShowError("Error: " + errorMessage);
        }

        public static readonly string[] TsvHeader = {
            "FeatureID", "MinScan", "MaxScan", "MinCharge", "MaxCharge",
            "MonoMass", "RepScan", "RepCharge", "RepMz", "Abundance",
            "ApexScanNum", "ApexIntensity",
            "MinElutionTime", "MaxElutionTime", "ElutionLength", "Envelope", "LikelihoodRatio"
        };

        public static readonly string[] TsvExtraScoreHeader = {
            "BestEvenCharge", "BestOddCharge",
            "CorrEvenCharge", "CorrOddCharge",
            "IntensityEvenCharge", "IntensityOddCharge",
            "SummedCorrEvenCharge", "SummedCorrOddCharge",
            "SummedIntensityEvenCharge", "SummedIntensityOddCharge",
            "XicCorrBetCharges1", "XicCorrBetCharges2",
            "AbundanceRatioEvenCharge", "AbundanceRatioOddCharge",
        };
    }
}
