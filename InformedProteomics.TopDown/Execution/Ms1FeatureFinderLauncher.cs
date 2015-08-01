using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Execution
{
    public class Ms1FeatureFinderLauncher
    {
        private LcMsFeatureLikelihood _likelihoodScorer;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="parameters"></param>
        public Ms1FeatureFinderLauncher(Ms1FeatureFinderInputParameter parameters = null)
        {
            Parameters = parameters ?? new Ms1FeatureFinderInputParameter();
        }

        public void Run()
        {
            //var scoreDataPath = AppDomain.CurrentDomain.BaseDirectory;
            try
            {
                _likelihoodScorer = new LcMsFeatureLikelihood(Parameters.LikelihoodScoreThreshold);
            }
            catch (FileNotFoundException fe)
            {
                ShowErrorMessage(fe.Message);
                return;
            }

            // Normalize the input path. Only affects paths to a file/folder in a folder-type dataset
            Parameters.InputPath = MassSpecDataReaderFactory.NormalizeDatasetPath(Parameters.InputPath);
            
            var attr = File.GetAttributes(Parameters.InputPath);

            if ((attr & FileAttributes.Directory) == FileAttributes.Directory &&
                !MassSpecDataReaderFactory.IsADirectoryDataset(Parameters.InputPath))
            {
                ProcessDirectory(Parameters.InputPath);
            }
            else
            {
                if (!MsRawFile(Parameters.InputPath) && !MsPbfFile(Parameters.InputPath))
                {
                    ShowErrorMessage(@"File extension not supported, " + Parameters.InputPath);
                }
                else
                {
                    ProcessFile(Parameters.InputPath);
                }
            }
        }

        public string GetHeaderString()
        {
            var header = ArrayUtil.ToString(TsvHeader);
            if (Parameters.ScoreReport) header = header + "\t" + ArrayUtil.ToString(TsvExtraScoreHeader);
            return header;
        }

        private void ProcessDirectory(string targetDirectory)
        {
            var fileEntries = Directory.GetFiles(targetDirectory);
            foreach (var fileName in fileEntries)
            {
                if ((MsRawFile(fileName) && !HasExistingPbfFile(fileName)) || MsPbfFile(fileName)) ProcessFile(fileName);
            }
        }

        private bool MsRawFile(string specFilePath)
        {
            //return (path.EndsWith(".raw") || path.EndsWith(".mzML"));
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
            return path.EndsWith(".pbf");
        }
        
        public const string FileExtension = "ms1ft";
        public readonly Ms1FeatureFinderInputParameter Parameters;

        private void ProcessFile(string rawFile)
        {
            var outDirectory = Parameters.OutputPath ?? Path.GetDirectoryName(Path.GetFullPath(rawFile));
            
            var baseName = Path.GetFileName(MassSpecDataReaderFactory.RemoveExtension(rawFile));
            var outTsvFilePath = Path.Combine(outDirectory, baseName + "." + FileExtension);
            var outCsvFilePath = Path.Combine(outDirectory, baseName + "_" + FileExtension + ".csv");
            
            if (File.Exists(outTsvFilePath))
            {
                Console.WriteLine(@"ProMex output already exists: {0}", outTsvFilePath);
                return;
            }

            if (!File.Exists(rawFile))
            {
                ShowErrorMessage(@"Cannot find input file: "+ rawFile);
                return;
            }

            var stopwatch = Stopwatch.StartNew();
            Console.WriteLine(@"Start loading MS1 data from {0}", rawFile);
            var run = PbfLcMsRun.GetLcMsRun(rawFile);

            var featureFinder = new LcMsPeakMatrix(run, _likelihoodScorer);
            Console.WriteLine(@"Complete loading MS1 data. Elapsed Time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds)/1000.0d);

            if (run.GetMs1ScanVector().Length == 0)
            {
                ShowErrorMessage(@"Data file has no MS1 spectra: " + Path.GetFileName(rawFile));
                return;
            }

            var container = new LcMsFeatureContainer(featureFinder.Ms1Spectra, _likelihoodScorer);
            var comparer = featureFinder.Comparer;
            var minSearchMassBin = comparer.GetBinNumber(Parameters.MinSearchMass);
            var maxSearchMassBin = comparer.GetBinNumber(Parameters.MaxSearchMass);
            double totalMassBin = maxSearchMassBin - minSearchMassBin + 1;

            Console.WriteLine(@"Start MS1 feature extraction.");
            stopwatch.Restart();
            for (var binNum = minSearchMassBin; binNum <= maxSearchMassBin; binNum++)
            {
                var clusters = featureFinder.FindFeatures(binNum);
                container.Add(clusters);

                if (binNum > minSearchMassBin && (binNum - minSearchMassBin)%1000 == 0)
                {
                    var elapsed = (stopwatch.ElapsedMilliseconds)/1000.0d;
                    var processedBins = binNum - minSearchMassBin;
                    var processedPercentage = ((double) processedBins/totalMassBin)*100;
                    Console.WriteLine(@"Processing {0:0.0}% of mass bins ({1:0.0} Da); elapsed time = {2:0.000} sec; # of features = {3}",
                        processedPercentage, featureFinder.Comparer.GetMzEnd(binNum), elapsed,
                        container.NumberOfFeatures);
                }
            }

            Console.WriteLine(@"Complete MS1 feature extraction.");
            Console.WriteLine(@" - Elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds)/1000.0d);
            Console.WriteLine(@" - Number of extracted features = {0}", container.NumberOfFeatures);
            Console.WriteLine(@"Start selecting mutually independent features from feature network graph");
            stopwatch.Restart();
            var connectedFeatures = container.GetAllConnectedFeatures();

            // write result files
            var tsvWriter = new StreamWriter(outTsvFilePath);
            tsvWriter.WriteLine(GetHeaderString());

            StreamWriter csvWriter = null;
            if (Parameters.CsvOutput)
            {
                csvWriter = new StreamWriter(outCsvFilePath);
                csvWriter.WriteLine("scan_num,charge,abundance,mz,fit,monoisotopic_mw,FeatureID");
            }

            var featureId = 0;
            foreach (var feature in container.GetFilteredFeatures(connectedFeatures))
            {
                featureId++;
                tsvWriter.WriteLine("{0}\t{1}", featureId, GetString(feature));

                var mostAbuIdx = feature.TheoreticalEnvelope.IndexOrderByRanking[0];

                if (csvWriter != null)
                {
                    foreach (var envelope in feature.EnumerateEnvelopes())
                    {
                        //var mostAbuIsotopeInternalIndex = cluster.IsotopeList.SortedIndexByIntensity[0];
                        var mostAbuPeak = envelope.Peaks[mostAbuIdx];
                        if (mostAbuPeak == null || !mostAbuPeak.Active) continue;

                        var fitscore = 1.0 - feature.BestCorrelationScore;
                        csvWriter.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6}", envelope.ScanNum, envelope.Charge, envelope.Abundance, mostAbuPeak.Mz, fitscore, envelope.MonoMass, featureId));
                    }
                }
            }
            tsvWriter.Close();

            Console.WriteLine(@"Complete feature filteration");
            Console.WriteLine(@" - Elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            Console.WriteLine(@" - Number of filtered features = {0}", featureId);
            Console.WriteLine(@" - ProMex output: {0}", outTsvFilePath);

            if (csvWriter != null)
            {
                csvWriter.Close();
                Console.WriteLine(@" - ProMex output in ICR2LS format: {0}", outCsvFilePath);
            }
        }

        private string GetString(LcMsPeakCluster feature)
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

            sb.AppendFormat("\t{0:0.0}", feature.MinElutionTime);
            sb.AppendFormat("\t{0:0.0}", feature.MaxElutionTime);
            sb.AppendFormat("\t{0:0.0}", feature.ElutionLength);

            sb.Append("\t");
            var intensity = feature.RepresentativeSummedEnvelop;
            var maxIntensity = intensity.Max();
            for (var i = 0; i < intensity.Length; i++)
            {
                if (i != 0) sb.Append(";");
                sb.AppendFormat("{0},{1:0.000}", feature.TheoreticalEnvelope.Isotopes[i].Index, intensity[i] / maxIntensity);
            }

            sb.Append(string.Format("\t{0:0.0}", feature.Score));
            if (Parameters.ScoreReport)
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
            Console.WriteLine(@"----------------------------------------------------------");
            Console.WriteLine(@"Error: " + errorMessage);
            Console.WriteLine(@"----------------------------------------------------------");

        }
        private static readonly string[] TsvHeader = new string[]
        {
            "FeatureID", "MinScan", "MaxScan", "MinCharge", "MaxCharge", 
            "MonoMass", "RepScan", "RepCharge", "RepMz", "Abundance",
            "MinElutionTime", "MaxElutionTime", "ElutionLength", 
            "Envelope", "LikelihoodRatio"
        };

        private static readonly string[] TsvExtraScoreHeader = new string[]
        {
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
