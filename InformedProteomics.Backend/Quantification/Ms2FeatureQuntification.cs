using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Quantification;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Ms2FeatureQuntification: Ms1FeatureMatrix
    {
        public Ms2FeatureQuntification(LcMsRun run, int minScanCharge = 2, int maxScanCharge = 60, int maxThreadCount = 0)
            : base(run, minScanCharge, maxScanCharge, maxThreadCount)
        {

        }

        public double GetMs1EvidenceScore(Ms2Feature ms2Feature)
        {
            SetQueryMass(ms2Feature.Mass);

            var totalElutionTime = Run.GetElutionTime(Run.MaxLcScan);
            var elutionTime = Run.GetElutionTime(ms2Feature.ScanNum);
            var charge = ms2Feature.Charge;
            var ms1ScanNums = Run.GetMs1ScanVector();

            var minElutionTime = elutionTime - totalElutionTime * 0.003;
            var maxElutionTime = elutionTime + totalElutionTime * 0.003;

            var ms1ScanIndex = Array.BinarySearch(ms1ScanNums, ms2Feature.ScanNum);
            if (ms1ScanIndex < 0) ms1ScanIndex = ~ms1ScanIndex; // next Ms1 scan num

            var bcDistances = new List<double>();
            var correlations = new List<double>();
            var enelopes    = new List<ObservedEnvelope>();
            var maxSearchScans = (int) Math.Max(ms1ScanNums.Length - ms1ScanIndex + 1, ms1ScanIndex);
            var reachRight = false;
            var reachLeft = false;
            
            for (var i = 0; i <= maxSearchScans; i++)
            {
                for (var j = 0; j < 2; j++)
                {
                    if (i == 0 && j > 0) continue;
                    
                    if (reachRight && reachLeft) break;

                    var col = (j < 1) ? ms1ScanIndex + i : ms1ScanIndex - i;
                    if (!reachRight && Run.GetElutionTime(ms1ScanNums[col]) > maxElutionTime)
                    {
                        reachRight = true;
                        continue;
                    }
                    if (!reachLeft && Run.GetElutionTime(ms1ScanNums[col]) < minElutionTime)
                    {
                        reachLeft = true;
                        continue;
                    }

                    var ms1Spec = Spectrums[col];
                    var observedPeaks = ms1Spec.GetAllIsotopePeaks(QueryMass, charge, TheoreticalEnvelope, MzTolerance);
                    var obsEnv = new ObservedEnvelope(charge - ChargeRange.Min, col, observedPeaks, TheoreticalEnvelope);
                    if (obsEnv.NumberOfPeaks < 3) continue;
                    var bcDist = obsEnv.GetBhattacharyyaDistance(TheoreticalEnvelope.EnvelopePdf);
                    var corrCoeff = obsEnv.GetPearsonCorrelation(TheoreticalEnvelope.Envelope);

                    if (bcDist < 0.1 || corrCoeff > 0.8)
                    {
                        return Math.Max(corrCoeff, 1 - bcDist);
                    }

                    bcDistances.Add(bcDist);
                    correlations.Add(corrCoeff);
                    enelopes.Add(obsEnv);                    
                }
            }

            if (bcDistances.Count < 1) return 0d;

            return Math.Max(1 - bcDistances.Min(), correlations.Max());
        }



        /*
        public override string GetFeatureFile(string rawFilePath, double minSearchMass = 3000, double maxSearchMass = 50000, bool scoreReport = false, bool csvOutput = false, bool tmpOutput = false)
        {
            var outTsvFilePath = Path.ChangeExtension(rawFilePath, FileExtension);
            if (File.Exists(outTsvFilePath)) return outTsvFilePath;

            var j = rawFilePath.LastIndexOf('.');
            var outPath = rawFilePath.Substring(0, j);

            var outCsvFilePath = string.Format("{0}_{1}.csv", outPath, FileExtension);
            var container = new Ms1FeatureContainer(Spectrums);

            if (File.Exists(outTsvFilePath)) return outTsvFilePath;

            MinSearchMassBin = Comparer.GetBinNumber(minSearchMass);
            MaxSearchMassBin = Comparer.GetBinNumber(maxSearchMass);

            double totalMassBin = MaxSearchMassBin - MinSearchMassBin + 1;

            Console.WriteLine("Start MS1 feature extraction.");
            if (csvOutput) Console.WriteLine("Csv Output File\t{0}", outCsvFilePath);

            var stopwatch = Stopwatch.StartNew();
            for (var binNum = MinSearchMassBin; binNum <= MaxSearchMassBin; binNum++)
            {
                var clusters = GetProbableClusters(binNum);
                container.Add(clusters);

                if (binNum > MinSearchMassBin && (binNum - MinSearchMassBin) % 3000 == 0)
                {
                    var elapsed = (stopwatch.ElapsedMilliseconds) / 1000.0d;
                    var processedBins = binNum - MinSearchMassBin;
                    var processedPercentage = ((double)processedBins / totalMassBin) * 100;
                    Console.WriteLine("Processing {0:0.0}% of mass bins ({1:0.0} Da); elapsed time = {2:0.000} sec; # of features = {3}", processedPercentage, Comparer.GetMzEnd(binNum), elapsed, container.NumberOfFeatures);
                }
            }

            Console.WriteLine("Complete MS1 feature extraction.");
            Console.WriteLine(" - Elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            Console.WriteLine(" - Number of extracted features = {0}", container.NumberOfFeatures);

            StreamWriter csvWriter = null;
            if (csvOutput)
            {
                csvWriter = new StreamWriter(outCsvFilePath);
                csvWriter.WriteLine("scan_num,charge,abundance,mz,fit,monoisotopic_mw,FeatureID");
            }

            Console.WriteLine("Start selecting mutually independent features from feature network graph");
            var connectedFeatures = container.GetAllConnectedFeatures();
            if (tmpOutput)
            {
                var tmpTsvFilePath = string.Format("{0}.tmp", outTsvFilePath);
                var tmpTsvWriter = new StreamWriter(tmpTsvFilePath);
                tmpTsvWriter.WriteLine(Ms1FeatureCluster.GetHeaderString(scoreReport));
                var clusterId = 0;
                foreach (var featureSet in connectedFeatures)
                {
                    clusterId++;
                    foreach (var feature in featureSet)
                    {
                        tmpTsvWriter.WriteLine("{0}\t{1}", clusterId, feature.GetString(scoreReport));
                    }
                }
                tmpTsvWriter.Close();
            }

            var filteredFeatures = container.GetFilteredFeatures(connectedFeatures);
            stopwatch.Stop();

            var featureId = 0;
            var tsvWriter = new StreamWriter(outTsvFilePath);
            tsvWriter.WriteLine(Ms1FeatureCluster.GetHeaderString(scoreReport));

            var tsvWriterQuant = new StreamWriter(Path.ChangeExtension(rawFilePath, QuantFileExtension));
            tsvWriterQuant.Write("FeatureId\t");
            for (var i = MinScanCharge; i <= MaxScanCharge; i++)
            {
                tsvWriterQuant.Write("{0}_BcDist\t", i);
            }
            for (var i = MinScanCharge; i <= MaxScanCharge; i++)
            {
                tsvWriterQuant.Write("{0}_Abundance", i);
                tsvWriterQuant.Write(i != MaxScanCharge ? "\t" : "\n");
            }

            foreach (var cluster in filteredFeatures)
            {
                featureId++;
                tsvWriter.WriteLine("{0}\t{1}", featureId, cluster.GetString(scoreReport));

                if (csvWriter != null)
                {
                    foreach (var envelope in cluster.Envelopes)
                    {
                        //var mostAbuIsotopeInternalIndex = cluster.IsotopeList.SortedIndexByIntensity[0];
                        var refInternalIndex = envelope.RefIsotopeInternalIndex;
                        var mostAbuPeak = envelope.Peaks[refInternalIndex];
                        var charge = envelope.Row + ChargeRange.Min;
                        var mass = Ion.GetMonoIsotopicMass(mostAbuPeak.Mz, charge, cluster.IsotopeList[refInternalIndex].Index);
                        var abundance = envelope.Abundance;
                        csvWriter.WriteLine(string.Format("{0},{1},{2},{3},{4},{5},{6}", Ms1ScanNums[envelope.Col], charge, abundance, mostAbuPeak.Mz, 1.0 - cluster.Probability, mass, featureId));
                    }
                }

                tsvWriterQuant.Write("{0}\t", featureId);
                for (var i = MinScanCharge; i <= MaxScanCharge; i++)
                {
                    var bc = (i >= cluster.MinCharge && i <= cluster.MaxCharge) ? cluster.BcDistPerCharge[i - cluster.MinCharge] : 0;
                    tsvWriterQuant.Write("{0}\t", bc);
                }
                for (var i = MinScanCharge; i <= MaxScanCharge; i++)
                {
                    var abu = (i >= cluster.MinCharge && i <= cluster.MaxCharge) ? cluster.AbundancePerCharge[i - cluster.MinCharge] : 0;
                    tsvWriterQuant.Write("{0}", abu);
                    tsvWriterQuant.Write(i != MaxScanCharge ? "\t" : "\n");
                }
            }
            tsvWriter.Close();
            tsvWriterQuant.Close();

            Console.WriteLine("Complete feature filteration");
            Console.WriteLine(" - Elapsed time = {0:0.000} sec", (stopwatch.ElapsedMilliseconds) / 1000.0d);
            Console.WriteLine(" - Number of filtered features = {0}", filteredFeatures.Count);
            Console.WriteLine(" - ProMex output: {0}", outTsvFilePath);
            if (csvOutput)
            {
                csvWriter.Close();
                Console.WriteLine(" - ProMex output in ICR2LS format: {0}", outCsvFilePath);
            }
            return outTsvFilePath;
        }
        */

        protected override double GetBcDistTh(double normalizedElutionLen)
        {
            if (QueryMass < 15000)
            {
                /*if (normalizedElutionLen < 0.005) return 0.6;
                if (normalizedElutionLen < 0.01) return 0.4;
                if (normalizedElutionLen < 0.02) return 0.2;
                return 0.1;*/

                return 0.25;
            }
            else if (QueryMass < 25000)
            {
                //if (normalizedElutionLen < 0.005) return 1.0;
                //if (normalizedElutionLen < 0.01) return 0.5;
                //if (normalizedElutionLen < 0.02) return 0.25;
                //return 0.2;

                return 0.5;
            }
            else // > 25K 
            {
                //if (normalizedElutionLen < 0.005) return 1.2;
                //if (normalizedElutionLen < 0.01) return 0.8;
                //if (normalizedElutionLen < 0.02) return 0.3;
                //return 0.2;
                return 0.8;
            }
        }

        protected override double GetCorrTh(double normalizedElutionLen)
        {
            if (QueryMass < 15000)
            {
                //if (normalizedElutionLen < 0.005) return 0.3;
                //if (normalizedElutionLen < 0.01) return 0.4;
                //if (normalizedElutionLen < 0.02) return 0.6;
                //return 0.7;
                return 0.4;
            }
            else if (QueryMass < 25000)
            {
                //if (normalizedElutionLen < 0.005) return 0;
                //if (normalizedElutionLen < 0.01) return 0.2;
                //if (normalizedElutionLen < 0.02) return 0.4;
                //return 0.6;
                return 0.2;
            }
            else // 25K
            {
                //if (normalizedElutionLen < 0.005) return -1;
                //if (normalizedElutionLen < 0.01) return 0.1;
                //if (normalizedElutionLen < 0.02) return 0.4;
                //return 0.5;
                return -1;
            }
        }

        /*
        protected override void EvaluateCluster(ref Ms1FeatureCluster cluster)
        {
            var minCol = cluster.MinCol;
            var maxCol = cluster.MaxCol;
            var minRow = cluster.MinRow;
            var maxRow = cluster.MaxRow;
            var massTol = MzTolerance.GetToleranceAsTh(cluster.RepresentativeMass);

            cluster.ClearMember();
            for (var col = minCol; col <= maxCol; col++)
            {
                for (var row = minRow; row <= maxRow; row++)
                {
                    var mass = AccurateMass[row][col];

                    if (mass > 0 && Math.Abs(cluster.RepresentativeMass - mass) < massTol)
                    {
                        var obsEnv = new ObservedEnvelope(row, col, FeatureMatrix[row][col], _isotopeList);
                        cluster.AddMember(obsEnv);
                    }
                    else
                    {
                        var observedPeaks = Spectrums[col].GetAllIsotopePeaks(cluster.RepresentativeMass, row + ChargeRange.Min, _isotopeList, MzTolerance);
                        var obsEnv = new ObservedEnvelope(row, col, observedPeaks, _isotopeList);
                        if (obsEnv.NumberOfPeaks < 3) continue;
                        cluster.AddMember(obsEnv);
                    }
                }
            }

            if (cluster.Envelopes.Count < 1) return;
            cluster.NormalizedElutionLength = (Run.GetElutionTime(cluster.MaxScanNum) -
                                               Run.GetElutionTime(cluster.MinScanNum)) /
                                              Run.GetElutionTime(Run.MaxLcScan);
            cluster.UpdateScores(Spectrums);
            if (cluster.GoodEnvelopeCount > 0)
            {
                CalculateXicCorrelationOverTimeBetweenIsotopes(cluster);
            }
        }
        */

    }
}

