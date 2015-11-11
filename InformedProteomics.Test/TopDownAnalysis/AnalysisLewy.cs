using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using MathNet.Numerics.Statistics;
using NUnit.Framework;



namespace InformedProteomics.Test.TopDownAnalysis
{
    class AnalysisLewy
    {
        public const string DataPath = @"\\Proto-5\VOrbiETD02\2014_3";
        public const string PbfPath = @"\\proto-11\MSXML_Cache\PBF_Gen_1_214\2014_3";
        public const string DatabaseFilePath = @"\\protoapps\UserData\Jungkap\Lewy\db\ID_005140_7A170668.fasta";

        public const string MsPfFolder = @"\\protoapps\UserData\Jungkap\Lewy\NewScore";
        public const string MsPfFolder2 = @"\\protoapps\UserData\Jungkap\Lewy\DMS_new\SYUA_Human_noVar";
        public const string Ms1FtFolder = @"\\protoapps\UserData\Jungkap\Lewy\DMS_new";

        public const int NdataSet = 51;

        public string GetDataSetNames(int index)
        {
            var i = index + 1;
            var dataName = "Lewy_intact_" + i.ToString("D2");
            return dataName;
        }

        private List<ProteinSpectrumMatch> MergePrsm(List<ProteinSpectrumMatch> targetList)
        {
            //var sortedList = targetList.OrderBy(prsm => prsm.ScanNum).ToList();
            //var minScan = sortedList.First().ScanNum;
            var maxScan = targetList.Max(prsm => prsm.ScanNum);

            var ret = new ProteinSpectrumMatch[maxScan + 1];
            foreach (var prsm in targetList)
            {
                var scan = prsm.ScanNum;

                if (ret[scan] == null)
                {
                    ret[scan] = prsm;
                }
                else
                {
                    if (prsm.SpectralEvalue < ret[scan].SpectralEvalue)
                        ret[scan] = prsm;
                }
            }

            return ret.Where(prsm => prsm != null).ToList();
        }


        [Test]
        public void TestFeatureAlignment()
        {
            const string outFilePath = @"\\protoapps\UserData\Jungkap\Lewy\aligned\promex_crosstab_temp.tsv";
            
            
            //CPTAC_Intact_CR32A_24Aug15_Bane_15-02-06-RZ
            var prsmReader = new ProteinSpectrumMatchReader();
            var tolerance = new Tolerance(10);
            var alignment = new LcMsFeatureAlignment(new AnalysisCompRef.CompRefFeatureComparer(tolerance));

            for (var i = 0; i < NdataSet; i++)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", PbfPath, GetDataSetNames(i));
                var mspFile = string.Format(@"{0}\{1}_IcTda.tsv", MsPfFolder, GetDataSetNames(i));
                var mspFile2 = string.Format(@"{0}\{1}_IcTda.tsv", MsPfFolder2, GetDataSetNames(i));
                var ms1FtFile = string.Format(@"{0}\{1}.ms1ft", Ms1FtFolder, GetDataSetNames(i));
                Console.WriteLine(rawFile);
                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                var prsmList1 = prsmReader.LoadIdentificationResult(mspFile, ProteinSpectrumMatch.SearchTool.MsPathFinder);
                var prsmList2 = prsmReader.LoadIdentificationResult(mspFile2, ProteinSpectrumMatch.SearchTool.MsPathFinder);
                prsmList1.AddRange(prsmList2);
                
                var prsmList = MergePrsm(prsmList1);
                var features = LcMsFeatureAlignment.LoadProMexResult(i, ms1FtFile, run);

                for (var j = 0; j < prsmList.Count; j++)
                {
                    var match = prsmList[j];
                    match.ProteinId = match.ProteinName;
                }

                // tag features by PrSMs
                for (var j = 0; j < features.Count; j++)
                {
                    //features[j].ProteinSpectrumMatches = new ProteinSpectrumMatchSet(i);
                    var massTol = tolerance.GetToleranceAsTh(features[j].Mass);
                    foreach (var match in prsmList)
                    {
                        if (features[j].MinScanNum < match.ScanNum && match.ScanNum < features[j].MaxScanNum && Math.Abs(features[j].Mass - match.Mass) < massTol)
                        {
                            features[j].ProteinSpectrumMatches.Add(match);
                        }
                    }
                }

                alignment.AddDataSet(i, features, run);
            }

            alignment.AlignFeatures();

            Console.WriteLine("{0} alignments ", alignment.CountAlignedFeatures);
            
            for (var i = 0; i < NdataSet; i++)
            {
                alignment.FillMissingFeatures(i);
                Console.WriteLine("{0} has been processed", GetDataSetNames(i));
            }
            
            OutputCrossTabWithId(outFilePath, alignment);
        }

        public void OutputCrossTabWithId(string outputFilePath, LcMsFeatureAlignment alignment)
        {

            var writer = new StreamWriter(outputFilePath);

            writer.Write("MonoMass");
            writer.Write("\t");
            writer.Write("MinElutionTime");
            writer.Write("\t");
            writer.Write("MaxElutionTime");

            for (var i = 0; i < NdataSet; i++)
            {
                var dataName = GetDataSetNames(i);
                writer.Write("\t");
                writer.Write(dataName + "_Abundance");
            }

            for (var i = 0; i < NdataSet; i++)
            {
                var dataName = GetDataSetNames(i);
                writer.Write("\t");
                writer.Write(dataName + "_Ms1Score");
            }

            writer.Write("\t");
            writer.Write("Pre");
            writer.Write("\t");
            writer.Write("Sequence");
            writer.Write("\t");
            writer.Write("Post");
            writer.Write("\t");
            writer.Write("Modifications");
            writer.Write("\t");
            writer.Write("ProteinName");
            writer.Write("\t");
            writer.Write("ProteinDesc");
            writer.Write("\t");
            writer.Write("ProteinLength");
            writer.Write("\t");
            writer.Write("Start");
            writer.Write("\t");
            writer.Write("End");

            for (var i = 0; i < NdataSet; i++)
            {
                var dataName = GetDataSetNames(i); 
                writer.Write("\t");
                writer.Write(dataName + "_SpectraCount");
            }
            writer.Write("\n");

            var alignedFeatureList = alignment.GetAlignedFeatures();
            for (var j = 0; j < alignedFeatureList.Count; j++)
            {
                var features = alignedFeatureList[j];
                var mass = features.Where(f => f != null).Select(f => f.Mass).Median();
                var minElutionTime = features.Where(f => f != null).Select(f => f.MinElutionTime).Median();
                var maxElutionTime = features.Where(f => f != null).Select(f => f.MaxElutionTime).Median();
                writer.Write(mass);
                writer.Write("\t");
                writer.Write(minElutionTime);
                writer.Write("\t");
                writer.Write(maxElutionTime);

                for (var i = 0; i < NdataSet; i++)
                {
                    writer.Write("\t");
                    writer.Write(features[i] == null ? 0 : features[i].Abundance);
                }

                for (var i = 0; i < NdataSet; i++)
                {
                    writer.Write("\t");
                    writer.Write(features[i] == null ? 0 : features[i].Score);
                }

                var prsm = (from f in features
                            where f != null && f.ProteinSpectrumMatches != null && f.ProteinSpectrumMatches.Count > 0
                            select f.ProteinSpectrumMatches[0]).FirstOrDefault();

                if (prsm == null)
                {
                    for (var k = 0; k < 9; k++)
                    {
                        writer.Write("\t");
                        writer.Write(" ");
                    }
                }
                else
                {
                    writer.Write("\t");
                    writer.Write(prsm.Pre);
                    writer.Write("\t");
                    writer.Write(prsm.Sequence);
                    writer.Write("\t");
                    writer.Write(prsm.Post);
                    writer.Write("\t");
                    writer.Write(prsm.Modifications);
                    writer.Write("\t");
                    writer.Write(prsm.ProteinName);
                    writer.Write("\t");
                    writer.Write(prsm.ProteinDesc);
                    writer.Write("\t");
                    writer.Write(prsm.ProteinLength);
                    writer.Write("\t");
                    writer.Write(prsm.FirstResidue);
                    writer.Write("\t");
                    writer.Write(prsm.LastResidue);
                }

                // spectral count from ms2
                for (var i = 0; i < NdataSet; i++)
                {
                    writer.Write("\t");
                    writer.Write(features[i] == null ? 0 : features[i].ProteinSpectrumMatches.Count);
                }
                writer.Write("\n");
            }

            writer.Close();
        }


        [Test]
        public void CountTagMatches()
        {
            for (var i = 1; i < 52; i++)
            {
                var dataName = "Lewy_intact_" + i.ToString("D2");
                var filePath = string.Format(@"{0}\{1}.pbf", PbfPath, dataName);
                var run = PbfLcMsRun.GetLcMsRun(filePath);
                var scans = run.GetScanNumbers(2);
                Console.WriteLine(scans.Count);
            }
        }

        [Test]
        public void CopyAllMSPFResult()
        {
            var destPath = @"\\protoapps\UserData\Jungkap\Lewy\DMS_old\201502";
            //var oldDestPath = @"\\protoapps\UserData\Jungkap\Lewy\DMS_old";

            for (var i = 1; i < 52; i++)
            {
                //if (i == 1 || i == 38 || i == 39 || i == 46) continue;
                var dataName = "Lewy_intact_" + i.ToString("D2");
                var dataDir = string.Format(@"{0}\{1}", DataPath, dataName);

                var directories = Directory.GetDirectories(dataDir);
                var directoryTimes = new List<DateTime>();

                foreach (var d in directories)
                {
                    directoryTimes.Add(Directory.GetCreationTime(d));
                }

                Array.Sort(directoryTimes.ToArray(), directories);
                var newMspDir = "";
                for (var j = directories.Length - 1; j >= 0; j--)
                {
                    if (directories[j].IndexOf("Auto10889") > 0)
                    {
                        newMspDir = directories[j];
                        break;
                    }
                }

                if (newMspDir.Length < 1)
                {
                    Console.WriteLine(dataName + " is mssing");
                    continue;
                }

                //if (i != 1 && i != 38 && i != 39 && i != 46)
                //{
                    Console.WriteLine(newMspDir);
                    var filePath = string.Format(@"{0}\{1}_IcTsv.zip", newMspDir, dataName);
                    var destFilePath = string.Format(@"{0}\{1}_IcTsv.zip", destPath, dataName);
                    File.Copy(filePath, destFilePath);

                    //filePath = string.Format(@"{0}\{1}.ms1ft", newMspDir, dataName);
                    //destFilePath = string.Format(@"{0}\{1}.ms1ft", destPath, dataName);
                    //File.Copy(filePath, destFilePath);
                //}

                /*
                var oldMspDir = "";
                for (var j = 0; j < directories.Length; j++)
                {
                    if (directories[j].IndexOf("MSP2014") > 0)
                    {
                        oldMspDir = directories[j];
                        break;
                    }
                }
                var filePath2 = string.Format(@"{0}\{1}_IcTsv.zip", oldMspDir, dataName);
                var destFilePath2 = string.Format(@"{0}\{1}_IcTsv.zip", oldDestPath, dataName);
                File.Copy(filePath2, destFilePath2);
                */
                //foreach(var directories
            }

        }

        [Test]
        public void CopyAllMSPFResultFor10Reps()
        {
            var destPath = @"\\protoapps\UserData\Jungkap\CPTAC_10reps\DMS_old";
            var oldDestPath = @"\\protoapps\UserData\Jungkap\Lewy\DMS_old";
            var dmsPath = @"\\Proto-5\VOrbiETD02\2015_1";

            for (var i = 1; i < 11; i++)
            {
                var dataName = string.Format("CPTAC_Intact_rep{0}_15Jan15_Bane_C2-14-08-02RZ", i);
                var dataDir = string.Format(@"{0}\{1}", dmsPath, dataName);

                var directories = Directory.GetDirectories(dataDir);
                var directoryTimes = new List<DateTime>();

                foreach (var d in directories)
                {
                    directoryTimes.Add(Directory.GetCreationTime(d));
                }

                Array.Sort(directoryTimes.ToArray(), directories);
                var newMspDir = "";
                for (var j = directories.Length - 1; j >= 0; j--)
                {
                    if (directories[j].IndexOf("MSP201502") > 0)
                    {
                        newMspDir = directories[j];
                        break;
                    }
                }

                if (newMspDir.Length < 1)
                {
                    Console.WriteLine(i + " is mssing");
                    continue;
                }

                Console.WriteLine(newMspDir);
                var filePath = string.Format(@"{0}\{1}_IcTsv.zip", newMspDir, dataName);
                var destFilePath = string.Format(@"{0}\{1}_IcTsv.zip", destPath, dataName);
                File.Copy(filePath, destFilePath);

                
            }

        }
        /*
        [Test]
        public void TestFindLowAbundanceFeature()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);
            var i = 27;
            var rawFile = string.Format(@"{0}\{1}.pbf", PbfPath, GetDataSetNames(i));

            if (!File.Exists(rawFile))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFile);
            }

            // var outTsvFilePath = MassSpecDataReaderFactory.ChangeExtension(rawFile, "ms1ft");
            //var scoreDataPath = @"D:\MassSpecFiles\training";
            var scorer = new LcMsFeatureLikelihood();
            var stopwatch = Stopwatch.StartNew();
            Console.WriteLine(@"Start loading MS1 data from {0}", rawFile);
            var run = PbfLcMsRun.GetLcMsRun(rawFile);
            var featureFinder = new LcMsPeakMatrix(run, scorer);

            var comparer = new MzComparerWithBinning(28);
            
            var minSearchMassBin = comparer.GetBinNumber(14458);
            var maxSearchMassBin = comparer.GetBinNumber(14584);
            double totalMassBin = maxSearchMassBin - minSearchMassBin + 1;

            Console.WriteLine(@"Start MS1 feature extraction.");

            stopwatch.Restart();
            for (var binNum = minSearchMassBin; binNum <= maxSearchMassBin; binNum++)
            {
                //var clusters = featureFinder.FindFeatures(binNum);
                var massCenter = comparer.GetMzAverage(binNum);
                var cluster = featureFinder.GetLcMsPeakClusterWithoutScoring(massCenter, 9, 22, 7885, 8057);
                var theoEnv = cluster.TheoreticalEnvelope;

                var summedEnvelope = new double[theoEnv.Size];
                var bestCorr = -1d;
                var bestCharge = 0;
                double[] bestChargeEnv = null;

                for (var c = 0; c < cluster.Envelopes.Length; c++)
                {
                    var summedAtCharge = new double[theoEnv.Size];
                    for (var t = 0; t < cluster.Envelopes[c].Length; t++)
                    {
                        if (cluster.Envelopes[c][t] == null) continue;
                        cluster.Envelopes[c][t].Peaks.SumEnvelopeTo(summedAtCharge);
                        cluster.Envelopes[c][t].Peaks.SumEnvelopeTo(summedEnvelope);
                    }
                    var corr1 = theoEnv.GetPearsonCorrelation(summedAtCharge);

                    if (corr1 > bestCorr)
                    {
                        bestCorr = corr1;
                        bestChargeEnv = summedAtCharge;
                        bestCharge = c + 9;
                    }
                }
                var corr2 = theoEnv.GetPearsonCorrelation(summedEnvelope);

                Console.Write(massCenter);
                Console.Write("\t");
                Console.Write(corr2);
                Console.Write("\t");
                Console.Write(bestCorr);
                Console.Write("\t");
                Console.Write(bestCharge);
                Console.Write("\t");
                Console.Write(string.Join("\t", summedEnvelope));
                Console.Write("\t");
                Console.Write(string.Join("\t", bestChargeEnv));
                Console.Write("\n");
            }            
        }*/
    }
}
