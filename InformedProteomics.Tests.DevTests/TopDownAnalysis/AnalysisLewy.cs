using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.FeatureFinding.Alignment;
using InformedProteomics.FeatureFinding.SpectrumMatching;
using MathNet.Numerics.Statistics;
using NUnit.Framework;

namespace InformedProteomics.Tests.DevTests.TopDownAnalysis
{
    internal class AnalysisLewy
    {
        // Ignore Spelling: Lewy, Jungkap, prsm, promex, Desc, dest, msp, foreach, Pbf, Tsv, Lc

        public const string DataPath = @"\\Proto-5\VOrbiETD02\2014_3";
        public const string PbfPath = @"\\proto-11\MSXML_Cache\PBF_Gen_1_214\2014_3";
        public const string DatabaseFilePath = @"\\protoapps\UserData\Jungkap\Lewy\db\ID_005140_7A170668.fasta";

        public const string MsPfFolder = @"\\protoapps\UserData\Jungkap\Lewy\NewScore";
        public const string MsPfFolder2 = @"\\protoapps\UserData\Jungkap\Lewy\DMS_new\SYUA_Human_noVar";
        public const string Ms1FtFolder = @"\\protoapps\UserData\Jungkap\Lewy\DMS_new";

        public const int DATASET_COUNT = 3;
        public const int DATASET_COUNT_FULL = 51;

        public string GetDataSetNames(int index)
        {
            var i = index + 1;
            var dataName = "Lewy_intact_" + i.ToString("D2");
            return dataName;
        }

        private List<ProteinSpectrumMatch> MergePrsm(IReadOnlyCollection<ProteinSpectrumMatch> targetList)
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
                    if (prsm.SpectralEValue < ret[scan].SpectralEValue)
                    {
                        ret[scan] = prsm;
                    }
                }
            }

            return ret.Where(prsm => prsm != null).ToList();
        }

        [Test]
        [Category("Local_Testing")]
        [Category("Long_Running")]
        [Ignore("Slow")]
        public void TestFeatureAlignment()
        {
            const string outFilePath = @"\\protoapps\UserData\Jungkap\Lewy\aligned\promex_crosstab_temp.tsv";

            //CPTAC_Intact_CR32A_24Aug15_Bane_15-02-06-RZ
            var prsmReader = new ProteinSpectrumMatchReader();
            var tolerance = new Tolerance(10);
            var alignment = new LcMsFeatureAlignment(new AnalysisCompRef.CompRefFeatureComparer(tolerance));

            for (var i = 0; i < DATASET_COUNT; i++)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", PbfPath, GetDataSetNames(i));
                var mspFile = string.Format(@"{0}\{1}_IcTda.tsv", MsPfFolder, GetDataSetNames(i));
                var mspFile2 = string.Format(@"{0}\{1}_IcTda.tsv", MsPfFolder2, GetDataSetNames(i));
                var ms1FtFile = string.Format(@"{0}\{1}.ms1ft", Ms1FtFolder, GetDataSetNames(i));

                if (!File.Exists(rawFile))
                {
                    Console.WriteLine("Skipping test since file not found: " + rawFile);
                    continue;
                }

                if (!File.Exists(ms1FtFile))
                {
                    Console.WriteLine("Skipping test since file not found: " + ms1FtFile);
                    continue;
                }

                if (!File.Exists(mspFile))
                {
                    Console.WriteLine("Skipping test since file not found: " + mspFile);
                    continue;
                }

                if (!File.Exists(mspFile2))
                {
                    Console.WriteLine("Skipping test since file not found: " + mspFile2);
                    continue;
                }

                Console.WriteLine(rawFile);

                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                var prsmList1 = prsmReader.LoadIdentificationResult(mspFile, ProteinSpectrumMatch.SearchTool.MsPathFinder);
                var prsmList2 = prsmReader.LoadIdentificationResult(mspFile2, ProteinSpectrumMatch.SearchTool.MsPathFinder);
                prsmList1.AddRange(prsmList2);

                var prsmList = MergePrsm(prsmList1);
                var features = LcMsFeatureAlignment.LoadProMexResult(i, ms1FtFile, run);

                foreach (var match in prsmList)
                {
                    match.ProteinId = match.ProteinName;
                }

                // tag features by PrSMs
                foreach (var feature in features)
                {
                    //features[j].ProteinSpectrumMatches = new ProteinSpectrumMatchSet(i);
                    var massTol = tolerance.GetToleranceAsMz(feature.Mass);
                    foreach (var match in prsmList)
                    {
                        if (feature.MinScanNum < match.ScanNum && match.ScanNum < feature.MaxScanNum && Math.Abs(feature.Mass - match.Mass) < massTol)
                        {
                            feature.ProteinSpectrumMatches.Add(match);
                        }
                    }
                }

                alignment.AddDataSet(i, features, run);
            }

            if (alignment.CountDatasets == 0)
            {
                Assert.Ignore("Skipping test since input data files were found");
            }

            alignment.AlignFeatures();

            Console.WriteLine("{0} alignments ", alignment.CountAlignedFeatures);

            for (var i = 0; i < DATASET_COUNT; i++)
            {
                alignment.FillMissingFeatures(i);
                Console.WriteLine("{0} has been processed", GetDataSetNames(i));
            }

            OutputCrossTabWithId(outFilePath, alignment);
        }

        private void OutputCrossTabWithId(string outputFilePath, LcMsFeatureAlignment alignment)
        {
            using var writer = new StreamWriter(outputFilePath);

            var headerLine = new List<string>
            {
                "MonoMass",
                "MinElutionTime",
                "MaxElutionTime"
            };

            for (var i = 0; i < DATASET_COUNT; i++)
            {
                var dataName = GetDataSetNames(i);
                headerLine.Add(dataName + "_Abundance");
            }

            for (var i = 0; i < DATASET_COUNT; i++)
            {
                var dataName = GetDataSetNames(i);
                headerLine.Add(dataName + "_Ms1Score");
            }

            headerLine.Add("Pre");
            headerLine.Add("Sequence");
            headerLine.Add("Post");
            headerLine.Add("Modifications");
            headerLine.Add("ProteinName");
            headerLine.Add("ProteinDesc");
            headerLine.Add("ProteinLength");
            headerLine.Add("Start");
            headerLine.Add("End");

            for (var i = 0; i < DATASET_COUNT; i++)
            {
                var dataName = GetDataSetNames(i);
                headerLine.Add(dataName + "_SpectraCount");
            }

            writer.WriteLine(string.Join("\t", headerLine));

            var alignedFeatureList = alignment.GetAlignedFeatures();
            foreach (var features in alignedFeatureList)
            {
                var mass = features.Where(f => f != null).Select(f => f.Mass).Median();
                var minElutionTime = features.Where(f => f != null).Select(f => f.MinElutionTime).Median();
                var maxElutionTime = features.Where(f => f != null).Select(f => f.MaxElutionTime).Median();

                var dataLine = new List<string>
                {
                    PRISM.StringUtilities.DblToString(mass, 4),
                    PRISM.StringUtilities.DblToString(minElutionTime, 3),
                    PRISM.StringUtilities.DblToString(maxElutionTime, 3)
                };

                for (var i = 0; i < DATASET_COUNT; i++)
                {
                    if (features[i] == null)
                    {
                        dataLine.Add("0");
                    }
                    else
                    {
                        dataLine.Add(PRISM.StringUtilities.DblToString(features[i].Abundance, 2));
                    }
                }

                for (var i = 0; i < DATASET_COUNT; i++)
                {
                    if (features[i] == null)
                    {
                        dataLine.Add("0");
                    }
                    else
                    {
                        if (features[i].Score <= float.MinValue)
                        {
                            dataLine.Add(PRISM.StringUtilities.DblToStringScientific(float.MinValue, 2));
                        }
                        else
                        {
                            dataLine.Add(PRISM.StringUtilities.DblToString(features[i].Score, 3));
                        }
                    }
                }

                var prsm = (from f in features
                            where f?.ProteinSpectrumMatches != null && f.ProteinSpectrumMatches.Count > 0
                            select f.ProteinSpectrumMatches[0]).FirstOrDefault();

                if (prsm == null)
                {
                    for (var k = 0; k < 9; k++)
                    {
                        dataLine.Add(" ");
                    }
                }
                else
                {
                    dataLine.Add(prsm.Pre);
                    dataLine.Add(prsm.Sequence);
                    dataLine.Add(prsm.Post);
                    dataLine.Add(prsm.Modifications);
                    dataLine.Add(prsm.ProteinName);
                    dataLine.Add(prsm.ProteinDesc);
                    dataLine.Add(prsm.ProteinLength.ToString());
                    dataLine.Add(prsm.FirstResidue.ToString());
                    dataLine.Add(prsm.LastResidue.ToString());
                }

                // spectral count from ms2
                for (var i = 0; i < DATASET_COUNT; i++)
                {
                    if (features[i] == null)
                    {
                        dataLine.Add("0");
                    }
                    else
                    {
                        dataLine.Add(features[i].ProteinSpectrumMatches.Count.ToString());
                    }
                }

                writer.WriteLine(string.Join("\t", dataLine));
            }

            Console.WriteLine("Results written to " + outputFilePath);
        }

        [Test]
        [Category("Local_Testing")]
        [Ignore("Local files")]
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
        [Category("PNL_Domain")]
        public void CopyAllMSPFResult()
        {
            const string destPath = @"\\protoapps\UserData\Jungkap\Lewy\DMS_old\201502";
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
                var newMspDir = string.Empty;
                for (var j = directories.Length - 1; j >= 0; j--)
                {
                    if (directories[j].IndexOf("Auto10889", StringComparison.OrdinalIgnoreCase) > 0)
                    {
                        newMspDir = directories[j];
                        break;
                    }
                }

                if (newMspDir.Length < 1)
                {
                    Console.WriteLine(dataName + " is missing");
                    continue;
                }

                //if (i != 1 && i != 38 && i != 39 && i != 46)
                //{
                Console.WriteLine(newMspDir);
                var sourceFile = new FileInfo(string.Format(@"{0}\{1}_IcTsv.zip", newMspDir, dataName));
                var destFile = new FileInfo(string.Format(@"{0}\{1}_IcTsv.zip", destPath, dataName));

                if (destFile.Exists && destFile.Length == sourceFile.Length)
                {
                    Console.WriteLine();
                    Console.WriteLine("Skipping file copy since destination file exists and matches the source file's size: " + destFile.FullName);
                    continue;
                }

                sourceFile.CopyTo(destFile.FullName, true);

                //filePath = string.Format(@"{0}\{1}.ms1ft", newMspDir, dataName);
                //destFilePath = string.Format(@"{0}\{1}.ms1ft", destPath, dataName);
                //File.Copy(filePath, destFilePath, true);
                //}

                /*
                var oldMspDir = string.Empty;
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
                File.Copy(filePath2, destFilePath2, true);
                */
                //foreach(var directories
            }
        }

        [Test]
        [Category("PNL_Domain")]
        public void CopyAllMSPFResultFor10Reps()
        {
            const string destPath = @"\\protoapps\UserData\Jungkap\CPTAC_10reps\DMS_old";
            //var oldDestPath = @"\\protoapps\UserData\Jungkap\Lewy\DMS_old";
            const string dmsPath = @"\\Proto-5\VOrbiETD02\2015_1";

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
                var newMspDir = string.Empty;
                for (var j = directories.Length - 1; j >= 0; j--)
                {
                    if (directories[j].IndexOf("MSP201502", StringComparison.OrdinalIgnoreCase) > 0)
                    {
                        newMspDir = directories[j];
                        break;
                    }
                }

                if (newMspDir.Length < 1)
                {
                    Console.WriteLine(i + " is missing");
                    continue;
                }

                Console.WriteLine(newMspDir);

                var sourceFile = new FileInfo(string.Format(@"{0}\{1}_IcTsv.zip", newMspDir, dataName));
                var destFile = new FileInfo(string.Format(@"{0}\{1}_IcTsv.zip", destPath, dataName));

                if (destFile.Exists && destFile.Length == sourceFile.Length)
                {
                    Console.WriteLine();
                    Console.WriteLine("Skipping file copy since destination file exists and matches the source file's size: " + destFile.FullName);
                    continue;
                }

                sourceFile.CopyTo(destFile.FullName, true);
            }
        }
        /*
        [Test]
        [Category("Local_Testing")]
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
                var theoreticalEnvelope = cluster.TheoreticalEnvelope;

                var summedEnvelope = new double[theoreticalEnvelope.Size];
                var bestCorr = -1d;
                var bestCharge = 0;
                double[] bestChargeEnvelope = null;

                for (var c = 0; c < cluster.Envelopes.Length; c++)
                {
                    var summedAtCharge = new double[theoreticalEnvelope.Size];
                    for (var t = 0; t < cluster.Envelopes[c].Length; t++)
                    {
                        if (cluster.Envelopes[c][t] == null) continue;
                        cluster.Envelopes[c][t].Peaks.SumEnvelopeTo(summedAtCharge);
                        cluster.Envelopes[c][t].Peaks.SumEnvelopeTo(summedEnvelope);
                    }
                    var corr1 = theoreticalEnvelope.GetPearsonCorrelation(summedAtCharge);

                    if (corr1 > bestCorr)
                    {
                        bestCorr = corr1;
                        bestChargeEnvelope = summedAtCharge;
                        bestCharge = c + 9;
                    }
                }
                var corr2 = theoreticalEnvelope.GetPearsonCorrelation(summedEnvelope);

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
                Console.Write(string.Join("\t", bestChargeEnvelope));
                Console.Write("\n");
            }
        }*/
    }
}
