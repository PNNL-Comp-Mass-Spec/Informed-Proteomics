using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.FeatureFinding;
using InformedProteomics.FeatureFinding.Alignment;
using InformedProteomics.FeatureFinding.FeatureDetection;
using InformedProteomics.FeatureFinding.Scoring;
using InformedProteomics.FeatureFinding.SpectrumMatching;
using NUnit.Framework;

namespace InformedProteomics.Tests.DevTests.TopDownAnalysis
{
    internal class AnalysisCompRefKelleherData
    {
        private List<string> GetDataSetNamesStudy3(int fraction1, int fraction2)
        {
            var i = fraction1;
            var j = fraction2;
            var ret = new List<string>();
            for (var k = 1; k <= 5; k++)
            {
                var dataName1 = string.Format("WHIM16_GFrep{0:D2}_GFfrac{1:D2}_inj{2:D2}", i, j, k); // Study-3
                var dataName2 = string.Format("WHIM2_GFrep{0:D2}_GFfrac{1:D2}_inj{2:D2}", i, j, k);  // Study-3
                ret.Add(dataName1);
                ret.Add(dataName2);
            }
            ret.Sort();
            return ret;
        }

        public const string DataPath = @"\\proto-4\NU_QExactive\2015_4";
        public const string PbfPath = @"\\proto-11\MSXML_Cache\PBF_Gen_1_214\2015_4";

        [Test]
        [Category("Local_Testing")]
        [Ignore("Local files")]
        public void CopyAllMSPFResult()
        {
            const string destPath = @"D:\MassSpecFiles\CompRef_Kelleher\Study2";

            if (!Directory.Exists(DataPath))
            {
                Assert.Ignore("Skipped since source directory not found: " + DataPath);
            }

            if (!Directory.Exists(destPath))
            {
                Assert.Ignore("Skipped since destination directory not found: " + destPath);
            }

            var filesCopied = 0;

            for (var i = 1; i <= 1; i++)
            {
                for (var j = 1; j <= 3; j++)
                {
                    for (var k = 1; k <= 6; k++)
                    {
                        //var dataName = string.Format(@"WHIM16_GFrep{0}_GFfrac{1}_inj{2}", i.ToString("D2"), j.ToString("D2"), k.ToString("D2")); // Study-3
                        var dataName = string.Format("RP4H_P33_WHIM16_biorep{0}_techrep{1}", j, k); // Study-2
                        var dataDir = string.Format(@"{0}\{1}", DataPath, dataName);

                        if (!Directory.Exists(dataDir))
                        {
                            Console.WriteLine("Skipping dataset since directory not found: " + dataDir);
                            continue;
                        }

                        var directories = Directory.GetDirectories(dataDir);

                        foreach (var dir in directories)
                        {
                            if (dir.Contains("MSP"))
                            {
                                var icTsvFile = string.Format(@"{0}\{1}_IcTsv.zip", dir, dataName);

                                if (File.Exists(icTsvFile))
                                {
                                    var destFilePath = string.Format(@"{0}\{1}_IcTsv.zip", destPath, dataName);
                                    File.Copy(icTsvFile, destFilePath);
                                    filesCopied++;
                                }

                                var featuresFile = string.Format(@"{0}\{1}.ms1ft", dir, dataName);
                                if (File.Exists(featuresFile))
                                {
                                    var destFilePath = string.Format(@"{0}\{1}.ms1ft", destPath, dataName);
                                    File.Copy(featuresFile, destFilePath);
                                    filesCopied++;
                                }
                            }
                        }
                    }
                }
            }

            if (filesCopied == 0)
            {
                Assert.Ignore("Skipping since data files were not found");
            }
        }

        [Test]
        [Category("Local_Testing")]
        [Ignore("Local files")]
        public void FindMissingLcMsFeatures()
        {
            const string mspfFolder = @"D:\MassSpecFiles\CompRef_Kelleher\Study3";
            const string ms1ftFolder = @"D:\MassSpecFiles\CompRef_Kelleher\Study3";

            const int Nfraction1 = 3;
            const int Nfraction2 = 5;
            var filesProcessed = 0;

            for (var frac1 = 1; frac1 <= Nfraction1; frac1++)
            {
                for (var frac2 = 1; frac2 <= Nfraction2; frac2++)
                {
                    var datasets = GetDataSetNamesStudy3(frac1, frac2);
                    //var outFilePath = string.Format(@"D:\MassSpecFiles\CompRef_Kelleher\study3_GFrep{0}_Gfrac{1}.tsv", frac1.ToString("D2"), frac2.ToString("D2"));
                    var nDataset = datasets.Count;
                    var prsmReader = new ProteinSpectrumMatchReader();
                    var tolerance = new Tolerance(12);

                    for (var i = 0; i < nDataset; i++)
                    {
                        var rawFile = string.Format(@"{0}\{1}.pbf", PbfPath, datasets[i]);
                        var mspFile = string.Format(@"{0}\{1}_IcTda.tsv", mspfFolder, datasets[i]);
                        var ms1FtFile = string.Format(@"{0}\{1}.ms1ft", ms1ftFolder, datasets[i]);
                        var outPath = string.Format(@"{0}\{1}.seqtag.ms1ft", ms1ftFolder, datasets[i]);

                        if (!File.Exists(rawFile))
                        {
                            Console.WriteLine("Skipping dataset since file not found: " + rawFile);
                            continue;
                        }

                        if (!File.Exists(ms1FtFile))
                        {
                            Console.WriteLine("Skipping dataset since file not found: " + ms1FtFile);
                            continue;
                        }

                        if (!File.Exists(mspFile))
                        {
                            Console.WriteLine("Skipping dataset since file not found: " + mspFile);
                            continue;
                        }

                        if (File.Exists(outPath))
                        {
                            Console.WriteLine("Skipping dataset since results file already exists: " + outPath);
                            continue;
                        }

                        filesProcessed++;

                        var run = PbfLcMsRun.GetLcMsRun(rawFile);
                        var features = LcMsFeatureAlignment.LoadProMexResult(i, ms1FtFile, run);
                        var prsmList = prsmReader.LoadIdentificationResult(mspFile, ProteinSpectrumMatch.SearchTool.MsPathFinder);
                        var prsmFeatureMatch = new bool[prsmList.Count];

                        foreach (var feature in features)
                        {
                            // feature.ProteinSpectrumMatches = new ProteinSpectrumMatchSet(i);
                            var massTol = tolerance.GetToleranceAsMz(feature.Mass);
                            for (var k = 0; k < prsmList.Count; k++)
                            {
                                var match = prsmList[k];
                                if (feature.MinScanNum < match.ScanNum && match.ScanNum < feature.MaxScanNum && Math.Abs(feature.Mass - match.Mass) < massTol)
                                {
                                    feature.ProteinSpectrumMatches.Add(match);
                                    prsmFeatureMatch[k] = true;
                                }
                            }
                        }

                        var missingPrsm = new List<ProteinSpectrumMatch>();
                        for (var k = 0; k < prsmList.Count; k++)
                        {
                            if (!prsmFeatureMatch[k])
                            {
                                missingPrsm.Add(prsmList[k]);
                            }
                        }

                        FeatureFind(missingPrsm, run, outPath);
                        Console.WriteLine(outPath);
                    }
                }
            }

            if (filesProcessed == 0)
            {
                Assert.Ignore("Skipping since data files were not found");
            }
        }

        private void FeatureFind(List<ProteinSpectrumMatch> prsms, LcMsRun run, string outTsvFilePath)
        {
            var featureFinder = new LcMsPeakMatrix(run, new LcMsFeatureLikelihood());
            // write result files
            var tsvWriter = new StreamWriter(outTsvFilePath);
            tsvWriter.WriteLine(LcMsFeatureFinderLauncher.GetHeaderString());

            var featureId = 1;
            foreach (var match in prsms)
            {
                var minScan = run.GetPrevScanNum(match.ScanNum, 1);
                var maxScan = run.GetNextScanNum(match.ScanNum, 1);
                var feature = featureFinder.GetLcMsPeakCluster(match.Mass, match.Charge, minScan, maxScan);

                if (feature == null)
                {
                    continue;
                }

                tsvWriter.WriteLine("{0}\t{1}", featureId, LcMsFeatureFinderLauncher.GetString(feature));
                featureId++;
            }

            tsvWriter.Close();
        }

        [Test]
        [Category("Local_Testing")]
        [Ignore("Local files")]
        public void AnalysisStudy3()
        {
            const string mspfFolder = @"D:\MassSpecFiles\CompRef_Kelleher\Study3";
            const string ms1ftFolder = @"D:\MassSpecFiles\CompRef_Kelleher\Study3";

            const int Nfraction1 = 3;
            const int Nfraction2 = 5;

            var filesProcessed = 0;

            for (var frac1 = 1; frac1 <= Nfraction1; frac1++)
            {
                for (var frac2 = 1; frac2 <= Nfraction2; frac2++)
                {
                    var datasets = GetDataSetNamesStudy3(frac1, frac2);
                    var outFilePath = string.Format(@"D:\MassSpecFiles\CompRef_Kelleher\study3_GFrep{0:D2}_Gfrac{1:D2}.tsv", frac1, frac2);
                    if (!File.Exists(outFilePath))
                    {
                        continue;
                    }

                    AlignFeatures(datasets, mspfFolder, ms1ftFolder, outFilePath);
                    Console.WriteLine("############ {0}-{1} has been completed", frac1, frac2);

                    filesProcessed++;
                }
            }

            if (filesProcessed == 0)
            {
                Assert.Ignore("Skipped since data files not found");
            }
        }

        private void AlignFeatures(List<string> datasets, string mspfFolder, string ms1ftFolder, string outFilePath)
        {
            var nDataset = datasets.Count;
            var prsmReader = new ProteinSpectrumMatchReader();
            var tolerance = new Tolerance(12);
            var alignment = new LcMsFeatureAlignment(new AnalysisCompRef.CompRefFeatureComparer(tolerance));
            for (var i = 0; i < nDataset; i++)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", PbfPath, datasets[i]);
                var mspFile = string.Format(@"{0}\{1}_IcTda.tsv", mspfFolder, datasets[i]);
                var ms1FtFile = string.Format(@"{0}\{1}.ms1ft", ms1ftFolder, datasets[i]);
                var ms1FtFile2 = string.Format(@"{0}\{1}.seqtag.ms1ft", ms1ftFolder, datasets[i]);

                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                var features = LcMsFeatureAlignment.LoadProMexResult(i, ms1FtFile, run);
                var features2 = LcMsFeatureAlignment.LoadProMexResult(i, ms1FtFile2, run);
                features.AddRange(features2);

                if (File.Exists(mspFile))
                {
                    var prsmList = prsmReader.LoadIdentificationResult(mspFile, ProteinSpectrumMatch.SearchTool.MsPathFinder);
                    //var prsmFeatureMatch = new bool[prsmList.Count];

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
                                //prsmFeatureMatch[k] = true;
                            }
                        }
                    }
                }

                alignment.AddDataSet(i, features, run);
            }

            alignment.AlignFeatures();

            Console.WriteLine("{0} alignments ", alignment.CountAlignedFeatures);

            for (var i = 0; i < nDataset; i++)
            {
                alignment.FillMissingFeatures(i);
                Console.WriteLine("{0} has been processed", datasets[i]);
            }

            AnalysisCompRef.OutputCrossTabWithId(outFilePath, alignment, datasets);
        }
    }
}
