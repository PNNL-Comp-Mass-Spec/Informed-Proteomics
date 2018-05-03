using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.FeatureFinding.Alignment;
using InformedProteomics.FeatureFinding.Data;
using InformedProteomics.FeatureFinding.SpectrumMatching;
using InformedProteomics.FeatureFinding.Util;
using MathNet.Numerics.Statistics;
using NUnit.Framework;

namespace InformedProteomics.Tests.DevTests.TopDownAnalysis
{
    public class AnalysisCompRef
    {
        //public const string OutFilePath = @"\\protoapps\UserData\Jungkap\CompRef\aligned\aligned_features.tsv";
        public const string Ms1FtFolder = @"\\protoapps\UserData\Jungkap\CompRef\ms1ft";
        public const string MsPfFolder = @"\\protoapps\UserData\Jungkap\CompRef\MSPF";
        public const string RawFolder = @"\\protoapps\UserData\Jungkap\CompRef\raw";

        internal class CompRefFeatureComparer : INodeComparer<LcMsFeature>
        {
            public CompRefFeatureComparer(Tolerance tolerance = null)
            {
                _tolerance = tolerance ?? new Tolerance(10);
            }

            public bool SameCluster(LcMsFeature f1, LcMsFeature f2)
            {
                if (f1.DataSetId == f2.DataSetId) return false;
                // tolerant in mass dimension?

                var massTol = Math.Min(_tolerance.GetToleranceAsMz(f1.Mass), _tolerance.GetToleranceAsMz(f2.Mass));
                if (Math.Abs(f1.Mass - f2.Mass) > massTol) return false;

                //if (!f1.CoElutedByNet(f2, 0.004)) return false; //e.g) 200*0.001 = 0.2 min = 30 sec
                if (f1.ProteinSpectrumMatches.Count > 0 && f2.ProteinSpectrumMatches.Count > 0)
                {
                    return f1.CoElutedByNet(f2, 0.03) && f1.ProteinSpectrumMatches.ShareProteinId(f2.ProteinSpectrumMatches);
                }

                return f1.CoElutedByNet(f2, 0.01);
            }

            private readonly Tolerance _tolerance;
        }

        [Test]
        [Category("Local_Testing")]
        [Ignore("Local files")]
        public void TestIMERFeatureAlignment()
        {
            const string outFilePath = @"D:\MassSpecFiles\IMER\promex_crosstab.tsv";
            const string rawFolder = @"D:\MassSpecFiles\IMER";
            var runLabels = new[] { "1", "2", "3", "4", "5", "6" };

            var nDataset = runLabels.Length;
            //CPTAC_Intact_CR32A_24Aug15_Bane_15-02-06-RZ
            var prsmReader = new ProteinSpectrumMatchReader();
            var tolerance = new Tolerance(10);
            var alignment = new LcMsFeatureAlignment(new CompRefFeatureComparer(tolerance));

            for (var i = 0; i < nDataset; i++)
            {
                var k = runLabels[i].Equals("2") || runLabels[i].Equals("3") ? 14 : 13;
                var rawFile = string.Format(@"{0}\Diabetes_iPSC_Beta_{1}_IMER_{2}May14_Alder_14-01-33.pbf", rawFolder, runLabels[i], k);
                var mspFile = string.Format(@"{0}\Diabetes_iPSC_Beta_{1}_IMER_{2}May14_Alder_14-01-33_msgfdb_syn.txt", rawFolder, runLabels[i], k);
                var ms1FtFile = string.Format(@"{0}\Diabetes_iPSC_Beta_{1}_IMER_{2}May14_Alder_14-01-33.ms1ft", rawFolder, runLabels[i], k);

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

                Console.WriteLine(rawFile);

                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                var features = LcMsFeatureAlignment.LoadProMexResult(i, ms1FtFile, run, 500, 15000);

                if (File.Exists(mspFile))
                {
                    var prsmList = prsmReader.LoadIdentificationResult(mspFile, ProteinSpectrumMatch.SearchTool.MsGfPlus);

                    foreach (var match in prsmList)
                    {
                        match.ProteinId = match.ProteinName;
                    }

                    // tag features by PrSMs
                    foreach (var item in features)
                    {
                        // item.ProteinSpectrumMatches = new ProteinSpectrumMatchSet(i);
                        var massTol = tolerance.GetToleranceAsMz(item.Mass);
                        foreach (var match in prsmList)
                        {
                            if (item.MinScanNum < match.ScanNum && match.ScanNum < item.MaxScanNum && Math.Abs(item.Mass - match.Mass) < massTol)
                            {
                                item.ProteinSpectrumMatches.Add(match);
                            }
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

            for (var i = 0; i < nDataset; i++)
            {
                alignment.FillMissingFeatures(i);
                Console.WriteLine("{0} has been processed", runLabels[i]);
            }

            OutputCrossTabWithId(outFilePath, alignment, runLabels);
        }

        [Test]
        [Category("PNL_Domain")]
        [Category("Long_Running")]
        [TestCase("32A, 32B, 33A")]
        [Ignore("Slow")]
        // Long test: [TestCase("32A, 32B, 32C, 32D, 32E, 32F, 32G, 33A, 33B, 33C, 33D, 33E, 33F, 33G")]
        public void TestFeatureAlignment(string runLabelList)
        {
            const string outFilePath = @"\\protoapps\UserData\Jungkap\CompRef\aligned\promex_crosstab_temp.tsv";

            var runLabels = new List<string>();
            foreach (var item in runLabelList.Split(','))
            {
                runLabels.Add(item.Trim());
            }

            var datasetCount = runLabels.Count;

            var prsmReader = new ProteinSpectrumMatchReader();
            var tolerance = new Tolerance(10);
            var alignment = new LcMsFeatureAlignment(new CompRefFeatureComparer(tolerance));

            var datasetId = 0;

            foreach (var runLabel in runLabels)
            {
                var rawFile = string.Format(@"{0}\CPTAC_Intact_CR{1}_24Aug15_Bane_15-02-06-RZ.pbf", RawFolder, runLabel);
                var mspFile = string.Format(@"{0}\CPTAC_Intact_CR{1}_24Aug15_Bane_15-02-06-RZ_IcTda.tsv", MsPfFolder, runLabel);
                var ms1FtFile = string.Format(@"{0}\CPTAC_Intact_CR{1}_24Aug15_Bane_15-02-06-RZ.ms1ft", Ms1FtFolder, runLabel);

                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                var features = LcMsFeatureAlignment.LoadProMexResult(datasetId, ms1FtFile, run);

                if (File.Exists(mspFile))
                {
                    var prsmList = prsmReader.LoadIdentificationResult(mspFile, ProteinSpectrumMatch.SearchTool.MsPathFinder);

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
                }

                alignment.AddDataSet(datasetId, features, run);
                datasetId++;
            }

            alignment.AlignFeatures();

            Console.WriteLine("{0} alignments ", alignment.CountAlignedFeatures);

            for (var i = 0; i < datasetCount; i++)
            {
                alignment.FillMissingFeatures(i);
                Console.WriteLine("{0} has been processed", runLabels[i]);
            }

            OutputCrossTabWithId(outFilePath, alignment, runLabels);
        }

        private void OutputCrossTabWithId(string outputFilePath, LcMsFeatureAlignment alignment, string[] runLabels)
        {
            OutputCrossTabWithId(outputFilePath, alignment, runLabels.ToList());
        }

        public static void OutputCrossTabWithId(string outputFilePath, LcMsFeatureAlignment alignment, List<string> runLabels)
        {
            var datasetCount = runLabels.Count;

            using (var writer = new StreamWriter(outputFilePath))
            {
                var headerLine = new List<string> {
                    "MonoMass",
                    "MinElutionTime",
                    "MaxElutionTime" };

                foreach (var dataName in runLabels)
                {
                    headerLine.Add(dataName + "_Abundance");
                }

                foreach (var dataName in runLabels)
                {
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

                foreach (var dataName in runLabels)
                {
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

                    for (var i = 0; i < datasetCount; i++)
                    {
                        if (features[i] == null)
                            dataLine.Add("0");
                        else
                            dataLine.Add(PRISM.StringUtilities.DblToString(features[i].Abundance, 2));
                    }

                    for (var i = 0; i < datasetCount; i++)
                    {
                        if (features[i] == null)
                            dataLine.Add("0");
                        else
                        {
                            if (features[i].Score <= float.MinValue)
                                dataLine.Add(PRISM.StringUtilities.DblToStringScientific(float.MinValue, 2));
                            else
                                dataLine.Add(PRISM.StringUtilities.DblToString(features[i].Score, 3));
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
                    for (var i = 0; i < datasetCount; i++)
                    {
                        if (features[i] == null)
                            dataLine.Add("0");
                        else
                            dataLine.Add(features[i].ProteinSpectrumMatches.Count.ToString());
                    }

                    writer.WriteLine(string.Join("\t", dataLine));
                }

            }

            Console.WriteLine("Results written to " + outputFilePath);
        }
    }
}
