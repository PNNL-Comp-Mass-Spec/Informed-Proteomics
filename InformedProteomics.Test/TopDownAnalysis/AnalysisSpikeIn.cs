using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Execution;
using MathNet.Numerics.Statistics;
using NUnit.Framework;

namespace InformedProteomics.Test.TopDownAnalysis
{
    public class AnalysisSpikeIn
    {
        private string[] datasets = new string[]
        {
            "CPTAC_Intact_Spike_1x_1_27Apr15_Bane_14-09-03RZ",
            "CPTAC_Intact_Spike_1x_2_27Apr15_Bane_14-09-03RZ",
            "CPTAC_Intact_Spike_1x_3_27Apr15_Bane_14-09-03RZ",
            "CPTAC_Intact_Spike_1x_4_27Apr15_Bane_14-09-03RZ",
            "CPTAC_Intact_Spike_1x_5_27Apr15_Bane_14-09-03RZ",
            "CPTAC_Intact_Spike_5x_1_27Apr15_Bane_14-09-03RZ",
            "CPTAC_Intact_Spike_5x_2_27Apr15_Bane_14-09-03RZ",
            "CPTAC_Intact_Spike_5x_3b_30Apr15_Bane_14-09-01RZ",
            "CPTAC_Intact_Spike_5x_4_27Apr15_Bane_14-09-01RZ",
            "CPTAC_Intact_Spike_5x_5_27Apr15_Bane_14-09-01RZ",
            "CPTAC_Intact_Spike_10x_1_27Apr15_Bane_14-09-01RZ",
            "CPTAC_Intact_Spike_10x_2_27Apr15_Bane_14-09-01RZ",
            "CPTAC_Intact_Spike_10x_3_2May15_Bane_14-09-01RZ",
            "CPTAC_Intact_Spike_10x_4_27Apr15_Bane_14-09-01RZ",
            "CPTAC_Intact_Spike_10x_5_27Apr15_Bane_14-09-01RZ"
        };

        public const string Ms1FtFolder = @"\\protoapps\UserData\Jungkap\Quant";
        public const string MsPfFolder = @"\\protoapps\UserData\Jungkap\Quant\MSPF";
        public const string RawFolder = @"\\protoapps\UserData\Jungkap\Quant";

        internal class SpikeInFeatureComparer : INodeComparer<LcMsFeature>
        {
            public SpikeInFeatureComparer(Tolerance tolerance = null)
            {
                _tolerance = tolerance ?? new Tolerance(10);
            }

            public bool SameCluster(LcMsFeature f1, LcMsFeature f2)
            {
                if (f1.DataSetId == f2.DataSetId) return false;
                // tolerant in mass dimension?

                var massTol = Math.Min(_tolerance.GetToleranceAsTh(f1.Mass), _tolerance.GetToleranceAsTh(f2.Mass));
                if (Math.Abs(f1.Mass - f2.Mass) > massTol) return false;

                //if (!f1.CoElutedByNet(f2, 0.004)) return false; //e.g) 200*0.001 = 0.2 min = 30 sec
                if (f1.ProteinSpectrumMatches.Count > 0 && f2.ProteinSpectrumMatches.Count > 0)
                {
                    return f1.CoElutedByNet(f2, 0.01) && f1.ProteinSpectrumMatches.ShareProteinId(f2.ProteinSpectrumMatches);
                }

                return f1.CoElutedByNet(f2, 0.01);
            }

            private readonly Tolerance _tolerance;
        }

        [Test]
        public void TestFeatureAlignment()
        {
            const string outFilePath = @"\\protoapps\UserData\Jungkap\Quant\aligned\promex_crosstab.tsv";
            //const string outFolder = @"\\protoapps\UserData\Jungkap\CompRef\aligned";
            var runLabels = new string[] { "1x1", "1x2", "1x3", "1x4", "1x5", "5x1", "5x2", "5x3", "5x4", "5x5", "10x1", "10x2", "10x3", "10x4", "10x5", };
            var nDataset = runLabels.Length;

            var prsmReader = new ProteinSpectrumMatchReader();
            var tolerance = new Tolerance(10);
            var alignment = new LcMsFeatureAlignment(new SpikeInFeatureComparer(tolerance));

            for (var i = 0; i < nDataset; i++)
            {
                var rawFile = string.Format(@"{0}\{1}.pbf", RawFolder, datasets[i]);
                var mspFile = string.Format(@"{0}\{1}_IcTda.tsv", MsPfFolder, datasets[i]);
                var ms1FtFile = string.Format(@"{0}\{1}.ms1ft", Ms1FtFolder, datasets[i]);

                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                var prsmList = prsmReader.LoadIdentificationResult(mspFile, ProteinSpectrumMatch.SearchTool.MsPathFinder);
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
            /*
            for (var i = 0; i < nDataset; i++)
            {
                alignment.FillMissingFeatures(i);
                Console.WriteLine("{0} has been processed", runLabels[i]);
            }
            */
            OutputCrossTabWithId(outFilePath, alignment, runLabels);
        }

        public void OutputCrossTabWithId(string outputFilePath, LcMsFeatureAlignment alignment, string[] runLabels)
        {
            var nDataset = runLabels.Length;
            var writer = new StreamWriter(outputFilePath);

            writer.Write("MonoMass");
            writer.Write("\t");
            writer.Write("MinElutionTime");
            writer.Write("\t");
            writer.Write("MaxElutionTime");

            foreach (var dataName in runLabels)
            {
                writer.Write("\t");
                writer.Write(dataName + "_Abundance");
            }

            foreach (var dataName in runLabels)
            {
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
            writer.Write("SequenceText");
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
            foreach (var dataName in runLabels)
            {
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

                for (var i = 0; i < nDataset; i++)
                {
                    writer.Write("\t");
                    writer.Write(features[i] == null ? 0 : features[i].Abundance);
                }

                for (var i = 0; i < nDataset; i++)
                {
                    writer.Write("\t");
                    writer.Write(features[i] == null ? 0 : features[i].Score);
                }

                var prsm = (from f in features
                            where f != null && f.ProteinSpectrumMatches != null && f.ProteinSpectrumMatches.Count > 0
                            select f.ProteinSpectrumMatches[0]).FirstOrDefault();

                if (prsm == null)
                {
                    for (var k = 0; k < 10; k++)
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
                    writer.Write(prsm.SequenceText);
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
                for (var i = 0; i < nDataset; i++)
                {
                    writer.Write("\t");
                    writer.Write(features[i] == null ? 0 : features[i].ProteinSpectrumMatches.Count);
                }
                writer.Write("\n");
            }
            writer.Close();
        }
    }
}
