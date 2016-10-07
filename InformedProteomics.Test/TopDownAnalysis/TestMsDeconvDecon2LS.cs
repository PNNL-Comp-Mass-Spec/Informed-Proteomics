using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.FeatureFinding.MassFeature;
using InformedProteomics.TopDown.Execution;
using NUnit.Framework;
namespace InformedProteomics.Test.TopDownAnalysis
{
    public class TestMsDeconvDecon2LS
    {
        private string[] spikeDatasets = new string[]
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

        [Test]
        public void TestCptacSpikeIn()
        {
            const string featureFolder = @"D:\MassSpecFiles\CPTAC_spike_in\promex";
            const string rawFolder = @"D:\MassSpecFiles\CPTAC_spike_in\raw";
            var outFilePath = string.Format(@"{0}\aligned_features.tsv", featureFolder);
            var align = new LcMsFeatureAlignment(new LcMsFeatureAlignComparer(new Tolerance(10)));

            for (var i = 0; i < spikeDatasets.Length; i++)
            {
                var featureFilePath = string.Format(@"{0}\{1}.ms1ft", featureFolder, spikeDatasets[i]);
                var rawFile = string.Format(@"{0}\{1}.pbf", rawFolder, spikeDatasets[i]);

                if (!File.Exists(rawFile))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", rawFile);
                    continue;
                }

                if (!File.Exists(featureFilePath))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", featureFilePath);
                    continue;
                }
                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                var s = 0d;
                foreach (var scanNum in run.GetMs1ScanVector())
                {
                    var spec = run.GetSpectrum(scanNum);
                    var summedIntensity = spec.Peaks.Sum(p => p.Intensity);
                    s += summedIntensity;
                }
                foreach (var scanNum in run.GetScanNumbers(2))
                {
                    var spec = run.GetSpectrum(scanNum);
                    var summedIntensity = spec.Peaks.Sum(p => p.Intensity);
                    s += summedIntensity;
                }

                Console.WriteLine("{0}\t{1}", i, s);
                //var features = LcMsFeatureAlignment.LoadProMexResult(i, featureFilePath, run);

                //align.AddDataSet(i, features, run);
            }
            //align.AlignFeatures();
            //Console.WriteLine("# of aligned features = {0}", align.CountAlignedFeatures);
            //align.RefineAbundance();
            //OutputAlignmentResult(align, outFilePath, spikeDatasets);
        }

        [Test]
        public void TestCptac10Replicates()
        {
            const string featureFolder = @"D:\MassSpecFiles\CPTAC_rep10\icr2ls";
            const string rawFolder = @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1";
            var outFilePath = string.Format(@"{0}\aligned_features.tsv", featureFolder);
            var align = new LcMsFeatureAlignment(new LcMsFeatureAlignComparer(new Tolerance(10)));

            var dataNames = new string[10];
            for (var i = 0; i < 10; i++)
            {
                dataNames[i] = string.Format(@"CPTAC_Intact_rep{0}_15Jan15_Bane_C2-14-08-02RZ", i+1);
                var featureFilePath = string.Format(@"{0}\{1}_isos.tsv", featureFolder, dataNames[i]);
                var rawFile = string.Format(@"{0}\{1}.pbf", rawFolder, dataNames[i]);

                if (!File.Exists(rawFile))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", rawFile);
                    continue;
                }

                if (!File.Exists(featureFilePath))
                {
                    Console.WriteLine(@"Warning: Skipping file not found: {0}", featureFilePath);
                    continue;
                }
                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                var features = LcMsFeatureAlignment.LoadProMexResult(i, featureFilePath, run);

                align.AddDataSet(i, features, run);
            }
            align.AlignFeatures();
            Console.WriteLine("# of aligned features = {0}", align.CountAlignedFeatures);
            //var tempOutPath = outFilePath + ".tmp";
            //OutputAlignmentResult(align, tempOutPath, rawFiles, true);
            //align.RefineAbundance();
            OutputAlignmentResult(align, outFilePath, dataNames);
        }

        private void OutputAlignmentResult(LcMsFeatureAlignment align, string outFilePath, string[] dataName)
        {
            var alignedFeatureList = align.GetAlignedFeatures();

            var writer = new StreamWriter(outFilePath);
            writer.Write("MonoMass\tMinElutionTime\tMaxElutionTime");
            for (var i = 0; i < align.CountDatasets; i++)
            {
                writer.Write("\t{0}", dataName[i]);
            }
            writer.Write("\n");

            for (var i = 0; i < align.CountAlignedFeatures; i++)
            {
                var features = alignedFeatureList[i];
                var minMaxNet = TestLcMsFeatureAlignment.GetMinMaxNet(features);
                writer.Write(@"{0}\t{1:0.00000}\t{2:0.00000}", minMaxNet.Item1, minMaxNet.Item3, minMaxNet.Item4);

                for (var j = 0; j < align.CountDatasets; j++)
                {
                    var feature = features[j];
                    writer.Write("\t");
                    writer.Write(feature != null ? feature.Abundance : 0d);
                }
                writer.Write("\n");
            }
            writer.Close();
        }
    }
}
