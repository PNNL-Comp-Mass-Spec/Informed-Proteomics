using System.IO;
using System;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.FeatureFinding.Alignment;
using InformedProteomics.FeatureFinding.SpectrumMatching;
using MathNet.Numerics.Statistics;
using PRISM;

namespace PromexAlign
{
    public class Program
    {
        public static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                ShowSyntax();
                return;
            }

            // Parse file
            var inputFilePath = args[0];

            if (!File.Exists(inputFilePath))
            {
                ConsoleMsgUtils.ShowError("File not found: " + inputFilePath);
                return;
            }

            var datasets = DatasetInfo.ParseDatasetInfoFile(inputFilePath);

            if (datasets.Count == 0)
            {
                ConsoleMsgUtils.ShowError("No valid data found in the dataset info file");
                ShowSyntax();
                return;
            }

            var fileName = Path.GetFileNameWithoutExtension(inputFilePath);
            var directory = Path.GetDirectoryName(inputFilePath);

            var crosstabFilename = string.Format("{0}_crosstab.tsv", fileName);

            string outputfilePath;
            if (string.IsNullOrWhiteSpace(directory))
                outputfilePath = crosstabFilename;
            else
                outputfilePath= Path.Combine(directory, crosstabFilename);

            var nDataset = datasets.Count;
            var prsmReader = new ProteinSpectrumMatchReader();
            var tolerance = new Tolerance(100);
            var alignment = new LcMsFeatureAlignment(new CompRefFeatureComparer(tolerance));

            var dataId = 0;
            foreach (var dataset in datasets)
            {
                if (!File.Exists(dataset.RawFilePath))
                {
                    ConsoleMsgUtils.ShowError("Instrument file not found: " + dataset.RawFilePath);
                    continue;
                }

                if (!File.Exists(dataset.Ms1FtFilePath))
                {
                    ConsoleMsgUtils.ShowError("ProMex results file not found: " + dataset.Ms1FtFilePath);
                    continue;
                }

                Console.WriteLine("Opening " + dataset.RawFilePath);
                var run = PbfLcMsRun.GetLcMsRun(dataset.RawFilePath, 0, 0);

                Console.WriteLine("Opening " + dataset.Ms1FtFilePath);
                var features = LcMsFeatureAlignment.LoadProMexResult(dataId++, dataset.Ms1FtFilePath, run);

                if (!string.IsNullOrWhiteSpace(dataset.MsPfIdFilePath) && File.Exists(dataset.MsPfIdFilePath))
                {
                    Console.WriteLine("Opening " + dataset.MsPfIdFilePath);
                    var prsmList = prsmReader.LoadIdentificationResult(dataset.MsPfIdFilePath, ProteinSpectrumMatch.SearchTool.MsPathFinder);

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

                alignment.AddDataSet(dataId, features, run);
            }

            alignment.AlignFeatures();

            Console.WriteLine("{0} alignments ", alignment.CountAlignedFeatures);

            var validResults = 0;
            for (var datasetIndex = 0; datasetIndex < nDataset; datasetIndex++)
            {
                if (datasetIndex >= alignment.CountDatasets)
                {
                    ConsoleMsgUtils.ShowWarning(string.Format("Could not align {0}; features not found", datasets[datasetIndex].Label));
                    continue;
                }

                alignment.FillMissingFeatures(datasetIndex);
                Console.WriteLine("{0} has been processed", datasets[datasetIndex].Label);
                validResults++;
            }

            if (validResults > 0)
            {
                OutputCrossTabWithId(outputfilePath, alignment, datasets.Select(ds => ds.Label).ToArray());
            }
        }

        public static void OutputCrossTabWithId(string outputFilePath, LcMsFeatureAlignment alignment, string[] runLabels)
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

            foreach (var dataName in runLabels)
            {
                writer.Write("\t");
                writer.Write(dataName + "_FeatureId");
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
            foreach (var dataName in runLabels)
            {
                writer.Write("\t");
                writer.Write(dataName + "_SpectraCount");
            }
            writer.Write("\n");

            var alignedFeatureList = alignment.GetAlignedFeatures();
            foreach (var features in alignedFeatureList)
            {
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

                for (var i = 0; i < nDataset; i++)
                {
                    writer.Write("\t");
                    writer.Write(features[i] == null ? 0 : features[i].FeatureId);
                }

                var prsm = (from f in features
                            where f?.ProteinSpectrumMatches != null && f.ProteinSpectrumMatches.Count > 0
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
                for (var i = 0; i < nDataset; i++)
                {
                    writer.Write("\t");
                    writer.Write(features[i] == null ? 0 : features[i].ProteinSpectrumMatches.Count);
                }
                writer.Write("\n");
            }

            writer.Close();
        }

        private static void ShowSyntax()
        {
            Console.WriteLine();
            Console.WriteLine("Syntax:");
            Console.WriteLine("PromexAlign.exe DatasetInfoFile");
            Console.WriteLine();
            Console.WriteLine("The dataset info file is a tab delimited text file with the following columns of info");
            Console.WriteLine("Label  RawFilePath  Ms1FtFilePath  MsPathfinderIdFilePath");
            Console.WriteLine();
            Console.WriteLine("The MsPathfinderIdFilePath column is optional");
        }

    }
}
