using System.IO;
using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.FeatureFinding.Alignment;
using InformedProteomics.FeatureFinding.Data;
using InformedProteomics.FeatureFinding.SpectrumMatching;
using MathNet.Numerics.Statistics;
using PRISM;
using ProMexAlign;

namespace PromexAlign
{
    public static class Program
    {
        // Ignore Spelling: pre

        public static void Main(string[] args)
        {
            try
            {
                if (args.Length == 0)
                {
                    ShowSyntax();
                    ConsoleMsgUtils.SleepSeconds(0.5);
                    return;
                }

                // Parse file
                var inputFilePath = args[0];

                var datasetInfoFile = new FileInfo(inputFilePath);
                if (!datasetInfoFile.Exists)
                {
                    ConsoleMsgUtils.ShowError("Dataset info file not found: " + inputFilePath);
                    if (!datasetInfoFile.FullName.Equals(inputFilePath, StringComparison.OrdinalIgnoreCase))
                    {
                        ConsoleMsgUtils.ShowWarning("Full path to the file: " + datasetInfoFile.FullName);
                    }

                    ConsoleMsgUtils.SleepSeconds(1.5);
                    return;
                }

                var datasets = DatasetInfo.ParseDatasetInfoFile(inputFilePath);

                if (datasets.Count == 0)
                {
                    ConsoleMsgUtils.ShowError("No valid data found in the dataset info file");
                    ConsoleMsgUtils.ShowWarning("See " + inputFilePath);
                    ShowSyntax();

                    ConsoleMsgUtils.SleepSeconds(1.5);
                    return;
                }

                AlignDatasets(datasetInfoFile, datasets);

                Console.WriteLine();
                Console.WriteLine("Processing Complete");
            }
            catch (Exception ex)
            {
                ConsoleMsgUtils.ShowError("Error in entry method", ex);
            }
        }

        private static bool AddPrsms(ProteinSpectrumMatchReader prsmReader, FileSystemInfo prsmFile, IEnumerable<LcMsFeature> features, Tolerance tolerance)
        {
            try
            {
                Console.WriteLine("Opening {0}", PathUtils.CompactPathString(prsmFile.FullName, 80));
                var prsmList = prsmReader.LoadIdentificationResult(prsmFile.FullName, ProteinSpectrumMatch.SearchTool.MsPathFinder);

                if (prsmList.Count == 0)
                    return false;

                foreach (var match in prsmList)
                {
                    match.ProteinId = match.ProteinName;
                }

                // Tag features by PrSMs
                foreach (var feature in features)
                {
                    var massTol = tolerance.GetToleranceAsMz(feature.Mass);
                    foreach (var match in prsmList)
                    {
                        if (feature.MinScanNum < match.ScanNum && match.ScanNum < feature.MaxScanNum && Math.Abs(feature.Mass - match.Mass) < massTol)
                        {
                            feature.ProteinSpectrumMatches.Add(match);
                        }
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                if (ex.Message.StartsWith("You can ignore bad data", StringComparison.OrdinalIgnoreCase))
                {
                    // The error message comes from CsvHelper; it won't be helpful to the end user so don't show it
                    ConsoleMsgUtils.ShowError("Error loading MSPathFinderResults file {0}; data format error", prsmFile.Name);
                    ConsoleMsgUtils.ShowWarning(StackTraceFormatter.GetExceptionStackTraceMultiLine(ex));
                }
                else
                {
                    ConsoleMsgUtils.ShowError(string.Format("Error loading MSPathFinderResults file {0}", prsmFile.Name), ex);
                }

                return false;
            }
        }

        private static void AlignDatasets(FileInfo datasetInfoFile, IReadOnlyList<DatasetInfo> datasets)
        {
            try
            {
                var baseName = Path.GetFileNameWithoutExtension(datasetInfoFile.Name);

                var workingDirectoryPath = datasetInfoFile.Directory == null ? "." : datasetInfoFile.Directory.FullName;

                var outputFilePath = Path.Combine(workingDirectoryPath, string.Format("{0}_crosstab.tsv", baseName));

                var datasetsWithoutSpectraData = new List<DatasetInfo>();
                var featureCountByDataset = new Dictionary<DatasetInfo, int>();

                var datasetCount = datasets.Count;

                var prsmReader = new ProteinSpectrumMatchReader();
                var tolerance = new Tolerance(100);
                var alignment = new LcMsFeatureAlignment(new CompRefFeatureComparer(tolerance));

                var dataId = 0;
                var prsmsDefined = false;

                foreach (var dataset in datasets)
                {
                    Console.WriteLine();

                    var datasetFile = FindInputFile(workingDirectoryPath, dataset.RawFilePath, "Instrument file", true);

                    var featuresFile = FindInputFile(workingDirectoryPath, dataset.Ms1FtFilePath, "ProMex results file");
                    if (featuresFile == null)
                        continue;

                    var prsmFile = FindInputFile(workingDirectoryPath, dataset.MsPfIdFilePath, "MSPathFinder results file", true);

                    LcMsRun run;
                    if (datasetFile == null)
                    {
                        ConsoleMsgUtils.ShowWarning("Instrument file not found: {0}\nWill align ProMex features but cannot add missing features found in other datasets", dataset.RawFilePath);

                        run = new VirtualLcMsRun(Path.GetFileName(dataset.RawFilePath));
                        datasetsWithoutSpectraData.Add(dataset);
                    }
                    else
                    {
                        Console.WriteLine("Opening {0}", PathUtils.CompactPathString(datasetFile.FullName, 80));
                        run = PbfLcMsRun.GetLcMsRun(datasetFile.FullName, 0, 0);
                    }

                    Console.WriteLine("Opening {0}", PathUtils.CompactPathString(featuresFile.FullName, 80));
                    var features = LcMsFeatureAlignment.LoadProMexResult(dataId++, featuresFile.FullName, run);

                    if (prsmFile != null)
                    {
                        if (prsmFile.Extension.Equals(".mzid", StringComparison.OrdinalIgnoreCase))
                        {
                            ConsoleMsgUtils.ShowWarning(
                                "The MSPathFinder results file should be a .tsv file with columns Scan, Pre, Sequence, etc. (typically the _IcTda.tsv file);\n" +
                                "{0} files are not supported:\nignoring {1}", prsmFile.Extension, prsmFile.Name);
                        }
                        else
                        {
                            if (!prsmFile.Extension.Equals(".tsv", StringComparison.OrdinalIgnoreCase))
                            {
                                ConsoleMsgUtils.ShowWarning(string.Format(
                                    "The MSPathFinder results file should be a .tsv file with columns Scan, Pre, Sequence, etc. (typically the _IcTda.tsv file);\n" +
                                    "suffix {0} indicates this may not be a compatible file:\n{1}", prsmFile.Extension, prsmFile.Name));
                            }

                            var prsmsFound = AddPrsms(prsmReader, prsmFile, features, tolerance);
                            if (prsmsFound)
                            {
                                prsmsDefined = true;
                            }
                        }
                    }

                    alignment.AddDataSet(dataId, features, run);

                    featureCountByDataset.Add(dataset, features.Count);
                }

                alignment.AlignFeatures();

                Console.WriteLine();
                Console.WriteLine("{0} alignments ", alignment.CountAlignedFeatures);

                var validResults = 0;
                for (var datasetIndex = 0; datasetIndex < datasetCount; datasetIndex++)
                {
                    var dataset = datasets[datasetIndex];

                    if (datasetIndex >= alignment.CountDatasets)
                    {
                        ConsoleMsgUtils.ShowWarning("Could not align {0}; features not found", dataset.Label);
                        continue;
                    }

                    if (datasetsWithoutSpectraData.Contains(dataset))
                    {
                        var featureCount = featureCountByDataset[dataset];
                        if (featureCount > 0)
                            validResults++;

                        continue;
                    }

                    Console.WriteLine();
                    Console.WriteLine("Processing {0}", dataset.Label);

                    alignment.FillMissingFeatures(datasetIndex);
                    validResults++;
                }

                if (validResults == 0)
                {
                    Console.WriteLine("None of the datasets had valid results; aborting");
                    return;
                }

                OutputCrossTabWithId(outputFilePath, alignment, datasets.Select(dataset => dataset.Label).ToList(), prsmsDefined);
            }
            catch (Exception ex)
            {
                ConsoleMsgUtils.ShowError("Error in AlignDatasets", ex);
            }
        }

        private static FileInfo FindInputFile(string workingDirectoryPath, string inputFilePath, string inputFileDescription, bool optionalFile = false)
        {
            if (string.IsNullOrWhiteSpace(inputFilePath))
            {
                if (!optionalFile)
                    ConsoleMsgUtils.ShowError("{0} not defined", inputFileDescription);

                return null;
            }

            if (File.Exists(inputFilePath))
            {
                return new FileInfo(inputFilePath);
            }

            // File not found in the executable's working directory; examine workingDirectoryPath
            var alternatePath = Path.Combine(workingDirectoryPath, Path.GetFileName(inputFilePath));
            if (File.Exists(alternatePath))
            {
                return new FileInfo(alternatePath);
            }

            if (!optionalFile)
                ConsoleMsgUtils.ShowError("{0} not found: {1}", inputFileDescription, inputFilePath);

            return null;
        }

        public static void OutputCrossTabWithId(string outputFilePath, LcMsFeatureAlignment alignment, List<string> runLabels, bool prsmsDefined)
        {
            Console.WriteLine();
            Console.WriteLine("Writing aligned data to {0}", PathUtils.CompactPathString(outputFilePath, 80));

            using var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

            var headerNames = new List<string>
            {
                "MonoMass",
                "MinElutionTime",
                "MaxElutionTime"
            };

            headerNames.AddRange(runLabels.Select(dataName => dataName + "_Abundance"));

            headerNames.AddRange(runLabels.Select(dataName => dataName + "_Ms1Score"));

            headerNames.AddRange(runLabels.Select(dataName => dataName + "_FeatureId"));

            if (prsmsDefined)
            {
                headerNames.Add("Pre");
                headerNames.Add("Sequence");
                headerNames.Add("Post");
                headerNames.Add("Modifications");
                headerNames.Add("ProteinName");
                headerNames.Add("ProteinDesc");
                headerNames.Add("ProteinLength");
                headerNames.Add("Start");
                headerNames.Add("End");
                headerNames.Add("BestEValue");

                headerNames.AddRange(runLabels.Select(dataName => dataName + "_SpectraCount"));
            }

            writer.WriteLine(string.Join("\t", headerNames));

            var datasetCount = runLabels.Count;
            var dataColumns = new List<string>();

            foreach (var features in alignment.GetAlignedFeatures())
            {
                dataColumns.Clear();

                var mass = features.Where(f => f != null).Select(f => f.Mass).Median();
                var minElutionTime = features.Where(f => f != null).Select(f => f.MinElutionTime).Median();
                var maxElutionTime = features.Where(f => f != null).Select(f => f.MaxElutionTime).Median();

                dataColumns.Add(StringUtilities.DblToString(mass, 5));
                dataColumns.Add(StringUtilities.DblToString(minElutionTime, 3));
                dataColumns.Add(StringUtilities.DblToString(maxElutionTime, 3));

                for (var i = 0; i < datasetCount; i++)
                {
                    dataColumns.Add(features[i] == null ? "0" : StringUtilities.ValueToString(features[i].Abundance, 7));
                }

                for (var i = 0; i < datasetCount; i++)
                {
                    dataColumns.Add(features[i] == null ? "0" : StringUtilities.ValueToString(features[i].Score, 7));
                }

                for (var i = 0; i < datasetCount; i++)
                {
                    dataColumns.Add(features[i] == null ? "0" : features[i].FeatureId.ToString());
                }

                if (!prsmsDefined)
                {
                    writer.WriteLine(string.Join("\t", dataColumns));
                    continue;
                }

                var matchingPrsms = (from f in features
                            where f?.ProteinSpectrumMatches != null && f.ProteinSpectrumMatches.Count > 0
                            select f.ProteinSpectrumMatches[0]).ToList();

                if (matchingPrsms.Count == 0)
                {
                    for (var k = 0; k < 10; k++)
                    {
                        dataColumns.Add(string.Empty);
                    }
                }
                else
                {
                    var firstPrsm = matchingPrsms[0];

                    dataColumns.Add(firstPrsm.Pre);
                    dataColumns.Add(firstPrsm.Sequence);
                    dataColumns.Add(firstPrsm.Post);
                    dataColumns.Add(firstPrsm.Modifications);
                    dataColumns.Add(firstPrsm.ProteinName);
                    dataColumns.Add(firstPrsm.ProteinDesc);
                    dataColumns.Add(firstPrsm.ProteinLength.ToString());
                    dataColumns.Add(firstPrsm.FirstResidue.ToString());
                    dataColumns.Add(firstPrsm.LastResidue.ToString());

                    var bestEvalue = matchingPrsms.Select(prsm => prsm.SpectralEValue).Min();

                    dataColumns.Add(StringUtilities.ValueToString(bestEvalue, 5));
                }

                // spectral count from ms2
                for (var i = 0; i < datasetCount; i++)
                {
                    dataColumns.Add(string.Format("{0}", features[i]?.ProteinSpectrumMatches.Count ?? 0));
                }

                writer.WriteLine(string.Join("\t", dataColumns));
            }
        }

        private static void ShowSyntax()
        {
            Console.WriteLine();
            Console.WriteLine("Syntax:");
            Console.WriteLine("ProMexAlign.exe DatasetInfoFile");
            Console.WriteLine();
            Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                "The dataset info file is a tab-delimited text file with either 3 or 4 columns of information."));

            Console.WriteLine();
            Console.WriteLine("Expected columns:");
            Console.WriteLine("Label  RawFilePath  Ms1FtFilePath  MsPathfinderIdFilePath");
            Console.WriteLine();
            Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                "If data in the Label column is an empty string, dataset labels will be auto-assigned as Dataset_1, Dataset_2, etc."));
            Console.WriteLine();
            Console.WriteLine("The raw files are either Thermo .raw files or .pbf files created by PbfGen");
            Console.WriteLine();
            Console.WriteLine("The MsPathfinderIdFilePath column is optional");
        }
    }
}
