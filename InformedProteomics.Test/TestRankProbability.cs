using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Scoring.LikelihoodScoring;
using InformedProteomics.Scoring.LikelihoodScoring.Config;
using InformedProteomics.Scoring.LikelihoodScoring.Data;
using InformedProteomics.Scoring.LikelihoodScoring.FileReaders;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestRankProbability
    {
        [Test]
        public void RankProbability()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string configurationFile = @"\\protoapps\UserData\Wilkins\BottomUp\RankProbabilityConfig.ini";
            if (!File.Exists(configurationFile))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, configurationFile);
            }

            var config = new TrainerConfiguration(configurationFile);

            var trainingParameters = new TrainingParameters(config);

            foreach (var dataSet in config.DataSets)
            {
                var tsvName = config.TsvPath.Replace("@", dataSet);
                var rawName = config.DataPath.Replace("@", dataSet);
                var txtFiles = new List<string>();
                var dataFiles = new List<string>();

                if (config.DataFormat == DataFileFormat.Mgf)
                {
                    dataFiles.Add(config.DataPath + dataSet + "." + config.DataFormat);
                }
                else
                {
                    // Read directory
                    var txtFilesTemp = Directory.GetFiles(tsvName).ToList();
                    txtFiles = txtFilesTemp.Where(txtFile => Path.GetExtension(txtFile) == ".tsv").ToList();
                    var dataFilesTemp = Directory.GetFiles(rawName).ToList();
                    dataFiles =
                        dataFilesTemp.Where(dataFile => Path.GetExtension(dataFile) == ".raw").ToList();
                    Assert.True(dataFiles.Count == txtFiles.Count);
                }

                // Read Data files
                for (int i = 0; i < dataFiles.Count; i++)
                {
                    SpectrumMatchList targets;

                    if (config.DataFormat == DataFileFormat.Mgf)
                    {
                        Console.WriteLine(dataSet + ".mgf");
                        targets = new SpectrumMatchList(dataFiles[i], config.PrecursorCharge);
                    }
                    else
                    {
                        string textFile = txtFiles[i];
                        string rawFile = dataFiles[i];
                        Console.WriteLine("{0}\t{1}", textFile, rawFile);
                        var lcms = new LazyLcMsRun(rawFile, 0, 0);
                        targets = new SpectrumMatchList(lcms, txtFiles[i], config.DataFormat, config.PrecursorCharge);
                    }

                    // Add PSMs
                    trainingParameters.AddMatches(targets);
                }

                // Write output file
                var outputFileName = config.OutputFileName.Replace("@", dataSet);
                trainingParameters.WriteToFile(outputFileName);
            }
        }
    }
}
