using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.LikelihoodScoring;
using InformedProteomics.Scoring.LikelihoodScoring.Config;
using InformedProteomics.Scoring.LikelihoodScoring.Data;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestRankProbability
    {
        [Test]
        public void RankProbability()
        {
            const string configurationFile = @"\\protoapps\UserData\Wilkins\BottomUp\RankProbabilityConfig.ini";
            var config = new TrainerConfiguration(configurationFile);

            var trainingParameters = new TrainingParameters(config);

            foreach (var dataSet in config.DataSets)
            {
                var tsvName = config.TsvPath.Replace("@", dataSet);
                var rawName = config.DataPath.Replace("@", dataSet);
                var txtFiles = new List<string>();
                var dataFiles = new List<string>();

                if (config.DataFormat == DataFileFormat.Msgf)
                {
                    // Read directory
                    var txtFilesTemp = Directory.GetFiles(tsvName).ToList();
                    txtFiles = txtFilesTemp.Where(txtFile => Path.GetExtension(txtFile) == ".tsv").ToList();
                    var dataFilesTemp = Directory.GetFiles(rawName).ToList();
                    dataFiles =
                        dataFilesTemp.Where(dataFile => Path.GetExtension(dataFile) == ".raw").ToList();
                    Assert.True(dataFiles.Count == txtFiles.Count);
                }
                else
                {
                    dataFiles.Add(config.DataPath + dataSet + "." + config.DataFormat);
                }

                // Read Data files
                for (int i = 0; i < dataFiles.Count; i++)
                {
                    var targets = new SpectrumMatchList(config.ActivationMethod, false, config.PrecursorCharge,
                        config.DataFormat);
                    var decoys = new SpectrumMatchList(config.ActivationMethod, true, config.PrecursorCharge,
                        config.DataFormat);

                    if (config.DataFormat == DataFileFormat.Msgf)
                    {
                        string textFile = txtFiles[i];
                        string rawFile = dataFiles[i];
                        Console.WriteLine("{0}\t{1}", textFile, rawFile);
                        var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0);
                        targets.AddMatchesFromTsvFile(lcms, new TsvFileParser(txtFiles[i]));
                    }
                    else
                    {
                        Console.WriteLine(dataSet + ".mgf");
                        targets.AddMatchesFromMgfFile(dataFiles[i]);
                    }

                    decoys.AddMatches(targets);

                    // Add PSMs
                    trainingParameters.AddMatches(targets, decoys);
                }

                // Write output file
                var outputFileName = config.OutputFileName.Replace("@", dataSet);
                trainingParameters.WriteToFile(outputFileName);
            }
        }
    }
}
