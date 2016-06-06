using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.LikelihoodScoring;
using InformedProteomics.Scoring.LikelihoodScoring.Config;
using InformedProteomics.Scoring.LikelihoodScoring.Data;
using InformedProteomics.Scoring.LikelihoodScoring.FileReaders;
using InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestOffsetFrequency
    {
        private string[] _names;
        private string _preTsv;
        private string _preRaw;
        private string _outPre;
        private string _outFileName;
        private double _noiseFiltration;
        private double _searchWidth;

        private int _precursorCharge;
        private double _binWidth;

        [Test]
        public void OffsetFrequencyFunction()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string configFilePath = @"C:\Users\wilk011\Documents\DataFiles\OffsetFreqConfig.ini";
            if (!File.Exists(configFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, configFilePath);
            }

            InitTest(new ConfigFileReader(configFilePath));

            foreach (var name in _names)
            {
                var tsvName = _preTsv.Replace("@", name);
                var rawName = _preRaw.Replace("@", name);
                var txtFiles = Directory.GetFiles(tsvName).ToList();
                var rawFilesTemp = Directory.GetFiles(rawName).ToList();
                var rawFiles = rawFilesTemp.Where(rawFile => Path.GetExtension(rawFile) == ".raw").ToList();

                Assert.True(rawFiles.Count == txtFiles.Count);

                int tableCount = 1;
                if (_precursorCharge > 0)
                    tableCount = _precursorCharge;

                var offsetFrequencyTables = new PrecursorOffsetFrequencyTable[tableCount];
                for (int i = 0; i < tableCount; i++)
                {
                    if (_precursorCharge > 0)
                        offsetFrequencyTables[i] = new PrecursorOffsetFrequencyTable(_searchWidth/(i+1), i+1, _binWidth/(i+1));
                    else
                        offsetFrequencyTables[i] = new PrecursorOffsetFrequencyTable(_searchWidth, i+1, _binWidth);
                }

                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", Path.GetFileName(textFile), Path.GetFileName(rawFile));
                    var lcms = new LazyLcMsRun(rawFile, _noiseFiltration, _noiseFiltration);
                    var matchList = new SpectrumMatchList(lcms, txtFiles[i], DataFileFormat.IcBottomUp);

                    for (int j = 0; j < tableCount; j++)
                    {
                        offsetFrequencyTables[j].AddMatches(
                            _precursorCharge > 0 ? matchList.GetCharge(j + 1) : matchList);
                    }
                }

                var offsetFrequencies = new List<Probability<double>>[tableCount];
                for (int i = 0; i < tableCount; i++)
                {
                    offsetFrequencies[i] = offsetFrequencyTables[i].GetProbabilities().ToList();
                }

                var outFileName = _outFileName.Replace("@", name);
                for (int i = 0; i < tableCount; i++)
                {
                    var chargeOutFileName = outFileName.Replace("*", (i + 1).ToString(CultureInfo.InvariantCulture));
                    using (var outFile = new StreamWriter(chargeOutFileName))
                    {
                        outFile.WriteLine("Offset\tFound");
                        foreach (var offsetFrequency in offsetFrequencies[i])
                        {
//                            int integerOffset = (int) Math.Round(offsetFrequency.Offset*_binWidth*(i + 1));
                            outFile.WriteLine("{0}\t{1}", offsetFrequency.Label, offsetFrequency.Prob);
                        }
                    }
                }
            }
        }

        [Test]
        public void PrecursorOffsetFrequencyFunction()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string configFilePath = @"C:\Users\wilk011\Documents\DataFiles\PrecursorOffsetFreqConfig.ini";
            if (!File.Exists(configFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, configFilePath);
            }

            InitTest(new ConfigFileReader(configFilePath));

            foreach (var name in _names)
            {
                var tsvName = _preTsv.Replace("@", name);
                var rawName = _preRaw.Replace("@", name);
                var txtFiles = Directory.GetFiles(tsvName).ToList();
                var rawFilesTemp = Directory.GetFiles(rawName).ToList();
                var rawFiles = rawFilesTemp.Where(rawFile => Path.GetExtension(rawFile) == ".raw").ToList();

                Assert.True(rawFiles.Count == txtFiles.Count);

                var offsetFrequencyTables = new List<PrecursorOffsetFrequencyTable>[_precursorCharge];

                for (int i = 0; i < _precursorCharge; i++)
                {
                    offsetFrequencyTables[i] = new List<PrecursorOffsetFrequencyTable>();
                    for (int j = 1; j <= (i + 1); j++)
                    {
                        offsetFrequencyTables[i].Add(new PrecursorOffsetFrequencyTable(_searchWidth/j, j, _binWidth/j));
                    }
                }

                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", Path.GetFileName(textFile), Path.GetFileName(rawFile));
                    var lcms = new LazyLcMsRun(rawFile, _noiseFiltration, _noiseFiltration);
                    var matchList = new SpectrumMatchList(lcms, txtFiles[i], DataFileFormat.IcBottomUp);

                    for (int j = 0; j < _precursorCharge; j++)
                    {
                        var chargeMatches = matchList.GetCharge(j + 1);

                        foreach (var offsetFrequencyTable in offsetFrequencyTables[j])
                        {
                            offsetFrequencyTable.AddMatches(chargeMatches);
                        }
                    }
                }

                var outFileName = _outFileName.Replace("@", name);
                for (int i = 0; i < _precursorCharge; i++)
                {
                    var chargeOutFileName = outFileName.Replace("*", (i + 1).ToString(CultureInfo.InvariantCulture));
                    using (var outFile = new StreamWriter(chargeOutFileName))
                    {
                        outFile.Write("Offset\t");
                        var offsetFrequencyList = new List<List<Probability<double>>>();
                        for (int j = offsetFrequencyTables[i].Count-1; j >=0; j--)
                        {
                            offsetFrequencyList.Add(offsetFrequencyTables[i][j].GetProbabilities().ToList());
                            outFile.Write("{0}", (j+1));
                            if (j != 0)
                                outFile.Write("\t");
                        }
                        outFile.WriteLine();
                        for (int j = 0; j < offsetFrequencyList.First().Count; j++)
                        {
                            var offset = offsetFrequencyList.First()[j].Label;
                            var integerOffset = Math.Round(offset*(1/_binWidth), 2);
                            outFile.Write(integerOffset+"\t");
                            for (int k = 0; k < offsetFrequencyList.Count; k++)
                            {
                                outFile.Write(offsetFrequencyList[k][j].Prob);
                                if (k != offsetFrequencyList.Count-1)
                                    outFile.Write("\t");
                            }
                            outFile.WriteLine();
                        }
                    }
                }
            }
        }


        // Read Configuration file
        private void InitTest(ConfigFileReader reader)
        {
            // Read program variables
            var config = reader.GetNodes("vars").First();
            var actStr = config.Contents["activationmethod"].ToLower();

            _noiseFiltration = Convert.ToDouble(config.Contents["noisefiltration"]);
            _precursorCharge = Convert.ToInt32(config.Contents["precursorcharge"]);
            _searchWidth = Convert.ToDouble(config.Contents["searchwidth"]);
            _binWidth = Convert.ToDouble(config.Contents["binwidth"]);

            // Read input and output file names
            var fileInfo = reader.GetNodes("fileinfo").First();
            _names = fileInfo.Contents["name"].Split(',');
            _preTsv = fileInfo.Contents["tsvpath"];
            _preRaw = fileInfo.Contents["rawpath"];
            var outPathtemp = fileInfo.Contents["outpath"];
            _outPre = outPathtemp;
            var outFiletemp = fileInfo.Contents["outfile"];
            _outFileName = _outPre + outFiletemp;
        }
    }
}
