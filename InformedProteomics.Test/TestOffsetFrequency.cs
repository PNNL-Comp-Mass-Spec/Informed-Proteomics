using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.LikelihoodScoring;
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
        private int _searchWidth;

        private ActivationMethod _act;
        private int _precursorCharge;
        private double _binWidth;

        [Test]
        public void OffsetFrequencyFunction()
        {
            InitTest(new ConfigFileReader(@"C:\Users\wilk011\Documents\DataFiles\OffsetFreqConfig.ini"));

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

                var offsetFrequencyTables = new OffsetFrequencyTable[tableCount];
                for (int i = 0; i < tableCount; i++)
                {
                    if (_precursorCharge > 0)
                        offsetFrequencyTables[i] = new OffsetFrequencyTable(_searchWidth, i+1, _binWidth/(i+1));
                    else
                        offsetFrequencyTables[i] = new OffsetFrequencyTable(_searchWidth, i+1, _binWidth);
                }


                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", Path.GetFileName(textFile), Path.GetFileName(rawFile));
                    var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, _noiseFiltration, _noiseFiltration);
                    var matchList = new SpectrumMatchList(lcms, new TsvFileParser(txtFiles[i]), _act);

                    for (int j = 0; j < tableCount; j++)
                    {
                        offsetFrequencyTables[j].AddMatches(
                            _precursorCharge > 0 ? matchList.GetCharge(j + 1) : matchList.Matches);
                    }
                }

                var offsetFrequencies = new List<OffsetProbability>[tableCount];
                for (int i = 0; i < tableCount; i++)
                {
                    offsetFrequencies[i] = offsetFrequencyTables[i].OffsetFrequencies;
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
                            outFile.WriteLine("{0}\t{1}", offsetFrequency.Offset, offsetFrequency.Prob);
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
            switch (actStr)
            {
                case "hcd":
                    _act = ActivationMethod.HCD;
                    break;
                case "cid":
                    _act = ActivationMethod.CID;
                    break;
                case "etd":
                    _act = ActivationMethod.ETD;
                    break;
            }

            _noiseFiltration = Convert.ToDouble(config.Contents["noisefiltration"]);
            _precursorCharge = Convert.ToInt32(config.Contents["precursorcharge"]);
            _searchWidth = Convert.ToInt32(config.Contents["searchwidth"]);
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
