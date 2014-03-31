using System;
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
        private double _precursorCharge;
        private int _binWidth;

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

                var offsetFrequencyTable = new OffsetFrequencyTable(_searchWidth);

                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", Path.GetFileName(textFile), Path.GetFileName(rawFile));
                    var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, _noiseFiltration, _noiseFiltration);
                    var matchList = new SpectrumMatchList(lcms, new TsvFileParser(txtFiles[i]), _act);

                    offsetFrequencyTable.AddMatches(matchList.GetCharge(1), new Tolerance(0.5, ToleranceUnit.Da));
                }

                var offsetFrequencies = offsetFrequencyTable.OffsetFrequencies;

                var outFileName = _outFileName.Replace("@", name);
                using (var outFile = new StreamWriter(outFileName))
                {
                    outFile.WriteLine("Offset\tFound");
                    foreach (var offsetFrequency in offsetFrequencies)
                    {
                        outFile.WriteLine("{0}\t{1}", offsetFrequency.Offset, offsetFrequency.Found);
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
            _precursorCharge = Convert.ToDouble(config.Contents["precursorcharge"]);
            _searchWidth = Convert.ToInt32(config.Contents["searchwidth"]);
            _binWidth = Convert.ToInt32(config.Contents["binwidth"]);

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
