using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;
using InformedProteomics.Scoring.LikelihoodScoring;

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

        private List<IonType> _ionTypes;
        private IonTypeFactory _ionTypeFactory;
        private ActivationMethod _act;
        double _relativeIntensityThreshold = 1.0;
//        private bool _precursorOff;
        private bool _combineCharges;
        private int _precursorCharge;

        private readonly Tolerance _defaultTolerance = new Tolerance(15, ToleranceUnit.Ppm);

        [Test]
        public void OffsetFrequencyFunction()
        {
            InitTest(new ConfigFileReader(@"C:\Users\wilk011\Documents\DataFiles\OffsetFreqConfig_Test.ini"));

            foreach (var name in _names)
            {
                var tsvName = _preTsv.Replace("@", name);
                var rawName = _preRaw.Replace("@", name);
                var txtFiles = Directory.GetFiles(tsvName).ToList();
                var rawFilesTemp = Directory.GetFiles(rawName).ToList();
                var rawFiles = rawFilesTemp.Where(rawFile => Path.GetExtension(rawFile) == ".raw").ToList();

                Assert.True(rawFiles.Count == txtFiles.Count);

                var offsetFrequencyFunctions = new OffsetFrequencyTable[_precursorCharge];
                for (int i = 0; i < _precursorCharge; i++)
                {
                    offsetFrequencyFunctions[i] = new OffsetFrequencyTable();
                }

                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", textFile, rawFile);
                    var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0);
                    var matchList = new SpectrumMatchList(lcms, new TsvFileParser(txtFiles[i]), _act);
                    for (int j = 1; j <= _precursorCharge; j++)
                    {
                        var chargeMatches = matchList.GetCharge(j);
                        foreach (var match in chargeMatches)
                        {
                            foreach (var ionType in _ionTypes)
                            {
                                var prob = match.ContainsCleavageIons(ionType, _defaultTolerance, _relativeIntensityThreshold);
                                offsetFrequencyFunctions[j-1].AddOffsetFrequency(prob);
                            }
                        }
                    }
                }

                // Print Offset Frequency tables to output file
                var outFile = _outFileName.Replace("@", name);
                for (int i = 0; i < _precursorCharge; i++)
                {
                    var outFileCharge = outFile.Replace("*", (i + 1).ToString());
                    var offsetFrequencies = new List<IonProbability>();
                    if (_combineCharges)
                        offsetFrequencies = offsetFrequencyFunctions[i].GetCombinedChargeTable(_ionTypeFactory).ToList();
                    else
                        offsetFrequencies = offsetFrequencyFunctions[i].IonProbabilityTable.ToList();
                    using (var finalOutputFile = new StreamWriter(outFileCharge))
                    {
                        foreach (var ionProbability in offsetFrequencies)
                        {
                            finalOutputFile.WriteLine("{0}\t{1}", ionProbability.Ion.Name, ionProbability.Probability);
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
            _precursorCharge = Convert.ToInt32(config.Contents["precursorcharge"]);
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

/*            _precursorOff = (config.Contents.ContainsKey("precursorfrequencies") &&
                             config.Contents["precursorfrequencies"].ToLower() == "true");*/

            _combineCharges = (config.Contents.ContainsKey("combinecharges") &&
                 config.Contents["combinecharges"].ToLower() == "true");

            _relativeIntensityThreshold = Convert.ToDouble(config.Contents["relativeintensitythreshold"]);

            // Read ion data
            var ionInfo = reader.GetNodes("ion").First();
            int totalCharges = Convert.ToInt32(ionInfo.Contents["totalcharges"]);
            var ionTypeStr = ionInfo.Contents["iontype"].Split(',');
            var ions = new BaseIonType[ionTypeStr.Length];
            for (int i = 0; i < ionTypeStr.Length; i++)
            {
                switch (ionTypeStr[i].ToLower())
                {
                    case "a":
                        ions[i] = BaseIonType.A;
                        break;
                    case "b":
                        ions[i] = BaseIonType.B;
                        break;
                    case "c":
                        ions[i] = BaseIonType.C;
                        break;
                    case "x":
                        ions[i] = BaseIonType.X;
                        break;
                    case "y":
                        ions[i] = BaseIonType.Y;
                        break;
                    case "z":
                        ions[i] = BaseIonType.Z;
                        break;
                }
            }
            var ionLossStr = ionInfo.Contents["losses"].Split(',');
            var ionLosses = new NeutralLoss[ionLossStr.Length];
            for (int i = 0; i < ionLossStr.Length; i++)
            {
                switch (ionLossStr[i].ToLower())
                {
                    case "noloss":
                        ionLosses[i] = NeutralLoss.NoLoss;
                        break;
                    case "nh3":
                        ionLosses[i] = NeutralLoss.NH3;
                        break;
                    case "h2o":
                        ionLosses[i] = NeutralLoss.H2O;
                        break;
                }
            }
            _ionTypeFactory = new IonTypeFactory(ions, ionLosses, totalCharges);
            _ionTypes = _ionTypeFactory.GetAllKnownIonTypes().ToList();
            var tempIonList = new List<IonType>();
            if (ionInfo.Contents.ContainsKey("exclusions"))
            {
                var ionExclusions = ionInfo.Contents["exclusions"].Split(',');
                foreach (var ionType in _ionTypes)
                {
                    if (!ionExclusions.Contains(ionType.Name))
                        tempIonList.Add(ionType);
                }
            }
            _ionTypes = tempIonList;

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
