using System;
using System.Collections.Generic;
using System.Globalization;
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
    public class TestIonFrequency
    {
        private string[] _names;
        private string _preTsv;
        private string _preRaw;
        private string _outPre;
        private string _outFileName;
        private const double NoiseFiltration = 0;

        private List<IonType> _ionTypes;
        private IonTypeFactory _ionTypeFactory;
        private ActivationMethod _act;
        double _relativeIntensityThreshold = 1.0;
        private bool _combineCharges;
        private bool _useDecoy;
        private int _precursorCharge;

        private readonly Tolerance _defaultTolerance = new Tolerance(0.5, ToleranceUnit.Th);

        [Test]
        public void IonFrequencyFunction()
        {
            InitTest(new ConfigFileReader(@"C:\Users\wilk011\Documents\DataFiles\IonFreqConfig.ini"));

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

                var ionFrequencyFunctions = new IonFrequencyTable[tableCount];
                var decoyionFrequencyFunctions = new IonFrequencyTable[tableCount];
                for (int i = 0; i < tableCount; i++)
                {
                    ionFrequencyFunctions[i] = new ProductIonFrequencyTable(_ionTypes,
                                _defaultTolerance, _relativeIntensityThreshold, _combineCharges);
                    decoyionFrequencyFunctions[i] = new ProductIonFrequencyTable(_ionTypes,
                                _defaultTolerance, _relativeIntensityThreshold, _combineCharges);
                }

                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", Path.GetFileName(textFile), Path.GetFileName(rawFile));
                    var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, NoiseFiltration, NoiseFiltration);
                    var matchList = new SpectrumMatchList(lcms, new TsvFileParser(txtFiles[i]), _act);
                    SpectrumMatchList decoyMatchList = null;
                    if (_useDecoy)
                        decoyMatchList = new SpectrumMatchList(lcms, new TsvFileParser(txtFiles[i]), _act, true);

                    for (int j = 0; j < tableCount; j++)
                    {
                        var matches = _precursorCharge > 0 ? matchList.GetCharge(j + 1) : matchList;

                        matches.FilterSpectra();
                        ionFrequencyFunctions[j].AddMatches(matches);

                        if (decoyMatchList != null)
                        {
                            var decoyMatches = _precursorCharge > 0 ? decoyMatchList.GetCharge(j + 1) : decoyMatchList;
                            decoyMatches.FilterSpectra();
                            decoyionFrequencyFunctions[j].AddMatches(decoyMatches);
                        }
                    }
                }

                // Print Offset Frequency tables to output file
                var outFile = _outFileName.Replace("@", name);
                for (int i = 0; i < tableCount; i++)
                {
                    string outFileName = outFile.Replace("*", (i + 1).ToString(CultureInfo.InvariantCulture));
                    var ionFrequencies = ionFrequencyFunctions[i].GetProbabilities();
                    var decoyionFrequencies = decoyionFrequencyFunctions[i].GetProbabilities();
                    using (var finalOutputFile = new StreamWriter(outFileName))
                    {
                        finalOutputFile.Write("Ion\tTarget");
                        if (_useDecoy) finalOutputFile.Write("\tDecoy");
                        finalOutputFile.WriteLine();
                        for (int j=0;j<ionFrequencies.Length;j++)
                        {
                            finalOutputFile.Write("{0}\t{1}", ionFrequencies[j].DataLabel.Name, ionFrequencies[j].Prob);
                            if (_useDecoy)
                            {
                                finalOutputFile.Write("\t{0}", decoyionFrequencies[j]);
                            }
                            finalOutputFile.WriteLine();
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

            _combineCharges = (config.Contents.ContainsKey("combinecharges") &&
                 config.Contents["combinecharges"].ToLower() == "true");

            _useDecoy = (config.Contents.ContainsKey("usedecoy") &&
                config.Contents["usedecoy"].ToLower() == "true");

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
                tempIonList.AddRange(_ionTypes.Where(ionType => !ionExclusions.Contains(ionType.Name)));
                _ionTypes = tempIonList;
            }

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
