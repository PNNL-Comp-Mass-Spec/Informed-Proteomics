using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Config;
using InformedProteomics.Scoring.LikelihoodScoring.Data;
using InformedProteomics.Scoring.LikelihoodScoring.FileReaders;
using InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

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
        private double _relativeIntensityThreshold = 1.0;
        private bool _combineCharges;
        private bool _useDecoy;
        private int _precursorCharge;

        private readonly Tolerance _defaultTolerance = new Tolerance(0.5, ToleranceUnit.Mz);

        [Test]
        [Category("Local_Testing")]
        public void IonFrequencyFunction()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string configFilePath = @"C:\Users\wilk011\Documents\DataFiles\IonFreqConfig.ini";

            if (!File.Exists(configFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, configFilePath);
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

                var tableCount = 1;
                if (_precursorCharge > 0)
                {
                    tableCount = _precursorCharge;
                }

                var ionFrequencyFunctions = new IonFrequencyTable[tableCount];
                var decoyionFrequencyFunctions = new IonFrequencyTable[tableCount];
                for (var i = 0; i < tableCount; i++)
                {
                    ionFrequencyFunctions[i] = new IonFrequencyTable(_ionTypes,
                                _defaultTolerance, _relativeIntensityThreshold, _combineCharges);
                    decoyionFrequencyFunctions[i] = new IonFrequencyTable(_ionTypes,
                                _defaultTolerance, _relativeIntensityThreshold, _combineCharges);
                }

                for (var i = 0; i < txtFiles.Count; i++)
                {
                    var textFile = txtFiles[i];
                    var rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", Path.GetFileName(textFile), Path.GetFileName(rawFile));
                    var lcms = new LazyLcMsRun(rawFile, NoiseFiltration, NoiseFiltration);
                    var matchList = new SpectrumMatchList(lcms, txtFiles[i], DataFileFormat.IcBottomUp);
                    SpectrumMatchList decoyMatchList = null;
                    if (_useDecoy)
                    {
                        decoyMatchList = new SpectrumMatchList(lcms, txtFiles[i], DataFileFormat.IcBottomUp, 0, true);
                    }

                    for (var j = 0; j < tableCount; j++)
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
                for (var i = 0; i < tableCount; i++)
                {
                    var outFileName = outFile.Replace("*", (i + 1).ToString(CultureInfo.InvariantCulture));
                    var ionFrequencies = ionFrequencyFunctions[i].GetProbabilities();
                    var decoyIonFrequencies = decoyionFrequencyFunctions[i].GetProbabilities();

                    using var finalOutputFile = new StreamWriter(outFileName);

                    finalOutputFile.Write("Ion\tTarget");
                    if (_useDecoy)
                    {
                        finalOutputFile.Write("\tDecoy");
                    }

                    finalOutputFile.WriteLine();
                    for (var j = 0; j < ionFrequencies.Length; j++)
                    {
                        finalOutputFile.Write("{0}\t{1}", ionFrequencies[j].Label.Name, ionFrequencies[j].Prob);
                        if (_useDecoy)
                        {
                            finalOutputFile.Write("\t{0}", decoyIonFrequencies[j]);
                        }
                        finalOutputFile.WriteLine();
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

            _combineCharges = (config.Contents.ContainsKey("combinecharges") &&
                 string.Equals(config.Contents["combinecharges"], "true", StringComparison.OrdinalIgnoreCase));

            _useDecoy = (config.Contents.ContainsKey("usedecoy") &&
                string.Equals(config.Contents["usedecoy"], "true", StringComparison.OrdinalIgnoreCase));

            _relativeIntensityThreshold = Convert.ToDouble(config.Contents["relativeintensitythreshold"]);

            // Read ion data
            var ionInfo = reader.GetNodes("ion").First();
            var totalCharges = Convert.ToInt32(ionInfo.Contents["totalcharges"]);
            var ionTypeStr = ionInfo.Contents["iontype"].Split(',');
            var ions = new BaseIonType[ionTypeStr.Length];
            for (var i = 0; i < ionTypeStr.Length; i++)
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
            for (var i = 0; i < ionLossStr.Length; i++)
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
