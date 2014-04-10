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
    class TestRankPeakProbability
    {
        private string[] _names;
        private string _preTsv;
        private string _preRaw;
        private string _outPre;
        private int _maxRanks;
        private double _relativeIntensityThreshold;

        private bool _writeIonProbabilities;
        private string _ionProbabilityOutFileName;

        private bool _writeRankProbabilities;
        private string _rankProbabilityOutFileName;

        private bool _writePrecursorOffsetProbabilities;
        private string _precursorOffsetProbabilityOutFileName;

        private bool _writeMassErrorProbabilities;
        private string _massErrorProbabilityOutFileName;

        private bool _writeIonPairProbabilities;
        private string _ionPairProbabilityOutFileName;
        
        private List<IonType> _ionTypes;
        private double _selectedIonThreshold;
        private List<IonType>[] _selectedIons;
        private List<IonType>[] _unselectedIons; 
        private IonTypeFactory _ionTypeFactory;
        private ActivationMethod _act;
        private int _precursorCharge;
        private const double BinWidth = 1.005;

        private readonly Tolerance _defaultTolerance = new Tolerance(15, ToleranceUnit.Ppm);

        private int _retentionCount;
        private int _searchWidth;
        private double _precursorOffsetThreshold;

        [Test]
        public void RankPeakProbability()
        {
            // read configuration settings
            InitTest(new ConfigFileReader(@"\\protoapps\UserData\Wilkins\BottomUp\RankPeakProbabilityConfig.ini"));

            foreach (var name in _names)
            {
                var tsvName = _preTsv.Replace("@", name);
                var rawName = _preRaw.Replace("@", name);
                var txtFiles = Directory.GetFiles(tsvName).ToList();
                var rawFilesTemp = Directory.GetFiles(rawName).ToList();
                var rawFiles = rawFilesTemp.Where(rawFile => Path.GetExtension(rawFile) == ".raw").ToList();

                Assert.True(rawFiles.Count == txtFiles.Count);

                // Initialize probability tables
                var rankTables = new RankTable[_precursorCharge];
                var ionFrequencyTables = new CleavageIonFrequencyTable[_precursorCharge];
                var offsetFrequencyTables = new List<PrecursorOffsetFrequencyTable>[_precursorCharge];
                for (int i = 0; i < _precursorCharge; i++)
                {
                    rankTables[i] = new RankTable(_ionTypes.ToArray(), _maxRanks);
                    ionFrequencyTables[i] = new CleavageIonFrequencyTable(_ionTypes, _defaultTolerance, _relativeIntensityThreshold);
                    offsetFrequencyTables[i] = new List<PrecursorOffsetFrequencyTable>();
                    for (int j = 1; j <= (i + 1); j++)
                    {
                        offsetFrequencyTables[i].Add(new PrecursorOffsetFrequencyTable(_searchWidth / j, j, BinWidth / j));
                    }
                }

                // Read files
                var matchList = new SpectrumMatchList(_act, false, _precursorCharge);
                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", textFile, rawFile);
                    var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0);
                    matchList.AddMatchesFromFile(lcms, new TsvFileParser(txtFiles[i]));
                }

                // Calculate probability tables
                for (int j = 0; j < _precursorCharge; j++)
                {
                    // Calculate ion probabilities
                    var chargeMatches = matchList.GetCharge(j + 1);
                    chargeMatches.FilterSpectra(_searchWidth, _retentionCount);
                    ionFrequencyTables[j].AddMatches(chargeMatches);

                    // Calculate precursor offset probabilities
                    for (int i = 0; i < _precursorCharge; i++)
                    {
                        foreach (var offsetFrequencyTable in offsetFrequencyTables[i])
                        {
                            offsetFrequencyTable.AddMatches(chargeMatches);
                        }
                    }
                }

                // Select ion types
                for (int j = 0; j < _precursorCharge; j++)
                {
                    var selected = ionFrequencyTables[j].SelectIons(_selectedIonThreshold);
                    foreach (var selectedIonProb in selected)
                    {
                        _selectedIons[j].Add(_ionTypeFactory.GetIonType(selectedIonProb.DataLabel));
                    }
                    _unselectedIons[j] = _ionTypes.Except(_selectedIons[j]).ToList();
                }

                // Create Mass Error tables and Ion Pair Frequency tables
                var selectedMassErrors = new List<List<Probability<double>>>[_precursorCharge];
                var selectedIonPairFrequencies = new List<List<Probability<IonPairFound>>>[_precursorCharge];
                var unselectedMassErrors = new List<List<Probability<double>>>[_precursorCharge];
                var unselectedIonPairFrequencies = new List<List<Probability<IonPairFound>>>[_precursorCharge];
                for (int i = 0; i < _precursorCharge; i++)
                {
                    var chargeMatches = matchList.GetCharge(i + 1);
                    // create tables for selected ions
                    selectedMassErrors[i] = new List<List<Probability<double>>>();
                    selectedIonPairFrequencies[i] = new List<List<Probability<IonPairFound>>>();
                    var selectedMassErrorTables = new MassErrorTable[_selectedIons[i].Count];
                    for (var j = 0; j < _selectedIons[i].Count; j++)
                    {
                        selectedMassErrorTables[j] = new MassErrorTable(new[] { _selectedIons[i][j] }, _defaultTolerance);
                        selectedMassErrorTables[j].AddMatches(chargeMatches);
                        selectedMassErrors[i].Add(selectedMassErrorTables[j].MassError);
                        selectedIonPairFrequencies[i].Add(selectedMassErrorTables[j].IonPairFrequency);
                    }
                    // create tables for unselected ions
                    unselectedMassErrors[i] = new List<List<Probability<double>>>();
                    unselectedIonPairFrequencies[i] = new List<List<Probability<IonPairFound>>>();
                    var unselectedMassErrorTables = new MassErrorTable[_unselectedIons[i].Count];
                    for (var j = 0; j < _unselectedIons[i].Count; j++)
                    {
                        unselectedMassErrorTables[j] = new MassErrorTable(new[] { _unselectedIons[i][j] }, _defaultTolerance);
                        unselectedMassErrorTables[j].AddMatches(chargeMatches);
                        unselectedMassErrors[i].Add(unselectedMassErrorTables[j].MassError);
                        unselectedIonPairFrequencies[i].Add(unselectedMassErrorTables[j].IonPairFrequency);
                    }
                }

                // Initialize precursor filter
                var precursorFilter = new PrecursorFilter(_precursorCharge, _defaultTolerance);
                for (int j = 0; j < _precursorCharge; j++)
                {
                    precursorFilter.SetChargeOffsets(new PrecursorOffsets(offsetFrequencyTables[j], j + 1, _precursorOffsetThreshold));
                }

                // Calculate rank probabilities
                for (int j = 0; j < _precursorCharge; j++)
                {
                    var chargeMatches = precursorFilter.FilterMatches(matchList.GetCharge(j+1));
                    rankTables[j].RankMatches(chargeMatches, _defaultTolerance);
                }

                // Write ion probability output files
                if (_writeIonProbabilities)
                {
                    var outFile = _ionProbabilityOutFileName.Replace("@", name);
                    for (int i = 0; i < _precursorCharge; i++)
                    {
                        string outFileName = outFile.Replace("*", (i + 1).ToString(CultureInfo.InvariantCulture));
                        var ionProbabilities = ionFrequencyTables[i].IonProbabilityTable;
                        using (var finalOutputFile = new StreamWriter(outFileName))
                        {
                            finalOutputFile.WriteLine("Ion\tTarget");
                            foreach (var ionProbability in ionProbabilities)
                            {
                                finalOutputFile.WriteLine("{0}\t{1}", ionProbability.DataLabel, ionProbability.Prob);
                            }
                        }
                    }
                }

                // Write rank probability output files
                if (_writeRankProbabilities)
                {
                    var outFileName = _rankProbabilityOutFileName.Replace("@", name);
                    for (int charge = 0; charge < _precursorCharge; charge++)
                    {
                        var chargeOutFileName = outFileName.Replace("*",
                            (charge + 1).ToString(CultureInfo.InvariantCulture));
                        var ionProbabilities = rankTables[charge].IonProbabilities;
                        WriteRankProbabilities(ionProbabilities, charge, rankTables[charge].TotalRanks,
                            chargeOutFileName);
                    }
                }

                // Write precursor offset probability output files
                if (_writePrecursorOffsetProbabilities)
                {
                    var outFileName = _precursorOffsetProbabilityOutFileName.Replace("@", name);
                    for (int charge = 0; charge < _precursorCharge; charge++)
                    {
                        var chargeOutFileName = outFileName.Replace("*",
                            (charge + 1).ToString(CultureInfo.InvariantCulture));
                        using (var outFile = new StreamWriter(chargeOutFileName))
                        {
                            for (int i = 0; i < offsetFrequencyTables[charge].Count; i++)
                            {
                                outFile.WriteLine("Charge\t{0}", i+1);
                                var offsetprob = offsetFrequencyTables[charge][i].OffsetFrequencies;
                                foreach (var prob in offsetprob)
                                {
                                    var integerOffset = Math.Round(prob.DataLabel * (1 / BinWidth)* (i+1));
                                    outFile.Write(integerOffset + "\t");
                                }
                                outFile.WriteLine();
                                foreach (var prob in offsetprob)
                                    outFile.Write(Math.Round(prob.Prob, 3)+"\t");
                                outFile.WriteLine();
                                outFile.WriteLine();
                            }
                        }
                    }
                }

                // Write mass error probability output files
                if (_writeMassErrorProbabilities)
                {
                    var outFileName = _massErrorProbabilityOutFileName.Replace("@", name);
                    for (var charge = 0; charge < _precursorCharge; charge++)
                    {
                        var chargeOutFileName = outFileName.Replace("*",
                            (charge + 1).ToString(CultureInfo.InvariantCulture));
                        using (var outFile = new StreamWriter(chargeOutFileName))
                        {
                            outFile.Write("Error\t");
                            foreach (var selectedIon in _selectedIons[charge]) outFile.Write(selectedIon.Name + "\t");
                            outFile.Write("Unexplained");
                            outFile.WriteLine();
                            var massErrorLength = selectedMassErrors[charge][0].Count;
                            for (var i = 0; i < massErrorLength; i++)
                            {
                                outFile.Write(Math.Round(selectedMassErrors[charge][0][i].DataLabel, 3)+"\t");
                                for (var j = 0; j < _selectedIons[charge].Count; j++)
                                {
                                    outFile.Write(selectedMassErrors[charge][j][i].Found + "\t");
                                }
                                var probTotal = 0;
                                for (var j = 0; j < _unselectedIons[charge].Count; j++)
                                {
                                    probTotal += unselectedMassErrors[charge][j][i].Found;
                                }
                                outFile.Write(Math.Round((double)probTotal / _unselectedIons.Length, 2));
                                outFile.WriteLine();
                            }
                        }
                    }
                }

                // Write ion pair probability table to output files
                if (_writeIonPairProbabilities)
                {
                    var outFileName = _ionPairProbabilityOutFileName.Replace("@", name);
                    for (var charge = 0; charge < _precursorCharge; charge++)
                    {
                        var chargeOutFileName = outFileName.Replace("*",
                            (charge + 1).ToString(CultureInfo.InvariantCulture));
                        using (var outFile = new StreamWriter(chargeOutFileName))
                        {
                            outFile.Write("Found\t");
                            foreach (var selectedIon in _selectedIons[charge]) outFile.Write(selectedIon.Name + "\t");
                            outFile.Write("Unexplained");
                            outFile.WriteLine();
                            var probLength = selectedIonPairFrequencies[charge][0].Count;
                            for (var i = 0; i < probLength; i++)
                            {
                                outFile.Write(selectedIonPairFrequencies[charge][0][i].DataLabel + "\t");
                                for (var j = 0; j < _selectedIons[charge].Count; j++)
                                {
                                    outFile.Write(Math.Round(selectedIonPairFrequencies[charge][j][i].Prob, 3) + "\t");
                                }
                                var probTotal = 0.0;
                                for (var j = 0; j < _unselectedIons[charge].Count; j++)
                                {
                                    probTotal += unselectedIonPairFrequencies[charge][j][i].Prob;
                                }
                                outFile.Write(Math.Round(probTotal / _unselectedIons.Length, 3));
                                outFile.WriteLine();
                            }
                        }
                    }
                }
            }
        }

        private void WriteRankProbabilities(IList<Dictionary<IonType, Probability<string>>> ionProbabilities, int charge, int totalRanks, string outFileCharge)
        {
            using (var outFile = new StreamWriter(outFileCharge))
            {
                // Write headers
                outFile.Write("Rank\t");
                foreach (var ionType in _selectedIons[charge])
                {
                    outFile.Write(ionType.Name + "\t");
                }
                outFile.Write("Unexplained");
                outFile.WriteLine();

                int maxRanks = _maxRanks;
                if (totalRanks < _maxRanks)
                    maxRanks = totalRanks;

                for (int i = 0; i < maxRanks; i++)
                {
                    outFile.Write(i + 1 + "\t");
                    // Write explained ion counts
                    for (int j = 0; j < _selectedIons[charge].Count; j++)
                    {
                        var ionFound = ionProbabilities[i][_selectedIons[charge][j]].Found;
                        outFile.Write("{0}\t", ionFound);
                    }

                    // Write average count for unexplained ions
                    int totalUnselected = 0;
                    int unselectedCount = 0;
                    for (int j = 0; j < _unselectedIons[charge].Count; j++)
                    {
                        totalUnselected += ionProbabilities[i][_unselectedIons[charge][j]].Found;
                        unselectedCount++;
                    }
                    outFile.Write(Math.Round((double)totalUnselected / unselectedCount, 2));
                    outFile.WriteLine();
                }
            }
        }

        // Read Configuration file
        private void InitTest(ConfigFileReader reader)
        {
            // Read program variables
            var config = reader.GetNodes("vars").First();
            _precursorCharge = Convert.ToInt32(config.Contents["precursorcharge"]);
            _precursorOffsetThreshold = Convert.ToDouble(config.Contents["precursoroffsetthreshold"]);
            _searchWidth = Convert.ToInt32(config.Contents["searchwidth"]);
            _retentionCount = Convert.ToInt32(config.Contents["retentioncount"]);
            _relativeIntensityThreshold = Convert.ToDouble(config.Contents["relativeintensitythreshold"]);
            _selectedIonThreshold = Convert.ToDouble(config.Contents["selectedionthreshold"]);
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

            _selectedIons = new List<IonType>[_precursorCharge];
            for (int i = 0; i < _precursorCharge; i++)
            {
                _selectedIons[i] = new List<IonType>();
            }
            _unselectedIons = new List<IonType>[_precursorCharge];
            for (int i = 0; i < _precursorCharge; i++)
            {
                _unselectedIons[i] = new List<IonType>();
            }

            _maxRanks = Convert.ToInt32(config.Contents["maxranks"]);

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

            _writeIonProbabilities = fileInfo.Contents.ContainsKey("ionprobabilityoutput");
            if (_writeIonProbabilities)
                _ionProbabilityOutFileName = _outPre + fileInfo.Contents["ionprobabilityoutput"];

            _writeRankProbabilities = fileInfo.Contents.ContainsKey("rankprobabilityoutput");
            if (_writeRankProbabilities)
                _rankProbabilityOutFileName = _outPre + fileInfo.Contents["rankprobabilityoutput"];

            _writePrecursorOffsetProbabilities = fileInfo.Contents.ContainsKey("precursoroffsetprobabilityoutput");
            if (_writePrecursorOffsetProbabilities)
                _precursorOffsetProbabilityOutFileName = _outPre + fileInfo.Contents["precursoroffsetprobabilityoutput"];

            _writeMassErrorProbabilities = fileInfo.Contents.ContainsKey("masserrorprobabilityoutput");
            if (_writeMassErrorProbabilities)
                _massErrorProbabilityOutFileName = _outPre + fileInfo.Contents["masserrorprobabilityoutput"];

            _writeIonPairProbabilities = fileInfo.Contents.ContainsKey("ionpairprobabilityoutput");
            if (_writeIonPairProbabilities)
                _ionPairProbabilityOutFileName = _outPre + fileInfo.Contents["ionpairprobabilityoutput"];
        }
    }
}
