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
        private IonTypeFactory _ionTypeFactory;
        private ActivationMethod _act;
        private int _precursorCharge;
        private const double BinWidth = 1.005;

        private readonly Tolerance _defaultTolerance = new Tolerance(0.5, ToleranceUnit.Th);

        private int _retentionCount;
        private double _windowWidth;
        private double _precursorOffsetThreshold;

        private const double MassErrorBinWidth = 0.01;
        private const double MassErrorWidth = 0.5;
        private const double MassErrorStartPoint = 0.005;

        [Test]
        public void RankPeakProbability()
        {
            // read configuration settings
            InitTest(new ConfigFileReader(@"\\protoapps\UserData\Wilkins\BottomUp\RankPeakProbabilityConfig.ini"));

            foreach (var name in _names)
            {
                // Read directory
                var tsvName = _preTsv.Replace("@", name);
                var rawName = _preRaw.Replace("@", name);
                var txtFiles = Directory.GetFiles(tsvName).ToList();
                var rawFilesTemp = Directory.GetFiles(rawName).ToList();
                var rawFiles = rawFilesTemp.Where(rawFile => Path.GetExtension(rawFile) == ".raw").ToList();
                Assert.True(rawFiles.Count == txtFiles.Count);

                // Initialize ion lists
                var selectedIons = new List<IonType>[_precursorCharge];
                var unselectedIons = new List<IonType>[_precursorCharge];
                for (int i = 0; i < _precursorCharge; i++)
                {
                    selectedIons[i] = new List<IonType>();
                    unselectedIons[i] = new List<IonType>();
                }

                // Initialize probability tables
                var rankTables = new RankTable[_precursorCharge];
                var decoyRankTables = new RankTable[_precursorCharge];
                var ionFrequencyTables = new ProductIonFrequencyTable[_precursorCharge];
                var offsetFrequencyTables = new List<PrecursorOffsetFrequencyTable>[_precursorCharge];
                for (int i = 0; i < _precursorCharge; i++)
                {
                    rankTables[i] = new RankTable(_ionTypes.ToArray(), _defaultTolerance, _maxRanks);
                    decoyRankTables[i] = new RankTable(_ionTypes.ToArray(), _defaultTolerance, _maxRanks);
                    ionFrequencyTables[i] = new ProductIonFrequencyTable(_ionTypes, _defaultTolerance, _relativeIntensityThreshold);
                    offsetFrequencyTables[i] = new List<PrecursorOffsetFrequencyTable>();
                    for (int j = 1; j <= (i + 1); j++)
                    {
                        offsetFrequencyTables[i].Add(new PrecursorOffsetFrequencyTable(_windowWidth / j, j, BinWidth / j));
                    }
                }

                // Read files
                var matchList = new SpectrumMatchList(_act, false, _precursorCharge);
                var decoyList = new SpectrumMatchList(_act, true, _precursorCharge);
                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", textFile, rawFile);
                    var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0);
                    matchList.AddMatchesFromFile(lcms, new TsvFileParser(txtFiles[i]));
                }
                decoyList.AddMatches(matchList);

                // Calculate probability tables
                for (int j = 0; j < _precursorCharge; j++)
                {
                    // Calculate ion probabilities
                    var chargeMatches = matchList.GetCharge(j + 1);
                    chargeMatches.FilterSpectra(_windowWidth, _retentionCount);
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
                        selectedIons[j].Add(_ionTypeFactory.GetIonType(selectedIonProb.DataLabel.Name));
                    }
                    unselectedIons[j] = _ionTypes.Except(selectedIons[j]).ToList();
                }

                // Create Mass Error tables and Ion Pair Frequency tables
                var selectedMassErrors = new List<List<Probability<double>>>[_precursorCharge];
                var decoySelectedMassErrors = new List<List<Probability<double>>>[_precursorCharge];
                var selectedIonPairFrequencies = new List<List<Probability<IonPairFound>>>[_precursorCharge];
                var decoySelectedIonPairFrequencies = new List<List<Probability<IonPairFound>>>[_precursorCharge];
                var unselectedMassErrors = new List<List<Probability<double>>>[_precursorCharge];
                var decoyUnselectedMassErrors = new List<List<Probability<double>>>[_precursorCharge];
                var unselectedIonPairFrequencies = new List<List<Probability<IonPairFound>>>[_precursorCharge];
                var decoyUnselectedIonPairFrequencies = new List<List<Probability<IonPairFound>>>[_precursorCharge];
                for (int i = 0; i < _precursorCharge; i++)
                {
                    var chargeMatches = matchList.GetCharge(i + 1);
                    var decoyChargeMatches = decoyList.GetCharge(i + 1);
                    // create tables for selected ions
                    selectedMassErrors[i] = new List<List<Probability<double>>>();
                    decoySelectedMassErrors[i] = new List<List<Probability<double>>>();
                    selectedIonPairFrequencies[i] = new List<List<Probability<IonPairFound>>>();
                    decoySelectedIonPairFrequencies[i] = new List<List<Probability<IonPairFound>>>();
                    var selectedMassErrorTables = new MassErrorTable[selectedIons[i].Count];
                    var decoySelectedMassErrorTables = new MassErrorTable[selectedIons[i].Count];
                    for (var j = 0; j < selectedIons[i].Count; j++)
                    {
                        selectedMassErrorTables[j] = new MassErrorTable(new[] { selectedIons[i][j] }, _defaultTolerance,
                                                        MassErrorWidth, MassErrorBinWidth, MassErrorStartPoint);
                        decoySelectedMassErrorTables[j] = new MassErrorTable(new[] { selectedIons[i][j] }, _defaultTolerance,
                                                        MassErrorWidth, MassErrorBinWidth, MassErrorStartPoint);
                        selectedMassErrorTables[j].AddMatches(chargeMatches);
                        decoySelectedMassErrorTables[j].AddMatches(decoyChargeMatches);
                        selectedMassErrors[i].Add(selectedMassErrorTables[j].GetProbabilities().ToList());
                        decoySelectedMassErrors[i].Add(decoySelectedMassErrorTables[j].GetProbabilities().ToList());
                        selectedIonPairFrequencies[i].Add(selectedMassErrorTables[j].IonPairFrequency);
                        decoySelectedIonPairFrequencies[i].Add(decoySelectedMassErrorTables[j].IonPairFrequency);
                    }
                    // create tables for unselected ions
                    unselectedMassErrors[i] = new List<List<Probability<double>>>();
                    decoyUnselectedMassErrors[i] = new List<List<Probability<double>>>();
                    unselectedIonPairFrequencies[i] = new List<List<Probability<IonPairFound>>>();
                    decoyUnselectedIonPairFrequencies[i] = new List<List<Probability<IonPairFound>>>();
                    var unselectedMassErrorTables = new MassErrorTable[unselectedIons[i].Count];
                    var decoyUnselectedMassErrorTables = new MassErrorTable[unselectedIons[i].Count];
                    for (var j = 0; j < unselectedIons[i].Count; j++)
                    {
                        unselectedMassErrorTables[j] = new MassErrorTable(new[] { unselectedIons[i][j] }, _defaultTolerance,
                                                        MassErrorWidth, MassErrorBinWidth, MassErrorStartPoint);
                        decoyUnselectedMassErrorTables[j] = new MassErrorTable(new[] { unselectedIons[i][j] }, _defaultTolerance,
                                                        MassErrorWidth, MassErrorBinWidth, MassErrorStartPoint);
                        unselectedMassErrorTables[j].AddMatches(chargeMatches);
                        decoyUnselectedMassErrorTables[j].AddMatches(decoyChargeMatches);
                        unselectedMassErrors[i].Add(unselectedMassErrorTables[j].GetProbabilities().ToList());
                        decoyUnselectedMassErrors[i].Add(decoyUnselectedMassErrorTables[j].GetProbabilities().ToList());
                        unselectedIonPairFrequencies[i].Add(unselectedMassErrorTables[j].IonPairFrequency);
                        decoyUnselectedIonPairFrequencies[i].Add(decoyUnselectedMassErrorTables[j].IonPairFrequency);
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
                    var decoyChargeMatches = precursorFilter.FilterMatches(decoyList.GetCharge(j + 1));
                    rankTables[j].AddMatches(chargeMatches);
                    decoyRankTables[j].AddMatches(decoyChargeMatches);
                    
                    // smooth ranks
                    rankTables[j].Smooth(2, 10, 30);
                    rankTables[j].Smooth(3, 30, 50);
                    rankTables[j].Smooth(5, 50, 90);
                    rankTables[j].Smooth(7, 90, 120);
                    rankTables[j].Smooth(10, 120);
                    decoyRankTables[j].Smooth(5, 10, 20);
                    decoyRankTables[j].Smooth(3, 20, 40);
                    decoyRankTables[j].Smooth(7, 40, 90);
                    decoyRankTables[j].Smooth(10, 90, 120);
                    decoyRankTables[j].Smooth(15, 120);
                }

                // Write ion probability output files
                if (_writeIonProbabilities)
                {
                    var outFile = _ionProbabilityOutFileName.Replace("@", name);
                    for (int i = 0; i < _precursorCharge; i++)
                    {
                        string outFileName = outFile.Replace("*", (i + 1).ToString(CultureInfo.InvariantCulture));
                        var ionProbabilities = ionFrequencyTables[i].GetProbabilities();
                        using (var finalOutputFile = new StreamWriter(outFileName))
                        {
                            finalOutputFile.WriteLine("Ion\tTarget");
                            foreach (var ionProbability in ionProbabilities)
                            {
                                finalOutputFile.WriteLine("{0}\t{1}", ionProbability.DataLabel.Name, ionProbability.Prob);
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
                        WriteRankProbabilities(rankTables[charge], decoyRankTables[charge], 
                                               selectedIons,
                                               charge, rankTables[charge].TotalRanks,
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
                                var offsetprob = offsetFrequencyTables[charge][i].GetProbabilities();
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
                            foreach (var selectedIon in selectedIons[charge]) 
                                outFile.Write("{0}\t{0}-De\t", selectedIon.Name);
                            outFile.Write("Unexplained\tUnexplained-De\tTotal");
                            outFile.WriteLine();
                            var massErrorLength = selectedMassErrors[charge][0].Count;
                            for (var i = 0; i < massErrorLength; i++)
                            {
                                outFile.Write(Math.Round(selectedMassErrors[charge][0][i].DataLabel, 3)+"\t");
                                for (var j = 0; j < selectedIons[charge].Count; j++)
                                {
                                    var prob = selectedMassErrors[charge][j][i].Found;
                                    var decoyProb = decoySelectedMassErrors[charge][j][i].Found;
                                    outFile.Write("{0}\t{1}\t", prob, decoyProb);
                                }
                                var probTotal = 0.0;
                                var decoyProbTotal = 0.0;
                                for (var j = 0; j < unselectedIons[charge].Count; j++)
                                {
                                    probTotal += unselectedMassErrors[charge][j][i].Found;
                                    decoyProbTotal += decoyUnselectedMassErrors[charge][j][i].Found;
                                }
                                outFile.Write("{0}\t{1}\t{2}", Math.Round(probTotal / unselectedIons.Length, 2),
                                                Math.Round(decoyProbTotal/unselectedIons.Length, 2),
                                                selectedMassErrors[charge][0][1].Total);
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
                            foreach (var selectedIon in selectedIons[charge]) outFile.Write(selectedIon.Name + "\t");
                            outFile.Write("Unexplained");
                            outFile.WriteLine();
                            var probLength = selectedIonPairFrequencies[charge][0].Count;
                            for (var i = 0; i < probLength; i++)
                            {
                                outFile.Write(selectedIonPairFrequencies[charge][0][i].DataLabel + "\t");
                                for (var j = 0; j < selectedIons[charge].Count; j++)
                                {
                                    outFile.Write(Math.Round(selectedIonPairFrequencies[charge][j][i].Prob, 3) + "\t");
                                }
                                var probTotal = 0.0;
                                for (var j = 0; j < unselectedIons[charge].Count; j++)
                                {
                                    probTotal += unselectedIonPairFrequencies[charge][j][i].Prob;
                                }
                                outFile.Write(Math.Round(probTotal / unselectedIons.Length, 3));
                                outFile.WriteLine();
                            }
                        }
                    }
                }
            }
        }

        private void WriteRankProbabilities(RankTable targetRanks,
                                            RankTable decoyRanks,
                                            List<IonType>[] selectedIons, 
                                            int charge, int totalRanks, string outFileCharge)
        {
            var targetProb = targetRanks.GetProbabilities();
            var decoyProb = decoyRanks.GetProbabilities();
            using (var outFile = new StreamWriter(outFileCharge))
            {
                // Write headers
                outFile.Write("Rank\t");
                foreach (var ionType in selectedIons[charge])
                {
                    outFile.Write("{0}\t{0}-De\t", ionType.Name);
                }
                outFile.Write("Unexplained\tUnexplained-De\tTotal");
                outFile.WriteLine();

                int maxRanks = _maxRanks;
                if (totalRanks < _maxRanks)
                    maxRanks = totalRanks;

                for (int i = 0; i < maxRanks+1; i++)
                {
                    if (i == maxRanks)
                        outFile.Write("None" + "\t");
                    else
                        outFile.Write(i + 1 + "\t");

                    var totalUnselected = 0.0;
                    var unselectedCount = 0.0;
                    var decoyTotalUnselected = 0.0;
                    var decoyUnselectedCount = 0.0;
                    for (int j = 0; j < _ionTypes.Count; j++)
                    {
                        var currProb = targetProb[i, j];
                        var currDecoyProb = decoyProb[i, j];
                        var currIonType = currProb.DataLabel;
                        if (selectedIons[charge].Contains(currIonType))
                        {
                            var currProbFound = Math.Round(currProb.Found, 2);
                            var currDecoyProbFound = Math.Round(currDecoyProb.Found, 2);
                            outFile.Write("{0}\t{1}\t", currProbFound, currDecoyProbFound);
                        }
                        else
                        {
                            totalUnselected += currProb.Found;
                            decoyTotalUnselected += currDecoyProb.Found;
                            unselectedCount++;
                            decoyUnselectedCount++;
                        }
                    }
                    outFile.Write("{0}\t{1}\t{2}",
                                Math.Round(totalUnselected / unselectedCount, 2),
                                Math.Round(decoyTotalUnselected / decoyUnselectedCount, 2),
                                targetProb[i, 0].Total);
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
            _windowWidth = Convert.ToInt32(config.Contents["searchwidth"]);
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
