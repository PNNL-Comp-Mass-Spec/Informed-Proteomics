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
    public class TestRankProbability
    {
        private string[] _dataSets;
        private string _preTsv;
        private string _preData;
        private string _outPre;
        private string _dataFormat;
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
        
        private List<IonType> _ionTypes;
        private double _selectedIonThreshold;
        private IonTypeFactory _ionTypeFactory;
        private ActivationMethod _act;
        private int _precursorCharge;
        private const double BinWidth = 1.005;

        private readonly Tolerance _defaultTolerancePpm = new Tolerance(10, ToleranceUnit.Ppm);
        private readonly Tolerance _defaultToleranceTh = new Tolerance(0.5, ToleranceUnit.Th);
        private Tolerance _tolerance;
        private readonly Tolerance _massErrorTolerance = new Tolerance(0.5, ToleranceUnit.Th);

        private int _retentionCount;
        private double _windowWidth;
        private double _precursorOffsetThreshold;

        private const double MassErrorBinWidth = 0.01;
        private const double MassErrorWidth = 0.5;
        private const double MassErrorStartPoint = 0.005;

        [Test]
        public void RankProbability()
        {
            InitTest(new ConfigFileReader(@"\\protoapps\UserData\Wilkins\BottomUp\RankProbabilityConfig.ini"));

            foreach (var dataSet in _dataSets)
            {
                var tsvName = _preTsv.Replace("@", dataSet);
                var rawName = _preData.Replace("@", dataSet);
                var txtFiles = new List<string>();
                var dataFiles = new List<string>();

                if (_dataFormat == "raw")
                {
                    // Read directory
                    var txtFilesTemp = Directory.GetFiles(tsvName).ToList();
                    txtFiles = txtFilesTemp.Where(txtFile => Path.GetExtension(txtFile) == ".tsv").ToList();
                    var dataFilesTemp = Directory.GetFiles(rawName).ToList();
                    dataFiles =
                        dataFilesTemp.Where(dataFile => Path.GetExtension(dataFile) == "." + _dataFormat).ToList();
                    Assert.True(dataFiles.Count == txtFiles.Count);
                }
                else
                {
                    dataFiles.Add(_preData+dataSet+"."+_dataFormat);
                }

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
                var massErrorTables = new MassErrorTable[_precursorCharge, _ionTypes.Count];
                var dMassErrorTables = new MassErrorTable[_precursorCharge, _ionTypes.Count];
                for (int chargeIndex = 0; chargeIndex < _precursorCharge; chargeIndex++)
                {
                    rankTables[chargeIndex] = new RankTable(_ionTypes.ToArray(), _tolerance, _maxRanks);
                    decoyRankTables[chargeIndex] = new RankTable(_ionTypes.ToArray(), _tolerance, _maxRanks);
                    ionFrequencyTables[chargeIndex] = new ProductIonFrequencyTable(_ionTypes, _tolerance,
                        _relativeIntensityThreshold);
                    offsetFrequencyTables[chargeIndex] = new List<PrecursorOffsetFrequencyTable>();
                    for (int j = 1; j <= (chargeIndex + 1); j++)
                    {
                        offsetFrequencyTables[chargeIndex].Add(new PrecursorOffsetFrequencyTable(_windowWidth/j, j, BinWidth/j));
                    }
                    for (int j = 0; j < _ionTypes.Count; j++)
                    {
                        massErrorTables[chargeIndex, j] = new MassErrorTable(new[] { _ionTypes[j] }, _massErrorTolerance,
                                                                             MassErrorWidth, MassErrorBinWidth, MassErrorStartPoint);
                        dMassErrorTables[chargeIndex, j] = new MassErrorTable(new[] { _ionTypes[j] }, _massErrorTolerance,
                                                        MassErrorWidth, MassErrorBinWidth, MassErrorStartPoint);
                    }
                }

                // Read files
                for (int i = 0; i < dataFiles.Count; i++)
                {
                    var matchList = new SpectrumMatchList(_act, false, _precursorCharge);
                    var decoyList = new SpectrumMatchList(_act, true, _precursorCharge);

                    if (_dataFormat == "raw")
                    {
                        string textFile = txtFiles[i];
                        string rawFile = dataFiles[i];
                        Console.WriteLine("{0}\t{1}", textFile, rawFile);
                        var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0);
                        matchList.AddMatchesFromTsvFile(lcms, new TsvFileParser(txtFiles[i]));
                    }
                    else
                    {
                        Console.WriteLine(dataSet+"."+_dataFormat);
                        matchList.AddMatchesFromMgfFile(dataFiles[i]);
                    }

                    decoyList.AddMatches(matchList);

                    for (int chargeIndex = 0; chargeIndex < _precursorCharge; chargeIndex++)
                    {
                        var chargeMatches = matchList.GetCharge(chargeIndex + 1);
                        var dChargeMatches = decoyList.GetCharge(chargeIndex + 1);

                        // Filter spectra
                        var filteredSpectra = new SpectrumMatchList(_act, false, _precursorCharge);
                        filteredSpectra.AddRange(chargeMatches);
                        filteredSpectra.FilterSpectra(_windowWidth, _retentionCount);
                        // Calculate ion probabilities
                        ionFrequencyTables[chargeIndex].AddMatches(filteredSpectra);

                        // Calculate precursor offset probabilities
                        foreach (var offsetFrequencyTable in offsetFrequencyTables[chargeIndex])
                        {
                            offsetFrequencyTable.AddMatches(chargeMatches);
                        }

                        // Calculate mass error tables
                        for (var j = 0; j < _ionTypes.Count; j++)
                        {
                            massErrorTables[chargeIndex, j].AddMatches(chargeMatches);
                            dMassErrorTables[chargeIndex, j].AddMatches(dChargeMatches);
                        }

                        // Initialize precursor filter
                        var precursorFilter = new PrecursorFilter(_precursorCharge, _tolerance);
                        precursorFilter.SetChargeOffsets(new PrecursorOffsets(offsetFrequencyTables[chargeIndex],
                            chargeIndex + 1,
                            _precursorOffsetThreshold));

                        // Filter precursor peaks
                        var precfilteredMatches = precursorFilter.FilterMatches(chargeMatches);
                        var dPrecFilteredMatches = precursorFilter.FilterMatches(dChargeMatches);
                        // Calculate rank probabilities
                        rankTables[chargeIndex].AddMatches(precfilteredMatches);
                        decoyRankTables[chargeIndex].AddMatches(dPrecFilteredMatches);
                    }
                }

                // Select ion types
                for (int chargeIndex = 0; chargeIndex < _precursorCharge; chargeIndex++)
                {
                    var selected = ionFrequencyTables[chargeIndex].SelectIons(_selectedIonThreshold);
                    foreach (var selectedIonProb in selected)
                    {
                        selectedIons[chargeIndex].Add(_ionTypeFactory.GetIonType(selectedIonProb.DataLabel.Name));
                    }
                    unselectedIons[chargeIndex] = _ionTypes.Except(selectedIons[chargeIndex]).ToList();
                }

                // Smooth ranks
                for (int chargeIndex = 0; chargeIndex < _precursorCharge; chargeIndex++)
                {
                    rankTables[chargeIndex].Smooth(2, 10, 30);
                    rankTables[chargeIndex].Smooth(3, 30, 50);
                    rankTables[chargeIndex].Smooth(5, 50, 90);
                    rankTables[chargeIndex].Smooth(7, 90, 120);
                    rankTables[chargeIndex].Smooth(10, 120);
                    decoyRankTables[chargeIndex].Smooth(5, 10, 20);
                    decoyRankTables[chargeIndex].Smooth(3, 20, 40);
                    decoyRankTables[chargeIndex].Smooth(7, 40, 90);
                    decoyRankTables[chargeIndex].Smooth(10, 90, 120);
                    decoyRankTables[chargeIndex].Smooth(15, 120);
                }

                // Write ion probability output files
                if (_writeIonProbabilities)
                {
                    var outFile = _ionProbabilityOutFileName.Replace("@", dataSet);
                    for (int chargeIndex = 0; chargeIndex < _precursorCharge; chargeIndex++)
                    {
                        string outFileName = outFile.Replace("*", (chargeIndex + 1).ToString(CultureInfo.InvariantCulture));
                        var file = new FileInfo(outFileName);
                        file.Directory.Create();
                        var ionProbabilities = ionFrequencyTables[chargeIndex].GetProbabilities();
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
                    var outFileName = _rankProbabilityOutFileName.Replace("@", dataSet);
                    for (int charge = 0; charge < _precursorCharge; charge++)
                    {
                        var chargeOutFileName = outFileName.Replace("*",
                            (charge + 1).ToString(CultureInfo.InvariantCulture));
                        var file = new FileInfo(chargeOutFileName);
                        file.Directory.Create();
                        WriteRankProbabilities(rankTables[charge], decoyRankTables[charge],
                                                selectedIons,
                                                charge, rankTables[charge].TotalRanks,
                                                chargeOutFileName);
                    }
                }

                // write precursor offset probability output files
                if (_writePrecursorOffsetProbabilities)
                {
                    var outFileName = _precursorOffsetProbabilityOutFileName.Replace("@", dataSet);
                    for (int charge = 0; charge < _precursorCharge; charge++)
                    {
                        var chargeOutFileName = outFileName.Replace("*",
                            (charge + 1).ToString(CultureInfo.InvariantCulture));
                        var file = new FileInfo(chargeOutFileName);
                        file.Directory.Create();
                        using (var outFile = new StreamWriter(chargeOutFileName))
                        {
                            for (int j = 0; j < offsetFrequencyTables[charge].Count; j++)
                            {
                                outFile.WriteLine("Charge\t{0}", j + 1);
                                var offsetprob = offsetFrequencyTables[charge][j].GetProbabilities();
                                foreach (var prob in offsetprob)
                                {
                                    var integerOffset = Math.Round(prob.DataLabel * (1 / BinWidth) * (j + 1));
                                    outFile.Write(integerOffset + "\t");
                                }
                                outFile.WriteLine();
                                foreach (var prob in offsetprob)
                                    outFile.Write(Math.Round(prob.Prob, 3) + "\t");
                                outFile.WriteLine();
                                outFile.WriteLine();
                            }
                        }
                    }
                }

                // Write mass error probability output files
                if (_writeMassErrorProbabilities)
                {
                    var outFileName = _massErrorProbabilityOutFileName.Replace("@", dataSet);
                    for (var chargeIndex = 0; chargeIndex < _precursorCharge; chargeIndex++)
                    {
                        var chargeOutFileName = outFileName.Replace("*",
                            (chargeIndex + 1).ToString(CultureInfo.InvariantCulture));
                        var file = new FileInfo(chargeOutFileName);
                        file.Directory.Create();
                        using (var outFile = new StreamWriter(chargeOutFileName))
                        {
                            outFile.Write("Error\t");
                            foreach (var selectedIon in selectedIons[chargeIndex])
                                outFile.Write("{0}\t{0}-De\t", selectedIon.Name);
                            outFile.Write("Unexplained\tUnexplained-De\tTotal");
                            outFile.WriteLine();

                            var selectedMassErrors = new List<Probability<double>[]>();
                            var dSelectedMassErrors = new List<Probability<double>[]>();
                            var unselectedMassErrors = new List<Probability<double>[]>();
                            var dUnselectedMassErrors = new List<Probability<double>[]>();
                            for (int i = 0; i < _ionTypes.Count; i++)
                            {
                                if (selectedIons[chargeIndex].Contains(_ionTypes[i]))
                                {
                                    selectedMassErrors.Add(massErrorTables[chargeIndex, i].GetProbabilities());
                                    dSelectedMassErrors.Add(dMassErrorTables[chargeIndex, i].GetProbabilities());
                                }
                                else
                                {
                                    unselectedMassErrors.Add(massErrorTables[chargeIndex, i].GetProbabilities());
                                    dUnselectedMassErrors.Add(dMassErrorTables[chargeIndex, i].GetProbabilities());
                                }
                            }
                            int massErrorLength = 0;
                            if (selectedMassErrors.Count > 0)
                                massErrorLength = selectedMassErrors[0].Length;
                            for (var i = 0; i < massErrorLength; i++)
                            {
                                outFile.Write(Math.Round(selectedMassErrors[0][i].DataLabel, 3) + "\t");
                                for (var j = 0; j < selectedIons[chargeIndex].Count; j++)
                                {
                                    var prob = selectedMassErrors[j][i].Found;
                                    var decoyProb = dSelectedMassErrors[j][i].Found;
                                    outFile.Write("{0}\t{1}\t", prob, decoyProb);
                                }
                                var probTotal = 0.0;
                                var decoyProbTotal = 0.0;
                                for (var j = 0; j < unselectedIons[chargeIndex].Count; j++)
                                {
                                    probTotal += unselectedMassErrors[j][i].Found;
                                    decoyProbTotal += dUnselectedMassErrors[j][i].Found;
                                }
                                outFile.Write("{0}\t{1}\t{2}", Math.Round(probTotal / unselectedIons.Length, 2),
                                                Math.Round(decoyProbTotal / unselectedIons.Length, 2),
                                                selectedMassErrors[0][1].Total);
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

                for (int i = 0; i < maxRanks + 1; i++)
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
                    _tolerance = _defaultTolerancePpm;
                    break;
                case "cid":
                    _act = ActivationMethod.CID;
                    _tolerance = _defaultToleranceTh;
                    break;
                case "etd":
                    _act = ActivationMethod.ETD;
                    _tolerance = _defaultTolerancePpm;
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
            _dataSets = fileInfo.Contents["name"].Split(',');
            _dataFormat = fileInfo.Contents["format"];
            _preTsv = fileInfo.Contents["tsvpath"];
            _preData = fileInfo.Contents["datapath"];
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
        }
    }
}