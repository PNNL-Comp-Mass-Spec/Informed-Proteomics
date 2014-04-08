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
        private string _outFileName;
        private int _maxRanks;
        private double _relativeIntensityThreshold;

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
                        _selectedIons[j].Add(_ionTypeFactory.GetIonType(selectedIonProb.IonName));
                    }
                    _unselectedIons[j] = _ionTypes.Except(_selectedIons[j]).ToList();
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

                // Write output files
                var outFileName = _outFileName.Replace("@", name);
                for (int charge = 0; charge < _precursorCharge; charge++)
                {
                    var outFileCharge = outFileName.Replace("*", (charge + 1).ToString(CultureInfo.InvariantCulture));
                    var ionProbabilities = rankTables[charge].IonProbabilities;
                    WriteRankProbabilities(ionProbabilities, charge, rankTables[charge].TotalRanks, outFileCharge);
                }
            }
        }

        private void WriteRankProbabilities(List<Dictionary<IonType, IonProbability>> ionProbabilities, int charge, int totalRanks, string outFileCharge)
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
            var outFiletemp = fileInfo.Contents["outfile"];
            _outFileName = _outPre + outFiletemp;
        }
    }
}
