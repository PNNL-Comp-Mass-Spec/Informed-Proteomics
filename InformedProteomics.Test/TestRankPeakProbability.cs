using System;
using System.Collections.Generic;
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

        private List<IonType> _ionTypes;
        private IonTypeFactory _ionTypeFactory;
        private ActivationMethod _act;
        private int _precursorCharge;

        private readonly Tolerance _defaultTolerance = new Tolerance(15, ToleranceUnit.Ppm);

        [Test]
        public void RankPeakProbability()
        {
            InitTest(new ConfigFileReader(@"\\protoapps\UserData\Wilkins\BottomUp\RankPeakProbabilityConfig.ini"));

            foreach (var name in _names)
            {
                var tsvName = _preTsv.Replace("@", name);
                var rawName = _preRaw.Replace("@", name);
                var txtFiles = Directory.GetFiles(tsvName).ToList();
                var rawFilesTemp = Directory.GetFiles(rawName).ToList();
                var rawFiles = rawFilesTemp.Where(rawFile => Path.GetExtension(rawFile) == ".raw").ToList();

                Assert.True(rawFiles.Count == txtFiles.Count);

                var rankTables = new RankTable[_precursorCharge];
                for (int i = 0; i < _precursorCharge; i++)
                {
                    rankTables[i] = new RankTable(_ionTypes.ToArray());
                }

                // Calculate Rank Probabilities
                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", textFile, rawFile);
                    var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0);
                    var matchList = new SpectrumMatchList(lcms, new TsvFileParser(txtFiles[i]), _act);
                    for (int j = 0; j < _precursorCharge; j++)
                    {
                        var chargeMatches = matchList.GetCharge(j+1);
                        rankTables[j].RankMatches(chargeMatches, _defaultTolerance);
                    }
                }

                // Write Rank Probability Table to file
                var outFileName = _outFileName.Replace("@", name);
                for (int charge = 0; charge < _precursorCharge; charge++)
                {
                    var outFileCharge = outFileName.Replace("*", (charge + 1).ToString());
                    var rankProbabilites = rankTables[charge].RankProbabilities;
                    using (var outFile = new StreamWriter(outFileCharge))
                    {
                        outFile.Write("Rank\t");
                        for (int i = 0; i < _ionTypes.Count; i++)
                        {
                            outFile.Write(_ionTypes[i].Name + "\t");
                        }
                        outFile.Write("Unexplained");
                        outFile.WriteLine();

                        int maxRanks = _maxRanks;
                        if (rankTables[charge].TotalRanks < _maxRanks)
                            maxRanks = rankTables[charge].TotalRanks;

                        for (int i = 0; i < maxRanks; i++)
                        {
                            double totalRankProbability = 0.0;
                            outFile.Write(i+1 + "\t");
                            for (int j = 0; j < _ionTypes.Count; j++)
                            {
                                var ionProbability = rankProbabilites[i, j].Prob;
                                totalRankProbability += ionProbability;
                                outFile.Write("{0}\t", ionProbability);
                            }
                            outFile.Write(1-totalRankProbability);
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
                foreach (var ionType in _ionTypes)
                {
                    if (!ionExclusions.Contains(ionType.Name))
                        tempIonList.Add(ionType);
                }
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
