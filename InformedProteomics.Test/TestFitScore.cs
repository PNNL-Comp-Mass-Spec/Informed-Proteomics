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
    public class TestFitScore
    {
        private string[] _names;
        private string _preTsv;
        private string _preRaw;
        private string _outPre;
        private string _outFileName;
        private int _intensityBins;
        private double _noiseFiltration;

        private List<IonType> _prefixionTypes;
        private List<IonType> _suffixionTypes; 
        private IonTypeFactory _ionTypeFactory;
        private ActivationMethod _act;
        private double _relativeIntensityThreshold = 1.0;
        private ScoreMethod _method;

        private readonly Tolerance _defaultTolerance = new Tolerance(15, ToleranceUnit.Ppm);

        private readonly double[] _binEdges =
        {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95};

        [Test]
        public void CosineScore()
        {
            InitTest(new ConfigFileReader(@"C:\Users\wilk011\Documents\DataFiles\CosineScoreConfig.ini"));

            foreach (var name in _names)
            {
                var tsvName = _preTsv.Replace("@", name);
                var rawName = _preRaw.Replace("@", name);
                var txtFiles = Directory.GetFiles(tsvName).ToList();
                var rawFilesTemp = Directory.GetFiles(rawName).ToList();
                var rawFiles = rawFilesTemp.Where(rawFile => Path.GetExtension(rawFile) == ".raw").ToList();

                Assert.True(rawFiles.Count == txtFiles.Count);

                var prefixTable = new ScoreTable(_method, _intensityBins, _binEdges);
                var suffixTable = new ScoreTable(_method, _intensityBins, _binEdges);

                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", Path.GetFileName(textFile), Path.GetFileName(rawFile));
                    var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, _noiseFiltration, _noiseFiltration);
                    var matchList = (new SpectrumMatchList(lcms, new TsvFileParser(txtFiles[i]), _act)).Matches;
                    prefixTable.AddMatches(matchList, _prefixionTypes.ToArray(), _defaultTolerance, _relativeIntensityThreshold);
                    suffixTable.AddMatches(matchList, _suffixionTypes.ToArray(), _defaultTolerance, _relativeIntensityThreshold);
                }

                var prefixhistograms = prefixTable.Histograms;
                var suffixhistograms = suffixTable.Histograms;
                var prefixEdges = prefixTable.IntensityBins;
                var suffixEdges = suffixTable.IntensityBins;
                var decoyprefixTable = new ScoreTable(_method, prefixEdges, _binEdges);
                var decoysuffixTable = new ScoreTable(_method, suffixEdges, _binEdges);

                for (int i = 0; i < txtFiles.Count; i++)
                {
                    string textFile = txtFiles[i];
                    string rawFile = rawFiles[i];
                    Console.WriteLine("{0}\t{1}", Path.GetFileName(textFile), Path.GetFileName(rawFile));
                    var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, _noiseFiltration, _noiseFiltration);
                    var decoyList = (new SpectrumMatchList(lcms, new TsvFileParser(txtFiles[i]), _act, true)).Matches;
                    decoyprefixTable.AddMatches(decoyList, _prefixionTypes.ToArray(), _defaultTolerance,
                        _relativeIntensityThreshold);
                    decoysuffixTable.AddMatches(decoyList, _suffixionTypes.ToArray(), _defaultTolerance,
                        _relativeIntensityThreshold);
                }

                var decoyprefixhistograms = decoyprefixTable.Histograms;
                var decoysuffixhistograms = decoysuffixTable.Histograms;

                var outFileName = _outFileName.Replace("@", name);
                var outFileNameIon = outFileName.Replace("*", "B");
                PrintOutput(outFileNameIon, prefixEdges, prefixhistograms, decoyprefixhistograms);

                outFileName = _outFileName.Replace("@", name);
                outFileNameIon = outFileName.Replace("*", "Y");
                PrintOutput(outFileNameIon, suffixEdges, suffixhistograms, decoysuffixhistograms);
            }
        }

        // Print tables to output file
        void PrintOutput(string fileName, double[] edges,
                List<Histogram<FitScore>> targethist, List<Histogram<FitScore>> decoyhist)
        {
            using (var outFile = new StreamWriter(fileName))
            {
                outFile.WriteLine("Intensity\tScore\tTarget\tDecoy");
                for (int i = 0; i < _intensityBins; i++)
                {
                    var prefixfrequencies = targethist[i].Frequencies;
                    var decoyprefixfrequencies = decoyhist[i].Frequencies;
                    for (int j = 0; j < prefixfrequencies.Count; j++)
                    {
                        outFile.WriteLine("{0}\t{1}\t{2}\t{3}", edges[i], _binEdges[j],
                            prefixfrequencies[j].Found, decoyprefixfrequencies[j].Found);
                    }
                }
            }
        }

        // Read Configuration file
        private void InitTest(ConfigFileReader reader)
        {
            // Read program variables
            var config = reader.GetNodes("vars").First();
            _intensityBins = Convert.ToInt32(config.Contents["intensitybins"]);
            _noiseFiltration = Convert.ToDouble(config.Contents["noisefiltration"]);
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

            var scoreMethod = config.Contents["activationmethod"].ToLower();
            switch (scoreMethod)
            {
                case "cosine":
                    _method = ScoreMethod.Cosine;
                    break;
                case "fitscore":
                    _method = ScoreMethod.FitScore;
                    break;
                case "pearson":
                    _method = ScoreMethod.Pearson;
                    break;
            }

            _relativeIntensityThreshold = Convert.ToDouble(config.Contents["relativeintensitythreshold"]);

            // Read ion configuration
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
            var ionTypes = _ionTypeFactory.GetAllKnownIonTypes();

            _prefixionTypes = new List<IonType>();
            _suffixionTypes = new List<IonType>();
            foreach (var ionType in ionTypes)
            {
                if (ionType.Name[0] == 'a' || ionType.Name[0] == 'b' || ionType.Name[0] == 'c')
                    _prefixionTypes.Add(ionType);
                else if (ionType.Name[0] == 'x' || ionType.Name[0] == 'y' || ionType.Name[0] == 'z')
                    _suffixionTypes.Add(ionType);
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
