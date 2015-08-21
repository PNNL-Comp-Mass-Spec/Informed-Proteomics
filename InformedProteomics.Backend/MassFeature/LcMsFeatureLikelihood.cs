using System;
using System.IO;
using System.Reflection;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsFeatureLikelihood
    {
        public LcMsFeatureLikelihood(double likelihoodThreshold = 0)
        {
            _massBins = new double[30];
            var idx = 0;
            for (var m = 800; m <= 30000; m += 1000)
            {
                _massBins[idx] = m;
                idx++;
            }

            _distScoreTable = LoadTableFromResource("InformedProteomics.Backend.MassFeature.ScoringData.DistScore.tsv");
            _corrScoreTable = LoadTableFromResource("InformedProteomics.Backend.MassFeature.ScoringData.CorrScore.tsv");
            _intScoreTable = LoadTableFromResource("InformedProteomics.Backend.MassFeature.ScoringData.IntScore.tsv");

            _distScoreTableSummed = LoadTableFromResource("InformedProteomics.Backend.MassFeature.ScoringData.SummedDistScore.tsv");
            _corrScoreTableSummed = LoadTableFromResource("InformedProteomics.Backend.MassFeature.ScoringData.SummedCorrScore.tsv");
            _intScoreTableSummed = LoadTableFromResource("InformedProteomics.Backend.MassFeature.ScoringData.SummedIntScore.tsv");

            _abuScoreTable = LoadTableFromResource("InformedProteomics.Backend.MassFeature.ScoringData.AbuScore.tsv");
            _xicScoreTable1 = LoadTableFromResource("InformedProteomics.Backend.MassFeature.ScoringData.XicCorrScore1.tsv");
            _xicScoreTable2 = LoadTableFromResource("InformedProteomics.Backend.MassFeature.ScoringData.XicCorrScore2.tsv");

            ScoreThreshold = likelihoodThreshold;
        }

     
        public readonly double ScoreThreshold;
        
        public double GetScore(LcMsPeakCluster feature)
        {
            var mi = (int)Math.Round((feature.Mass - _massBins[0]) / (_massBins[1] - _massBins[0]));
            mi = (int) Math.Min(Math.Max(mi, 0), _massBins.Length - 1);
            var score = 0d;
            var abundance = feature.AbundanceDistributionAcrossCharge;

            for (var i = 0; i < 2; i++)
            {
                //score += _chargeScoreTable[mi][charge - 1];
                var abuScore = abundance[i];
                var k = (int)Math.Min(Math.Max(Math.Round(abuScore / 0.001), 0), NumberOfBins - 1);
                score += _abuScoreTable[mi][k];

                //if (!(abuScore > 0)) continue;
                var distScore = Math.Min(feature.EnvelopeDistanceScoreAcrossCharge[i], 1.0d);
                k = (int)Math.Min(Math.Max(Math.Round(distScore / 0.001), 0), NumberOfBins - 1);
                score += _distScoreTableSummed[mi][k];

                var corrScore = Math.Min(feature.EnvelopeCorrelationScoreAcrossCharge[i], 1.0d);
                k = (int)Math.Min(Math.Max(Math.Round(corrScore / 0.001), 0), NumberOfBins - 1);
                score += _corrScoreTableSummed[mi][k];

                var intScore = Math.Min(feature.EnvelopeIntensityScoreAcrossCharge[i], 1.0d);
                k = (int)Math.Min(Math.Max(Math.Round(intScore / 0.001), 0), NumberOfBins - 1);
                score += _intScoreTableSummed[mi][k];

                distScore = Math.Min(feature.BestDistanceScoreAcrossCharge[i], 1.0d);
                k = (int)Math.Min(Math.Max(Math.Round(distScore / 0.001), 0), NumberOfBins - 1);
                score += _distScoreTable[mi][k];

                corrScore = Math.Min(feature.BestCorrelationScoreAcrossCharge[i], 1.0d);
                k = (int)Math.Min(Math.Max(Math.Round(corrScore / 0.001), 0), NumberOfBins - 1);
                score += _corrScoreTable[mi][k];

                intScore = Math.Min(feature.BestIntensityScoreAcrossCharge[i], 1.0d);
                k = (int)Math.Min(Math.Max(Math.Round(intScore / 0.001), 0), NumberOfBins - 1);
                score += _intScoreTable[mi][k];

                var xicScore = Math.Min(feature.XicCorrelationBetweenBestCharges[i], 1.0d);
                k = (int)Math.Min(Math.Max(Math.Round(xicScore / 0.001), 0), NumberOfBins - 1);
                score += (i == 0) ? _xicScoreTable1[mi][k] : _xicScoreTable2[mi][k];
            }
            return score;
        }

        private const int NumberOfBins = 1001;

        private double[][] LoadTableFromResource(string resourceName)
        {
            //if (!File.Exists(fname)) throw new FileNotFoundException("Missing score datafile: " + fname);
            //var parser = new TsvFileParser(fname);
            var assembly = Assembly.GetExecutingAssembly();
            var textStreamReader = new StreamReader(assembly.GetManifestResourceStream(resourceName));
            var table = new double[_massBins.Length][];

            for (var i = 0; i < _massBins.Length; i++)
            {
                table[i] = new double[NumberOfBins];
                var line = textStreamReader.ReadLine();
                var token = line.Split('\t');
                for (var k = 0; k < NumberOfBins; k++)
                {
                    //var colData = parser.GetData(string.Format("{0}", k));
                    //table[i][k] = double.Parse(colData[i]);
                    table[i][k] = double.Parse(token[k]);
                }
            }
            return table;
        }
        
        private double[][] LoadTable(string fname)
        {
            if (!File.Exists(fname)) throw new FileNotFoundException("Missing score datafile: " + fname);
            
            var parser = new TsvFileParser(fname);
            var table = new double[_massBins.Length][];
            
            for (var i = 0; i < _massBins.Length; i++)
            {
                table[i] = new double[NumberOfBins];

                for (var k = 0; k < NumberOfBins; k++)
                {
                    var colData = parser.GetData(string.Format("{0}", k));
                    table[i][k] = double.Parse(colData[i]);
                }
            }
            return table;
        }
        
        private readonly double[][] _distScoreTable;
        private readonly double[][] _corrScoreTable;
        private readonly double[][] _intScoreTable;

        private readonly double[][] _distScoreTableSummed;
        private readonly double[][] _corrScoreTableSummed;
        private readonly double[][] _intScoreTableSummed;

        private readonly double[][] _xicScoreTable1;
        private readonly double[][] _xicScoreTable2;

        private readonly double[][] _abuScoreTable;

        private readonly double[] _massBins;


    }
}
