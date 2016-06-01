using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net;
using InformedProteomics.Backend.Data.Spectrometry;
using MathNet.Numerics.Distributions;

namespace InformedProteomics.TopDown.Scoring
{
    internal class DiscretizedNumber
    {
        internal DiscretizedNumber(double minBinCenter, double maxBinCenter, double binWidth)
        {
            BinWidth = binWidth;
            MaxBinCenter = maxBinCenter;
            MinBinCenter = minBinCenter;
            BinCount = (int)Math.Floor((MaxBinCenter - MinBinCenter + BinWidth)/BinWidth);
        }

        public int GetBinIndex(double val)
        {
            if (val >= MaxBinCenter) return BinCount - 1;
            if (val <= MinBinCenter) return 0;
            var binNum = (int)Math.Round((val - MinBinCenter) / BinWidth);
            return binNum;
        }

        public double GetBinCenter(int binNumber)
        {
            if (binNumber >= BinCount) return MaxBinCenter;
            if (binNumber <= 0) return MinBinCenter;
            return MinBinCenter + binNumber * BinWidth;
        }

        public double LowerBound { get { return MinBinCenter - BinWidth * 0.5; } }
        public double UpperBound { get { return MaxBinCenter + BinWidth * 0.5; } }

        public readonly double BinWidth;
        public readonly double MaxBinCenter;
        public readonly double MinBinCenter;
        public readonly int BinCount;
    }

    public class LikelihoodScoringModel
    {
        public LikelihoodScoringModel(string dataPath)
        {
            /*
            var scorePath = Path.Combine(dataPath, "IonTypeScore.txt");
            var reader = new StreamReader(scorePath);
            LoadIonTypeScoreTable(reader);
            reader.Close();

            scorePath = Path.Combine(dataPath, "MassScore.txt");
            reader = new StreamReader(scorePath);
            LoadMassScoreTable(reader);
            reader.Close();

            scorePath = Path.Combine(dataPath, "ChargeScore.txt");
            reader = new StreamReader(scorePath);
            LoadChargeScoreTable(reader);
            reader.Close();

            scorePath = Path.Combine(dataPath, "IntScore.txt");
            reader = new StreamReader(scorePath);
            LoadPeakScoreTable(reader, IntensityScoreTable);
            reader.Close();

            scorePath = Path.Combine(dataPath, "CorrScore.txt");
            reader = new StreamReader(scorePath);
            LoadPeakScoreTable(reader, CorrScoreTable);
            reader.Close();

            scorePath = Path.Combine(dataPath, "DistScore.txt");
            reader = new StreamReader(scorePath);
            LoadPeakScoreTable(reader, DistScoreTable);
            reader.Close();
            */
        }

        public double GetNodeScoreWithoutMassError(ActivationMethod activationMethod, bool isPrefix, DeconvolutedPeak matchedPeak, double refIntensity)
        {
            return GetNodeScoreWithoutMassError(activationMethod, isPrefix, matchedPeak.Mass, matchedPeak.Charge, matchedPeak.Corr, matchedPeak.Dist, matchedPeak.Intensity / refIntensity);
        }

        public double GetNodeScoreWithoutMassError(ActivationMethod activationMethod, bool isPrefix, double fragmentIonMass, int fragmentIonCharge,
            double corrScore, double distScore, double intensityScore)
        {
            var m = MassBinning.GetBinIndex(fragmentIonMass);
            var c = GetChargeIndex(fragmentIonCharge);
            var a = GetActivationMethodIndex(activationMethod);
            var t = GetIonTypeIndex(isPrefix);

            var ci = CorrBinning.GetBinIndex(corrScore);
            var di = DistBinning.GetBinIndex(distScore);
            var ii = IntensityBinning.GetBinIndex(intensityScore);

            var score = IonTypeScoreTable[a][t];
            score += MassScoreTable[a][t][m];
            score += ChargeScoreTable[m][c];
            score += CorrScoreTable[c][m][ci];
            score += DistScoreTable[c][m][di];
            score += IntensityScoreTable[c][m][ii];

            return score;
        }

        public double GetNodeScore(ActivationMethod activationMethod, bool isPrefix, double fragmentIonMass, DeconvolutedPeak matchedPeak, double refIntensity)
        {
            var massErrorPpm = 1e6 * ((matchedPeak.Mass - fragmentIonMass) / fragmentIonMass);
            return GetNodeScore(activationMethod, isPrefix, fragmentIonMass, matchedPeak.Charge, matchedPeak.Corr, matchedPeak.Dist, matchedPeak.Intensity / refIntensity, massErrorPpm);
        }

        public double GetNodeScore(ActivationMethod activationMethod, bool isPrefix, double fragmentIonMass, int fragmentIonCharge,
            double corrScore, double distScore, double intensityScore, double massErrorPpm)
        {
            var score = GetNodeScoreWithoutMassError(activationMethod, isPrefix, fragmentIonMass, fragmentIonCharge, corrScore,
                distScore, intensityScore);
            score += GetMassErrorScore(massErrorPpm);
            return score;
        }

        public double GetEdgeScore(ActivationMethod activationMethod, int chargeDiff, double massErrorPpm)
        {
            return 0d;
        }

        public double GetMassErrorScore(double massErrorPpm)
        {
            //y = normpdf(bin_Z, 5.0140, 3.1534);
            //y2 = normpdf(bin_Z, 6.3361, 4.0447);
            var p = Normal.PDF(5.014, 3.1534, massErrorPpm);
            var pnull = Normal.PDF(7.3361, 4.0447, massErrorPpm);
            return Math.Log(p/pnull);
        }

        public double GetEdgeScore(ActivationMethod activationMethod, int charge1, double mass1, int charge2, double mass2)
        {
            var mzErrorPpm = (Math.Abs(mass1 - mass2)/mass1)*1e6;
            return GetEdgeScore(activationMethod, Math.Abs(charge1 - charge2), mzErrorPpm);
        }

        public int GetActivationMethodIndex(ActivationMethod actMethod)
        {
            return (actMethod != ActivationMethod.ETD) ? 0 : 1;
        }

        public int GetActivationMethodIndex(string activationMethodName)
        {
            return activationMethodName.Equals("ETD") ? 1 : 0;
        }

        public int GetIonTypeIndex(bool isPrefix)
        {
            return isPrefix ? 0 : 1;
        }

        public int GetIonTypeIndex(BaseIonType ionType)
        {
            return GetIonTypeIndex(ionType.IsPrefix);
        }

        public int GetIonTypeIndex(string ionType)
        {
            return ionType.Equals("prefix") ? 0 : 1;
        }

        public int GetChargeDiffIndex(int chargeDiff)
        {
            if (chargeDiff < 0) return 0;
            if (chargeDiff >= ChargeBinLength) return ChargeBinLength - 1;
            return chargeDiff;
        }

        public int GetChargeIndex(int charge)
        {
            if (charge < 1) return 1;
            if (charge > ChargeBinLength) return ChargeBinLength - 1;
            return charge - 1;
        }


        private void LoadIonTypeScoreTable(StreamReader reader)
        {
            for (var i = 0; i < ActivationBinLength; i++)
            {
                var line = reader.ReadLine();
                if (string.IsNullOrEmpty(line)) break;

                var a = GetActivationMethodIndex(line);
                line = reader.ReadLine();
                var token = line.Split('\t');

                for (var j = 0; j < IonTypeBinLength; j++) IonTypeScoreTable[a][j] = GetBoundedLkScore(token[j]);
            }
        }

        private void LoadMassScoreTable(StreamReader reader)
        {
            for(var i = 0; i < ActivationBinLength; i++)
            {
                var line = reader.ReadLine();
                if (string.IsNullOrEmpty(line)) break;

                var tokens = line.Split(',');

                if (tokens.Length != 2) break;

                var a = GetActivationMethodIndex(tokens[0]);
                var t = GetIonTypeIndex(tokens[1]);

                line = reader.ReadLine();
                var token = line.Split('\t');

                for (var j = 0; j < MassBinning.BinCount; j++) MassScoreTable[a][t][j] = GetBoundedLkScore(token[j]);
            }
        }

        private void LoadChargeScoreTable(StreamReader reader)
        {
            while(true)
            {
                var line = reader.ReadLine();
                if (string.IsNullOrEmpty(line)) break;

                //var tokens = line.Split(',');
                var m = MassBinning.GetBinIndex(double.Parse(line));

                line = reader.ReadLine();
                var token = line.Split('\t');

                for (var j = 0; j < ChargeBinLength; j++) ChargeScoreTable[m][j] = GetBoundedLkScore(token[j]);
            }
        }

        private void LoadPeakScoreTable(StreamReader reader, double[][][] table)
        {
            var binLen = table[0][0].Length;
            while (true)
            {
                var line = reader.ReadLine();

                if (string.IsNullOrEmpty(line)) break;

                var tokens = line.Split(',');
                if (tokens.Length != 2) break;

                var c = GetChargeIndex(int.Parse(tokens[0]));
                var m = MassBinning.GetBinIndex(double.Parse(tokens[1]));

                line = reader.ReadLine();
                var token = line.Split('\t');

                for (var j = 0; j < binLen; j++) table[c][m][j] = GetBoundedLkScore(token[j]);
            }
        }

        private double GetBoundedLkScore(string scoreStr)
        {
            return Math.Max(Math.Min(double.Parse(scoreStr), ScoreUpperBound), ScoreLowerBound);
        }

        private const double ScoreLowerBound = -1.5;
        private const double ScoreUpperBound = 1.5;

        private const int ActivationBinLength = 2;
        private const int ChargeBinLength = 20;
        private const int IonTypeBinLength = 2;

        private static readonly DiscretizedNumber MassBinning;
        private static readonly DiscretizedNumber IntensityBinning;
        private static readonly DiscretizedNumber CorrBinning;
        private static readonly DiscretizedNumber DistBinning;
        //private static readonly DiscretizedNumber MzErrorBinning;

        private static readonly double[][] IonTypeScoreTable;
        private static readonly double[][][] MassScoreTable;

        private static readonly double[][] ChargeScoreTable;

        private static readonly double[][][] CorrScoreTable;
        private static readonly double[][][] DistScoreTable;
        private static readonly double[][][] IntensityScoreTable;

        static LikelihoodScoringModel()
        {
            MassBinning = new DiscretizedNumber(100, 15000, 100);
            IntensityBinning = new DiscretizedNumber(0.2, 30, 0.4);
            CorrBinning = new DiscretizedNumber(0.5, 1.0, 0.01);
            DistBinning = new DiscretizedNumber(0, 0.5, 0.01);
            //MzErrorBinning = new DiscretizedNumber(0.1, 20, 0.1);

            IonTypeScoreTable = new double[ActivationBinLength][];
            MassScoreTable = new double[ActivationBinLength][][];

            ChargeScoreTable = new double[MassBinning.BinCount][];
            for (var i = 0; i < MassBinning.BinCount; i++)
            {
                ChargeScoreTable[i] = new double[ChargeBinLength];
            }

            CorrScoreTable = new double[ChargeBinLength][][];
            DistScoreTable = new double[ChargeBinLength][][];
            IntensityScoreTable = new double[ChargeBinLength][][];

            for (var i = 0; i < ActivationBinLength; i++)
            {
                IonTypeScoreTable[i] = new double[IonTypeBinLength];
                ChargeScoreTable[i] = new double[MassBinning.BinCount];

                MassScoreTable[i] = new double[IonTypeBinLength][];
                for (var j = 0; j < IonTypeBinLength; j++)
                {
                    MassScoreTable[i][j] = new double[MassBinning.BinCount];
                }
            }

            for (var i = 0; i < ChargeBinLength; i++)
            {
                CorrScoreTable[i] = new double[MassBinning.BinCount][];
                DistScoreTable[i] = new double[MassBinning.BinCount][];
                IntensityScoreTable[i] = new double[MassBinning.BinCount][];

                for (var j = 0; j < MassBinning.BinCount; j++)
                {
                    CorrScoreTable[i][j] = new double[CorrBinning.BinCount];
                    DistScoreTable[i][j] = new double[DistBinning.BinCount];
                    IntensityScoreTable[i][j] = new double[IntensityBinning.BinCount];
                }
            }
        }
    }


}
