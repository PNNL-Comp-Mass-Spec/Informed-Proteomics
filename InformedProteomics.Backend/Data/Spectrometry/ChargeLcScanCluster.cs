using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public static class ChargeLcScanScore
    {
        public const int EnvelopeCorrelation = 0;
        public const int EnvelopeCorrelationSummed = 1;

        public const int RankSum = 2;
        public const int RankSumMedian = 3;

        public const int Poisson = 4;
        public const int PoissonMedian = 5;
        public const int PoissonSummed = 6;

        public const int BhattacharyyaDistance = 7;
        public const int BhattacharyyaDistanceSummed = 8;

        public const int KullbackLeiblerDivergence = 9;
        public const int KullbackLeiblerDivergenceSummed = 10;

        public const int MzError = 11;
        public const int MzErrorSummed = 12;
        public const int XicCorrelation = 13;

        public const int Count = 14;
    }

    public class ChargeLcScanCluster : Ms1Feature
    {
        public bool Active;
        public IsotopeList IsotopeList { get; private set; }

        public double Probability
        {
            get;
            set;
        }

        public override int MinScanNum
        {
            get
            {
                if ((ScanLength < 16) && MinCol - 1 >= 0) return _ms1ScanNums[MinCol - 1];
                
                return _ms1ScanNums[MinCol];
            }
        }

        public override int MaxScanNum
        {
            get
            {
                if ((ScanLength < 16) && MaxCol + 1 < _ms1ScanNums.Length) return _ms1ScanNums[MaxCol + 1];
                
                return _ms1ScanNums[MaxCol];
            }
        }

        public override int MinCharge { get { return _minCharge + MinRow; } }
        public override int MaxCharge { get { return _minCharge + MaxRow; } }

        public int ScanLength { get { return MaxCol - MinCol + 1; } }
        public int ChargeLength { get { return MaxRow - MinRow + 1; } }

        public int EnvelopeCount { get; internal set; }
     
        /*
        private static string[] _tsvHeaders = new string[]
        {
            "min_scan_num", "max_scan_num", "min_charge", "max_charge", "monoisotopic_mw", "rep_scan_num", "rep_charge", "rep_mz", "abundance",
            "summed_envelope_count", "scan_length",
            "envelope_corr", "summed_envelope_corr", 
            "ranksum_score", "avg_ranksum_score",
            "poisson_score","avg_poisson_score", "summed_poisson_score",
            "bc_distance", "summed_bc_distance", 
            "kl_divergence", "summed_kl_divergence", 
            "mz_error", "summed_mz_error", "xic_corr",
            "isotopic_envelope", "probability", "flag"
        };*/

        public static readonly string[] TsvHeaderWithScore = new string[]
        {
            "FeatureID", "MinScan", "MaxScan", "MinCharge", "MaxCharge", "MonoMass", "RepScan", "RepCharge", "RepMz", "Abundance",
            "ObservedEnvelopes", "ScanLength",
            "BestCorr", "SummedCorr", 
            "RanksumScore", "AvgRanksumScore",
            "PoissonScore","AvgPoissonScore", "SummedPoissonScore",
            "BcDist", "SummedBcDist", 
            "KlDiv", "SummedKlDiv", 
            "MzDiffPpm", "AvgMzDiff", "XicCorr",
            "Envelope", "Probability", "GoodEnough"
        };
        /*
        public static readonly string[] OldTsvHeaderWithScore = new string[]
        {
            "FeatureID", "min_scan_num", "max_scan_num", "min_charge", "max_charge", "monoisotopic_mw", "rep_scan_num", "rep_charge",
            "rep_mz", "abundance", "summed_envelope_count", "scan_length", "envelope_corr",
            "summed_envelope_corr", "ranksum_score", "avg_ranksum_score", "poisson_score",
            "avg_poisson_score", "summed_poisson_score", "bc_distance", "summed_bc_distance",
            "kl_divergence", "summed_kl_divergence", "mz_error",
            "summed_mz_error", "xic_corr", "isotopic_envelope", "probability", "good_enough"
        };
        
        public static readonly string[] OldTsvHeader = new string[]
        {
            "FeatureID", "min_scan_num", "max_scan_num", "min_charge", "max_charge", "monoisotopic_mw", "rep_scan_num", "rep_charge", "rep_mz", "abundance", 
            "envelope_corr","summed_envelope_corr", 
            "mz_error", "xic_corr", 
            "isotopic_envelope", "probability", "good_enough"
        };
        */
        private static readonly string[] TsvHeader = new string[]
        {
            "FeatureID", "MinScan", "MaxScan", "MinCharge", "MaxCharge", "MonoMass", "RepScan", "RepCharge", "RepMz", "Abundance",
            "BestCorr", "SummedCorr", 
            "MzDiffPpm", "XicCorr", 
            "Envelope", "Probability", "GoodEnough"
        };

        public static string GetHeaderString(bool withScore = false)
        {
            return (withScore) ? ArrayUtil.ToString(TsvHeaderWithScore) : ArrayUtil.ToString(TsvHeader);        
        }

        public override string ToString()
        {
            return GetString();
        }

        public string GetString(bool withScore = false, bool withEnvelope = true)
        {
            var intensity = GetSummedEnvelope();

            var sb = new StringBuilder(string.Format("{0}\t{1}\t{2}\t{3}\t{4:0.0000}\t{5}\t{6}\t{7:0.0000}\t{8:0.00}",
                                        MinScanNum, MaxScanNum, MinCharge, MaxCharge, RepresentativeMass, RepresentativeScanNum, RepresentativeCharge,
                                        RepresentativeMz, intensity.Sum()));
            
            if (withScore)
            {
                sb.AppendFormat("\t{0}", EnvelopeCount);
                sb.AppendFormat("\t{0}", ScanLength);
                foreach (var score in _scores) sb.AppendFormat("\t{0}", score);
            }
            else
            {
                sb.AppendFormat("\t{0:0.00000}", GetScore(ChargeLcScanScore.EnvelopeCorrelation));
                sb.AppendFormat("\t{0:0.00000}", GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed));
                sb.AppendFormat("\t{0:0.00000}", GetScore(ChargeLcScanScore.MzError));
                sb.AppendFormat("\t{0:0.00000}", GetScore(ChargeLcScanScore.XicCorrelation));
            }

            if (withEnvelope)
            {
                sb.Append("\t");
                var maxIntensity = intensity.Max();
                for (var i = 0; i < intensity.Length; i++)
                {
                    if (i != 0) sb.Append(";");
                    sb.AppendFormat("{0},{1:0.000}", IsotopeList[i].Index, intensity[i] / maxIntensity);
                }                
            }

            sb.AppendFormat("\t{0:0.0000}\t{1}", Probability, GoodEnough ? 1 : 0);

            return sb.ToString();
        }

        public double Abundance
        {
            get { return MemberEnvelope.Sum(e => e.Sum()); }
        }

        public double GetScore(int scoreType)
        {
            return _scores[scoreType];
        }
        
        private static readonly Tolerance MergeTolerance = new Tolerance(4);
        public static bool Merge(ChargeLcScanCluster c1, ChargeLcScanCluster c2)
        {
            var scanLen = Math.Max(c1.ScanLength, c2.ScanLength);
            var overlapLength = 0;
            if (c1.MinCol <= c2.MinCol && c2.MinCol < c1.MaxCol)
            {
                if (c2.MaxCol <= c1.MaxCol) overlapLength = scanLen;
                else overlapLength = c1.MaxCol - c2.MinCol + 1;
            }
            else if (c2.MinCol <= c1.MinCol && c1.MinCol < c2.MaxCol)
            {
                if (c1.MaxCol <= c2.MaxCol) overlapLength = scanLen;
                else overlapLength = c2.MaxCol - c1.MinCol + 1;
            }
            else
            {
                return false;
            }

            if (overlapLength < 0.75*scanLen) return false;

            if (Math.Abs(c1.RepresentativeMass - c2.RepresentativeMass) > MergeTolerance.GetToleranceAsTh(c1.RepresentativeMass))
                return false;

            ChargeLcScanCluster src = null;
            ChargeLcScanCluster dest = null;
            if (c1.GoodEnough == c2.GoodEnough)
            {
                if (c1.Probability > c2.Probability)
                {
                    dest = c1;
                    src = c2;
                }
                else
                {
                    dest = c2;
                    src = c1;
                }
            }
            else
            {
                if (c1.GoodEnough)
                {
                    dest = c1;
                    src = c2;
                }
                else
                {
                    dest = c2;
                    src = c1;
                }                
            }

            if (dest.Active == false)
            {
                return Merge(src, dest._mergedTo);
            }
            
            //Console.Write("{0:0.00} {1} {2}", src.RepresentativeMass, src.MinScanNum, src.MaxScanNum);
            //Console.WriteLine("  ->>  {0:0.00} {1} {2}", dest.RepresentativeMass, dest.MinScanNum, dest.MaxScanNum);

            src.Active = false;
            src._mergedTo = dest;

            dest.Active = true;
            dest._mergedTo = null;
            return true;
        }

        private ChargeLcScanCluster _mergedTo;
        public double ClusteringScore { get; private set; }
        internal List<ChargeLcScanCell> Members { get; private set; }
        internal List<double[]> MemberEnvelope { get; private set; }
        internal double[] ClusteringEnvelope { get; private set; }
        internal int MaxCol { get; private set; }
        internal int MinCol { get; private set; }
        internal int MaxRow { get; private set; }
        internal int MinRow { get; private set; }

        private readonly double[] _scores;
        private readonly int[] _ms1ScanNums;
        private readonly int _minCharge;
        internal double[] SummedEnvelope;


        public double ClusteringScore2;

        internal ChargeLcScanCluster(int minCharge, int[] ms1ScanNums, IsotopeList isoList)
        {
            ClusteringScore2 = 0; 

            Active = true;
            _ms1ScanNums = ms1ScanNums;
            _minCharge = minCharge;

            IsotopeList = isoList;
            Members = new List<ChargeLcScanCell>();
            MemberEnvelope = new List<double[]>();
            ClusteringEnvelope = new double[IsotopeList.Count];
            ClusteringScore = 0;

            RepresentativeMass = 0d;
            RepresentativeCharge = 0;
            RepresentativeMz = 0d;
            RepresentativeScanNum = 0;

            _scores = new double[ChargeLcScanScore.Count];
            MaxCol = 0;
            MinCol = _ms1ScanNums.Length - 1;
            MaxRow = 0;
            MinRow = 99999;
            SummedEnvelope = null;
            EnvelopeCount = 0;

            _mergedTo = null;
        }

        public void SetScore(int scoreType, double scoreValue)
        {
            _scores[scoreType] = scoreValue;
        }

        internal void AddMember(ChargeLcScanCell member, double[] newObservedEnvelope, double updatedClusteringScore = -1)
        {
            if (member.Col > MaxCol) MaxCol = member.Col;
            if (member.Col < MinCol) MinCol = member.Col;
            if (member.Row > MaxRow) MaxRow = member.Row;
            if (member.Row < MinRow) MinRow = member.Row;

            Members.Add(member);
            MemberEnvelope.Add(newObservedEnvelope);

            if (updatedClusteringScore < 0) return;
            for (var i = 0; i < ClusteringEnvelope.Length; i++) ClusteringEnvelope[i] += newObservedEnvelope[i];
            ClusteringScore = updatedClusteringScore;
        }

        internal void ClearMember()
        {
            MaxCol = -1;
            MinCol = _ms1ScanNums.Length - 1;
            MaxRow = -1;
            MinRow = 99999;
            Members.Clear();
            MemberEnvelope.Clear();
        }

        internal void SetBoundary(IEnumerable<int> selectedMembers)
        {
            MaxCol = 0;
            MinCol = _ms1ScanNums.Length - 1;
            MaxRow = 0;
            MinRow = 99999;
            var tempMembers = Members;
            var tempMemberEnvelop = MemberEnvelope;
            Members = new List<ChargeLcScanCell>();
            MemberEnvelope = new List<double[]>();

            foreach(var i in selectedMembers)
            {
                var member = tempMembers[i];
                if (member.Col > MaxCol) MaxCol = member.Col;
                if (member.Col < MinCol) MinCol = member.Col;
                if (member.Row > MaxRow) MaxRow = member.Row;
                if (member.Row < MinRow) MinRow = member.Row;

                Members.Add(member);
                MemberEnvelope.Add(tempMemberEnvelop[i]);
            }

        }

        internal void SetRepresentativeMass(double mass, double mz, int charge, int scanNum)
        {
            RepresentativeMass = mass;
            RepresentativeMz = mz;
            RepresentativeScanNum = scanNum;
            RepresentativeCharge = charge;
        }

        internal double[] GetSummedEnvelope()
        {
            if (SummedEnvelope != null) return SummedEnvelope;

            var ret = new double[IsotopeList.Count];
            foreach (var e in MemberEnvelope)
            {
                for (var i = 0; i < ret.Length; i++) ret[i] += e[i];
            }

            return ret;
        }

        private static readonly double[] LogisticRegressionBetaVector = new double[]
        {
          -9.80322248576848,
            0.0177345049138546,
            0.0182034423264027,
            1.34377709499630,
            1.22056790970138,
            -0.0121680010258230,
            0.212804311237213,
            0.0438596045871836,
            -0.0537174416494930,
            0.0627816486522023,
            11.7161939412475,
            97.3958203478190,
            -28.4567577808012,
            -87.3263081327653,
            -0.0594784589790504,
            0.0707772936171733,
            2.25809384761825
        };
        
        internal double GetProbabilityByLogisticRegression()
        {
            var eta = LogisticRegressionBetaVector[0];

            eta += EnvelopeCount * LogisticRegressionBetaVector[1];
            eta += ScanLength * LogisticRegressionBetaVector[2];

            for (var i = 3; i < LogisticRegressionBetaVector.Length; i++)
                eta += _scores[i - 3] * LogisticRegressionBetaVector[i];
                
            var pi = Math.Exp(eta);
            return pi/(pi + 1);
        }

        internal bool TooBad
        {
            get
            {
                if (RepresentativeMass < 2000)
                {
                    if (GetScore(ChargeLcScanScore.EnvelopeCorrelation) < 0.5) return true;
                    if (GetScore(ChargeLcScanScore.RankSum) < 4.3219) return true;
                    if (GetScore(ChargeLcScanScore.Poisson) < 4.3219) return true;
                    return false;
                }

                if (EnvelopeCount < 10)
                {
                    if (GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) < 0.5) return true;
                    if (GetScore(ChargeLcScanScore.BhattacharyyaDistance) > 0.08) return true;

                    if (GetScore(ChargeLcScanScore.RankSum) < 9 && GetScore(ChargeLcScanScore.Poisson) < 15 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.01) return true;

                    if (GetScore(ChargeLcScanScore.RankSum) < 8 && GetScore(ChargeLcScanScore.Poisson) < 15 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.005) return true;

                    if (GetScore(ChargeLcScanScore.RankSum) < 7 && GetScore(ChargeLcScanScore.Poisson) < 10)
                        return true;

                    if (GetScore(ChargeLcScanScore.RankSum) < 6 || GetScore(ChargeLcScanScore.Poisson) < 6.6439)
                        return true;

                    if (GetScore(ChargeLcScanScore.RankSum) < 8 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.005)
                        return true;

                    if (GetScore(ChargeLcScanScore.Poisson) < 15 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.005)
                        return true;

                    if (GetScore(ChargeLcScanScore.MzError) > 5) return true;
                    if (GetScore(ChargeLcScanScore.MzErrorSummed) > 9) return true;

                    if (GetScore(ChargeLcScanScore.XicCorrelation) < 0.4 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.005) return true;

                    if (GetScore(ChargeLcScanScore.PoissonSummed) < 10 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.005) return true;
                }
                else
                {
                    if (GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) < 0.6) return true;
                    if (GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.08) return true;
                    if (GetScore(ChargeLcScanScore.BhattacharyyaDistance) > 1) return true;
                    if (GetScore(ChargeLcScanScore.EnvelopeCorrelation) < 0.4) return true;
                    if (GetScore(ChargeLcScanScore.RankSum) < 4.3219) return true;
                    if (GetScore(ChargeLcScanScore.Poisson) < 6.6439) return true;

                    if (GetScore(ChargeLcScanScore.RankSum) < 6 && GetScore(ChargeLcScanScore.Poisson) < 20 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.015 &&
                        GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) < 0.9) return true;

                    if (GetScore(ChargeLcScanScore.RankSum) < 8 && GetScore(ChargeLcScanScore.Poisson) < 18 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.015 &&
                        GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) < 0.9) return true;

                    if (GetScore(ChargeLcScanScore.RankSum) < 6 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.01 &&
                        GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) < 0.95) return true;

                    if (GetScore(ChargeLcScanScore.Poisson) < 10 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.01 &&
                        GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) < 0.95) return true;

                    if (GetScore(ChargeLcScanScore.RankSumMedian) < 2 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.01 &&
                        GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) < 0.95) return true;

                    if (GetScore(ChargeLcScanScore.RankSumMedian) < 4.5 &&
                        GetScore(ChargeLcScanScore.PoissonMedian) < 6 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.01 &&
                        GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) < 0.95) return true;

                    if (GetScore(ChargeLcScanScore.RankSumMedian) > 7 &&
                        GetScore(ChargeLcScanScore.PoissonMedian) < 8 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.01 &&
                        GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) < 0.95) return true;

                    if (GetScore(ChargeLcScanScore.PoissonMedian) < 9 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.015 &&
                        GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) < 0.9) return true;

                    if (GetScore(ChargeLcScanScore.PoissonSummed) < 16) return true;

                    if (GetScore(ChargeLcScanScore.PoissonSummed) < 45 &&
                        GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) > 0.01) return true;
                }
                return false;
            }
        }

        internal bool GoodEnough
        {
            get
            {
               if (RepresentativeMass < 13000)
               {
                   return GetScore(ChargeLcScanScore.EnvelopeCorrelation) > 0.8 &&
                          GetScore(ChargeLcScanScore.RankSum) > 6 &&
                          GetScore(ChargeLcScanScore.RankSumMedian) > 3 &&
                          GetScore(ChargeLcScanScore.Poisson) > 12 &&
                          GetScore(ChargeLcScanScore.PoissonMedian) > 8 &&
                          GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) < 0.01 &&
                          GetScore(ChargeLcScanScore.MzError) < 3 &&
                          GetScore(ChargeLcScanScore.MzErrorSummed) < 5 &&
                          GetScore(ChargeLcScanScore.XicCorrelation) > 0.2;
                }
                else
                {
                    return GetScore(ChargeLcScanScore.EnvelopeCorrelation) > 0.6 &&
                           GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed) > 0.95 &&
                           GetScore(ChargeLcScanScore.RankSum) > 8 &&
                           GetScore(ChargeLcScanScore.RankSumMedian) > 3 &&
                           GetScore(ChargeLcScanScore.Poisson) > 10 &&
                           GetScore(ChargeLcScanScore.PoissonMedian) > 8 &&
                           GetScore(ChargeLcScanScore.BhattacharyyaDistanceSummed) < 0.01 &&
                           GetScore(ChargeLcScanScore.MzError) < 3 &&
                           GetScore(ChargeLcScanScore.MzErrorSummed) < 5;
                }
            }
        }
    }


    class ChargeLcScanCell : Tuple<int, int>
    {
        internal ChargeLcScanCell(int row, int col) : base(row, col) { }
        internal int Row { get { return Item1; } }
        internal int Col { get { return Item2; } }
    }
}
