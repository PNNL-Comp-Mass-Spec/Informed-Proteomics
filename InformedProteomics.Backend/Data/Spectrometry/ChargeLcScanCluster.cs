using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public static class ChargeLcScanScore
    {
        public const int EnvelopeCorrelation = 0;
        public const int EnvelopeCorrelationSummed = 1;
        public const int RankSum = 2;
        public const int RankSumSummed = 3;
        public const int HyperGeometric = 4;
        public const int HyperGeometricSummed = 5;
        

        public const int FinalScore = 6;

        public const int Count = 7;
    }
    
    public class ChargeLcScanCluster : Ms1Feature
    {
        public bool Active;
        public IsotopeList IsotopeList { get; private set; }
        public override int MaxScanNum { get { return _ms1ScanNums[MaxCol]; } }
        public override int MinScanNum { get { return _ms1ScanNums[MinCol]; } }
        public override int MinCharge { get { return _minCharge + MinRow; } }
        public override int MaxCharge { get { return _minCharge + MaxRow; } }

        public int ScanLength { get { return MaxCol - MinCol + 1; } }
        public int ChargeLength { get { return MaxRow - MinRow + 1; } }

        public static string GetHeaderString()
        {
            return string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}",
                            "min_scan_num", "max_scan_num", "min_charge", "max_charge", 
                            "monoisotopic_mw", "rep_scan_num", "rep_charge", "rep_mz", 
                            "abundance",
                            "envelope_corr", "summed_envelope_corr", 
                            "ranksum_score", "summed_ranksum_score",
                            "hypergeometric_score", "summed_hypergeometric_score", "isotopic_envelope");
        }
        
        public override string ToString()
        {
            //var intensity = GetSummedEnvelope();
            //var summedCorr = IsotopeList.GetPearsonCorrelation(intensity);
            var intensity = SummedEnvelope;
            var sb = new StringBuilder(
                        string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9:0.00000}\t{10:0.00000}\t{11:0.000}\t{12:0.000}\t{13:0.000}\t{14:0.000}\t",
                        MinScanNum, MaxScanNum, MinCharge, MaxCharge,
                        RepresentativeMass, RepresentativeScanNum, RepresentativeCharge, RepresentativeMz, 
                        Abundance,
                        GetScore(ChargeLcScanScore.EnvelopeCorrelation), GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed),
                        GetScore(ChargeLcScanScore.RankSum), GetScore(ChargeLcScanScore.RankSumSummed),
                        GetScore(ChargeLcScanScore.HyperGeometric), GetScore(ChargeLcScanScore.HyperGeometricSummed)));

            var maxIntensity = intensity.Max();
            for (var i = 0; i < intensity.Length; i++)
            {
                if (i != 0) sb.Append(";");
                sb.AppendFormat("{0},{1:0.000}", IsotopeList[i].Index, intensity[i] / maxIntensity);
            }
            return sb.ToString();
        }

        public double Abundance
        {
            get { return MemberEnvelope.Sum(m => m.Sum()); }
        }

        public double GetScore(int scoreType = ChargeLcScanScore.FinalScore)
        {
            return _scores[scoreType];
        }

        public double[] SummedEnvelope { get; internal set; }
        /*
        public double[] GetSummedEnvelope()
        {
            var ret = new double[IsotopeList.Count];
            for (var i = 0; i < IsotopeList.Count; i++)
            {
                foreach (var m in MemberEnvelope)
                    ret[i] += m[i];
            }
            return ret;
        }
        */
        public bool Overlaps(ChargeLcScanCluster other)
        {
            const int colBuf = 1;
            const int rowBuf = 0;

            var ret = ((MinCol - colBuf <= other.MinCol && other.MinCol <= MaxCol + colBuf) ||
                       (MinCol - colBuf <= other.MaxCol && other.MaxCol <= MaxCol + colBuf)) &&
                       ((MinRow - rowBuf <= other.MinRow && other.MinRow <= MaxRow + rowBuf) ||
                       (MinRow - rowBuf <= other.MaxRow && other.MaxRow <= MaxRow + rowBuf));
            return ret;
        }

        public void Merge(ChargeLcScanCluster other)
        {
            other.Active = false;

            if (other.MaxCol > MaxCol) MaxCol = other.MaxCol;
            if (other.MinCol < MinCol) MinCol = other.MinCol;
            if (other.MaxRow > MaxRow) MaxRow = other.MaxRow;
            if (other.MinRow < MinRow) MinRow = other.MinRow;

            for (var i = 0; i < ClusteringEnvelope.Length; i++) ClusteringEnvelope[i] += other.ClusteringEnvelope[i];

            Members.AddRange(other.Members);
            MemberEnvelope.AddRange(other.MemberEnvelope);
        }

        public double ClusteringScore { get; private set; }

        public double[] GetInputValueForSvm()
        {
            var svmFeatures = new double[]
            {
                MinCharge,
                MaxCharge,
                GetScore(ChargeLcScanScore.EnvelopeCorrelation),
                GetScore(ChargeLcScanScore.EnvelopeCorrelationSummed),
                GetScore(ChargeLcScanScore.HyperGeometric),
                GetScore(ChargeLcScanScore.HyperGeometricSummed),
                GetScore(ChargeLcScanScore.RankSum),
                GetScore(ChargeLcScanScore.RankSumSummed)
            };

            return svmFeatures;
        }

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


        internal ChargeLcScanCluster(int minCharge, int[] ms1ScanNums, IsotopeList isoList)
        {
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
            //RepresentativeMemberIndex = new HashSet<int>();
        }

        internal void SetScore(int scoreType, double scoreValue)
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

        internal void SetRepresentativeMass(double mass, double mz, int charge, int scanNum)
        {
            RepresentativeMass = mass;
            RepresentativeMz = mz;
            RepresentativeScanNum = scanNum;
            RepresentativeCharge = charge;
        }
    }
    

    class ChargeLcScanCell : Tuple<int, int>
    {
        internal ChargeLcScanCell(int row, int col) : base(row, col) { }
        internal int Row { get { return Item1; } }
        internal int Col { get { return Item2; } }
    }
}
