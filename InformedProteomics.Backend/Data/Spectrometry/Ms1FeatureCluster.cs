using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public static class Ms1FeatureScore
    {
        public const byte EnvelopeCorrelation = 0;
        public const byte EnvelopeCorrelationSummed = 1;
        public const byte RankSum = 2;
        public const byte Poisson = 3;

        public const byte BhattacharyyaDistance = 4;
        public const byte BhattacharyyaDistanceSummed = 5;
        public const byte BhattacharyyaDistanceSummedOverCharges = 6;
        public const byte BhattacharyyaDistanceSummedOverTimes = 7;

        public const byte XicCorrMean = 8;
        public const byte XicCorrMin = 9;
        public const byte MzError = 10;
        public const byte TotalMzError = 11;
        public const byte Count = 12;
    }

    public class Ms1FeatureCluster : Ms1Feature
    {
        public IsotopeList IsotopeList { get; private set; }
        //public double ClusteringScore { get; internal set; }
        //public double ClusteringScore2 { get; internal set; }
        public double Probability { get; private set; }

        public override int MinScanNum
        {
            get
            {
                if ((ScanLength < 11) && MinCol - 1 >= 0) return _ms1ScanNums[MinCol - 1];
                return _ms1ScanNums[MinCol];
            }
        }

        public override int MaxScanNum
        {
            get
            {
                if ((ScanLength < 11) && MaxCol + 1 < _ms1ScanNums.Length) return _ms1ScanNums[MaxCol + 1];
                return _ms1ScanNums[MaxCol];
            }
        }

        public override int MinCharge { get { return _minCharge + MinRow; } }
        public override int MaxCharge { get { return _minCharge + MaxRow; } }

        public int ScanLength { get { return MaxCol - MinCol + 1; } }
        public int ChargeLength { get { return MaxRow - MinRow + 1; } }

        public int GoodEnvelopeCount { get; internal set; }

        public double Abundance { get { return SummedEnvelope.Sum(); } }
        public double GetScore(byte scoreType) { return _scores[scoreType]; }
        public double[] SummedEnvelope { get; private set; }
        public double NormalizedElutionLength { get; internal set; }

        public IEnumerable<Ms1Peak> GetMajorPeaks()
        {
            foreach (var envelope in Envelopes)
            {
                if (!envelope.GoodEnough) continue;
                for (var i = 0; i < IsotopeList.Count; i++)
                {
                    var j = IsotopeList.SortedIndexByIntensity[i];
                    var peak = envelope.Peaks[j];
                    if (peak == null) continue;
                    yield return peak;
                }
            }
        }

        public IEnumerable<Ms1Peak> GetMinorPeaks()
        {
            foreach (var envelope in Envelopes)
            {
                if (envelope.GoodEnough) continue;

                for (var i = IsotopeList.Count - 1; i >= 0; i--)
                {
                    var j = IsotopeList.SortedIndexByIntensity[i];
                    var peak = envelope.Peaks[j];
                    if (peak == null) continue;
                    
                    yield return peak;
                }
            }
        }

        public void InActivateSignificantPeaks()
        {
            foreach (var peak in GetMajorPeaks()) peak.InActivate();
        }

        public static readonly string[] TsvHeaderWithScore = new string[]
        {
            "FeatureID", "MinScan", "MaxScan", "MinCharge", "MaxCharge", "MonoMass", "RepScan", "RepCharge", "RepMz", "Abundance",
            "ObservedEnvelopes", "ScanLength",
            "BestCorr", "SummedCorr", 
            "RanksumScore", "PoissonScore",
            "BcDist", "SummedBcDist", 
            "SummedBcDistOverCharges", "SummedBcDistOverTimes", 
            "XicDist", "MinXicDist",
            "MzDiffPpm", "TotalMzDiffPpm",
            "Envelope", 
            "Probability", "GoodEnough"
        };
       
        private static readonly string[] TsvHeader = new string[]
        {
            "FeatureID", "MinScan", "MaxScan", "MinCharge", "MaxCharge", "MonoMass", "RepScan", "RepCharge", "RepMz", "Abundance",
            "BestCorr", "SummedCorr", 
            "XicDist", 
            "MzDiffPpm", "Envelope",
            "Probability", "GoodEnough"
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
            var intensity = SummedEnvelope;

            var sb = new StringBuilder(string.Format("{0}\t{1}\t{2}\t{3}\t{4:0.0000}\t{5}\t{6}\t{7:0.0000}\t{8:0.00}",
                                        MinScanNum, MaxScanNum, MinCharge, MaxCharge, RepresentativeMass, RepresentativeScanNum, RepresentativeCharge,
                                        RepresentativeMz, intensity.Sum()));
            
            if (withScore)
            {
                sb.AppendFormat("\t{0}", GoodEnvelopeCount);
                sb.AppendFormat("\t{0}", ScanLength);
                foreach (var score in _scores) sb.AppendFormat("\t{0}", score);
            }
            else
            {
                sb.AppendFormat("\t{0:0.00000}", GetScore(Ms1FeatureScore.EnvelopeCorrelation));
                sb.AppendFormat("\t{0:0.00000}", GetScore(Ms1FeatureScore.EnvelopeCorrelationSummed));
                sb.AppendFormat("\t{0:0.00000}", GetScore(Ms1FeatureScore.MzError));
                sb.AppendFormat("\t{0:0.00000}", GetScore(Ms1FeatureScore.XicCorrMean));
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

        internal List<ObservedEnvelope> Envelopes { get; private set; }
        
        internal double[] ClusteringEnvelope { get; private set; }
        internal int MaxCol { get; private set; }
        internal int MinCol { get; private set; }
        internal int MaxRow { get; private set; }
        internal int MinRow { get; private set; }
        //internal readonly HashSet<Ms1FeatureCluster> OverlappedFeatures;

        internal HashSet<Ms1FeatureCluster> OverlappedFeatures
        {
            // definition
            // Ovelapped features = features whose scores are affected by inactivating major peaks of this feature;
            //                    = features that share any major peaks in this features as their peaks
            //                     + features that share any minor peaks in this  feature as their major peaks
            get
            {
                var ret = new HashSet<Ms1FeatureCluster>();
                foreach (var peak in GetMajorPeaks())
                {
                    foreach (var f in peak.GetTaggedAllFeatures()) ret.Add(f);
                }
                foreach (var peak in GetMinorPeaks())
                {
                    foreach (var f in peak.GetTaggedMajorFeatures()) ret.Add(f);
                }

                return ret;
            }
        }

        private readonly double[] _scores;
        private readonly int[] _ms1ScanNums;
        private readonly int _minCharge;
        private static readonly double LogP1 = -Math.Log(0.01, 2);
        private static readonly double LogP2 = -Math.Log(0.02, 2);
        private static readonly double LogP5 = -Math.Log(0.05, 2);

        internal Ms1FeatureCluster(int minCharge, int[] ms1ScanNums, IsotopeList isoList)
        {
            _ms1ScanNums = ms1ScanNums;
            _minCharge = minCharge;

            IsotopeList = isoList;
            Envelopes = new List<ObservedEnvelope>();
            
            ClusteringEnvelope = new double[IsotopeList.Count];
            //ClusteringScore = 0;
            //ClusteringScore2 = 0; 
            _scores = new double[Ms1FeatureScore.Count];
            
            MaxCol = 0;
            MinCol = _ms1ScanNums.Length - 1;
            MaxRow = 0;
            MinRow = 99999;
            SummedEnvelope = new double[isoList.Count];

            Flag = 0;
        }

        internal void UpdateScores(bool firstTime = false)
        {
            if (Envelopes.Count < 1) return;

            var bestEnvelopeCorrelation = 0.0d;
            var bestBhattacharyyaDistance = 1.0d;
            var bestRankSumScore = 1.0d;
            var bestPoissonScore = 1.0d;

            ObservedEnvelope repEnvelope = null;
            ObservedEnvelope repEnvelope2 = null;
            var repEnvelope2Bc = 10.0d;

            var envelopePerTime = new double[ScanLength][];
            var envelopePerCharge = new double[ChargeLength][];
            for (var i = 0; i < ScanLength; i++) envelopePerTime[i] = new double[IsotopeList.Count];
            for (var i = 0; i < ChargeLength; i++) envelopePerCharge[i] = new double[IsotopeList.Count];
            GoodEnvelopeCount = 0;

            var memberCorr = new double[Envelopes.Count];
            var memberBc = new double[Envelopes.Count];

            for (var e = 0; e < Envelopes.Count; e++)
            {
                var envelope = Envelopes[e];
                var mostAbuPeak = envelope.Peaks[IsotopeList.SortedIndexByIntensity[0]];
                if (!mostAbuPeak.Active)
                {
                    memberCorr[e] = 0d;
                    memberBc[e] = 1.0d;
                    continue;
                }

                var spectrum = mostAbuPeak.Spectrum;
                var statSigTestResult = spectrum.TestStatisticalSignificance(IsotopeList, envelope.Peaks);
                var envCorr = envelope.GetPearsonCorrelation(IsotopeList.Envelope);
                var bcDistance = envelope.GetBhattacharyyaDistance(IsotopeList.EnvelopePdf);

                if (bcDistance < repEnvelope2Bc)
                {
                    repEnvelope2Bc = bcDistance;
                    repEnvelope2 = envelope;
                }

                for (var i = 0; i < IsotopeList.Count; i++)
                {
                    if (envelope.Peaks[i] == null || !envelope.Peaks[i].Active) continue;
                    envelopePerTime[envelope.Col - MinCol][i] += envelope.Peaks[i].Intensity;
                    envelopePerCharge[envelope.Row - MinRow][i] += envelope.Peaks[i].Intensity;
                }

                if (envCorr > 0.6 && bcDistance < 0.3 && statSigTestResult.PoissonScore > LogP2 && statSigTestResult.RankSumScore > LogP2)
                {
                    GoodEnvelopeCount++;
                    if (firstTime) envelope.GoodEnough = true;
                    
                    if (bestPoissonScore < statSigTestResult.PoissonScore) bestPoissonScore = statSigTestResult.PoissonScore;
                    if (bestRankSumScore < statSigTestResult.RankSumScore) bestRankSumScore = statSigTestResult.RankSumScore;

                    if (bestEnvelopeCorrelation < envCorr) bestEnvelopeCorrelation = envCorr;
                    if (bcDistance < bestBhattacharyyaDistance)
                    {
                        bestBhattacharyyaDistance = bcDistance;
                        repEnvelope = envelope;
                    }
                }

                memberCorr[e] = envCorr;
                memberBc[e] = bcDistance;
            }

            if (RepresentativeCharge < 1)
            {
                if (repEnvelope == null) repEnvelope = repEnvelope2;
                if (repEnvelope == null) throw new Exception("Cannot find representative envelope");

                var mostAbuInternalIndex = IsotopeList.SortedIndexByIntensity[0];
                var repPeak = repEnvelope.Peaks[mostAbuInternalIndex];
                RepresentativeCharge = _minCharge + repEnvelope.Row;
                RepresentativeMass = Ion.GetMonoIsotopicMass(repPeak.Mz, RepresentativeCharge, IsotopeList[mostAbuInternalIndex].Index);
                RepresentativeMz = repPeak.Mz;
                RepresentativeScanNum = _ms1ScanNums[repEnvelope.Col];
            }
            
            var bestBcPerTime = 10d;
            for (var col = MinCol; col <= MaxCol; col++)
            {
                var bc = IsotopeList.GetBhattacharyyaDistance(envelopePerTime[col - MinCol]);
                if (bc < bestBcPerTime) bestBcPerTime = bc;
            }
            var bestBcPerCharge = 10d;
            for (var row = MinRow; row <= MaxRow; row++)
            {
                var bc = IsotopeList.GetBhattacharyyaDistance(envelopePerCharge[row - MinRow]);
                if (bc < bestBcPerCharge) bestBcPerCharge = bc;
            }

            var bestMzErrorPpm = 10d;
            var totalMzErrorPpm = 0d;
            var totalMzPairCount = 0;
            foreach (var envelope in Envelopes.Where(e => e.GoodEnough))
            {
                var mzErrorPpm = 0d;
                var n = 0;
                var charge = _minCharge + envelope.Row;
                for (var j = 0; j < IsotopeList.Count; j++)
                {
                    if (envelope.Peaks[j] == null || !envelope.Peaks[j].Active) continue;
                    var theoreticalMz = Ion.GetIsotopeMz(RepresentativeMass, charge, IsotopeList[j].Index);
                    mzErrorPpm += (Math.Abs(envelope.Peaks[j].Mz - theoreticalMz) * 1e6) / theoreticalMz;
                    n++;
                }

                totalMzErrorPpm += mzErrorPpm;
                totalMzPairCount += n;
                mzErrorPpm /= n;
                if (mzErrorPpm < bestMzErrorPpm) bestMzErrorPpm = mzErrorPpm;                
            }
            
            SetScore(Ms1FeatureScore.TotalMzError, totalMzErrorPpm/totalMzPairCount);
            SetScore(Ms1FeatureScore.EnvelopeCorrelation, bestEnvelopeCorrelation);
            SetScore(Ms1FeatureScore.RankSum, bestRankSumScore);
            SetScore(Ms1FeatureScore.Poisson, bestPoissonScore);
            SetScore(Ms1FeatureScore.BhattacharyyaDistance, bestBhattacharyyaDistance);
            SetScore(Ms1FeatureScore.BhattacharyyaDistanceSummedOverCharges, bestBcPerTime);
            SetScore(Ms1FeatureScore.BhattacharyyaDistanceSummedOverTimes, bestBcPerCharge);
            SetScore(Ms1FeatureScore.MzError, bestMzErrorPpm);
            SetSummedEnvelope(memberCorr, memberBc);
          
            Probability = GetProbabilityByLogisticRegression();
        }

        internal void SetScore(byte scoreType, double scoreValue)
        {
            _scores[scoreType] = scoreValue;
        }

        internal void AddMember(ObservedEnvelope envelope)
        {
            if (envelope.Col > MaxCol) MaxCol = envelope.Col;
            if (envelope.Col < MinCol) MinCol = envelope.Col;
            if (envelope.Row > MaxRow) MaxRow = envelope.Row;
            if (envelope.Row < MinRow) MinRow = envelope.Row;
            Envelopes.Add(envelope);

            //if (updatedClusteringScore < 0) return;
            //envelope.SumEnvelopeTo(ClusteringEnvelope);
            //ClusteringScore = updatedClusteringScore;
        }

        internal void ClearMember()
        {
            MaxCol = -1;
            MinCol = _ms1ScanNums.Length - 1;
            MaxRow = -1;
            MinRow = 99999;
            Envelopes.Clear();
        }
      
        private void SetSummedEnvelope(double[] cellCorr, double[] cellBc)
        {
            if (Envelopes.Count == 1)
            {
                Array.Clear(SummedEnvelope, 0, IsotopeList.Count);
                Envelopes[0].SumEnvelopeTo(SummedEnvelope);
                SetScore(Ms1FeatureScore.BhattacharyyaDistanceSummed, GetScore(Ms1FeatureScore.BhattacharyyaDistance));
                SetScore(Ms1FeatureScore.EnvelopeCorrelationSummed, GetScore(Ms1FeatureScore.EnvelopeCorrelation));
                return;
            }

            Array.Clear(SummedEnvelope, 0, IsotopeList.Count);
            var tempEnvelop = new double[IsotopeList.Count];
            var summedBc = 99999d;
            var index2 = Enumerable.Range(0, Envelopes.Count).ToArray();
            Array.Sort(cellBc, index2);
            foreach (var j in index2)
            {
                Array.Copy(SummedEnvelope, tempEnvelop, tempEnvelop.Length);
                Envelopes[j].SumEnvelopeTo(tempEnvelop);

                var tempBc = IsotopeList.GetBhattacharyyaDistance(tempEnvelop);

                if (tempBc > summedBc) continue;
                summedBc = tempBc;

                //seleted[j] = true;
                Array.Copy(tempEnvelop, SummedEnvelope, tempEnvelop.Length);
            }
            SetScore(Ms1FeatureScore.BhattacharyyaDistanceSummed, summedBc);

            var index = Enumerable.Range(0, Envelopes.Count).ToArray();
            Array.Sort(cellCorr, index);
            Array.Reverse(index);
            var summedCorr = 0d;
            Array.Clear(SummedEnvelope, 0, SummedEnvelope.Length);
            foreach (var j in index)
            {
                Array.Clear(tempEnvelop, 0, tempEnvelop.Length);
                Array.Copy(SummedEnvelope, tempEnvelop, tempEnvelop.Length);
                Envelopes[j].SumEnvelopeTo(tempEnvelop);

                var tempCorr = IsotopeList.GetPearsonCorrelation(tempEnvelop);
                if (tempCorr < summedCorr) continue;

                summedCorr = tempCorr;
                //seleted[j] = true;
                Array.Copy(tempEnvelop, SummedEnvelope, tempEnvelop.Length);
            }
            SetScore(Ms1FeatureScore.EnvelopeCorrelationSummed, summedCorr);
            
        }

        internal bool GoodEnough
        {
            get
            {
                if (GoodEnvelopeCount < 1) return false;
                
                // basic filtering conditions
                if (GetScore(Ms1FeatureScore.RankSum) < LogP2) return false;
                if (GetScore(Ms1FeatureScore.Poisson) < LogP2) return false;
                if (GetScore(Ms1FeatureScore.EnvelopeCorrelation) < 0.6) return false;
                if (GetScore(Ms1FeatureScore.EnvelopeCorrelationSummed) < 0.7) return false;
                if (GetScore(Ms1FeatureScore.BhattacharyyaDistance) > 0.3) return false;
                if (GetScore(Ms1FeatureScore.BhattacharyyaDistanceSummed) > 0.1) return false;
                if (GetScore(Ms1FeatureScore.BhattacharyyaDistanceSummedOverCharges) > 0.12) return false;
                if (GetScore(Ms1FeatureScore.BhattacharyyaDistanceSummedOverTimes) > 0.12) return false;
                if (GetScore(Ms1FeatureScore.XicCorrMean) > 0.3) return false;
                if (GetScore(Ms1FeatureScore.XicCorrMin) > 0.2) return false;
                if (GetScore(Ms1FeatureScore.MzError) > 4) return false;
                if (GetScore(Ms1FeatureScore.TotalMzError) > 5) return false;
                
                /** further considerations: 
                  * short features must have a strong evidence
                  * high charge molecules should have high "summed" scores
                ***/
               
                return true;
            }
        }
        internal byte Flag;

        private static readonly double[] LogisticRegressionBetaVector = new double[]
        {
            -6.32017825981311,
            1.81729164993594,
            0.666260087669721,
            0.0300213257125810,
            0.0351494822314954,
            3.44370112481914,
            32.6628108691280,
            -28.9890069475280,
            -71.9737221672874,
            -1.10305154889920,
            3.73686326096962,
            0.122237630801168,
            0.496631413971031              
        };
        
        internal double GetProbabilityByLogisticRegression()
        {
            var eta = LogisticRegressionBetaVector[0];
            for (var i = 1; i < LogisticRegressionBetaVector.Length; i++)
                eta += _scores[i - 1] * LogisticRegressionBetaVector[i];
                
            var pi = Math.Exp(eta);
            return pi/(pi + 1);
        }
    }


    class ObservedEnvelope
    {
        internal ObservedEnvelope(int row, int col, Ms1Peak[] peaks, int nIsotopes)
        {
            Row = row;
            Col = col;
            Peaks = new Ms1Peak[nIsotopes];
            Array.Copy(peaks, Peaks, nIsotopes);
            GoodEnough = false;

            for (var i = 0; i < Peaks.Length; i++)
            {
                if (Peaks[i] != null && Peaks[i].Active)
                {
                    MinMzPeak = Peaks[i];
                    break;
                }
            }
            for (var i = nIsotopes - 1; i >= 0; i--)
            {
                if (Peaks[i] != null && Peaks[i].Active)
                {
                    MaxMzPeak = Peaks[i];
                    break;
                }
            }
        }


        internal int Row { get; private set; }
        internal int Col { get; private set; }
        internal Ms1Peak[] Peaks { get; private set; }
        internal bool GoodEnough { get; set; }

        internal int NumberOfPeaks
        {
            get { return Peaks.Count(x => x != null && x.Active); }
        }

        internal double Abundance
        {
            get
            {
                return Peaks.Where(p => p != null && p.Active).Sum(p => p.Intensity);
            }
        }

        internal void SumEnvelopeTo(double[] targetEnvelope)
        {
            Peaks.SumEnvelopeTo(targetEnvelope);
        }

        internal double GetPearsonCorrelation(double[] theoreticalEnvelope)
        {
            return Peaks.GetPearsonCorrelation(theoreticalEnvelope);
        }

        internal double GetBhattacharyyaDistance(double[] theoreticalEnvelope)
        {
            return Peaks.GetBhattacharyyaDistance(theoreticalEnvelope);
        }
        
        internal Ms1Peak MinMzPeak { get; private set; }
        internal Ms1Peak MaxMzPeak { get; private set; }

        internal bool Overlap(ObservedEnvelope other)
        {
            if (Col != other.Col) return false;
            if (MinMzPeak.Mz > other.MaxMzPeak.Mz) return false;
            if (MaxMzPeak.Mz < other.MinMzPeak.Mz) return false;
            
            return true;
        }
    }
}
