using System;
using System.Collections.Generic;
using System.Linq;
using System.Linq.Expressions;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsPeakCluster : LcMsFeature
    {
        /*
        public LcMsPeakCluster(LcMsRun run, double monoMass, int charge, double mz, double scanNum)
            : base(monoMass, charge, mz, scanNum, 0)
        {
            Run = run;
            TheoreticalEnvelope = observedEnvelope.TheoreticalEnvelope;
            Flag = 0;

            DetectableMaxCharge = (int)Math.Min(Math.Floor(Mass / Run.MinMs1Mz), LcMsPeakMatrix.MaxScanCharge);
            DetectableMinCharge = (int)Math.Max(Math.Ceiling(Mass / Run.MaxMs1Mz), LcMsPeakMatrix.MinScanCharge);

            RepresentativeSummedEnvelop = new double[TheoreticalEnvelope.Size];

            AbundanceDistributionAcrossCharge = new double[2];
            BestCorrelationScoreAcrossCharge = new double[2];
            BestDistanceScoreAcrossCharge = new double[2];
            BestIntensityScoreAcrossCharge = new double[2];

            EnvelopeDistanceScoreAcrossCharge = new double[2];
            EnvelopeCorrelationScoreAcrossCharge = new double[2];
            EnvelopeIntensityScoreAcrossCharge = new double[2];
            BestCharge = new int[2];

            XicCorrelationBetweenBestCharges = new double[2];
            _initScore = false;
        }*/
        
        
        public LcMsPeakCluster(LcMsRun run, ObservedIsotopeEnvelope observedEnvelope)
            : base(observedEnvelope.MonoMass, observedEnvelope.Charge, observedEnvelope.RepresentativePeak.Mz, observedEnvelope.ScanNum, observedEnvelope.Abundance)
        {
            Run = run;
            TheoreticalEnvelope = observedEnvelope.TheoreticalEnvelope;
            Flag = 0;

            DetectableMaxCharge = (int)Math.Min(Math.Floor(Mass / Run.MinMs1Mz), LcMsPeakMatrix.MaxScanCharge);
            DetectableMinCharge = (int)Math.Max(Math.Ceiling(Mass / Run.MaxMs1Mz), LcMsPeakMatrix.MinScanCharge);

            RepresentativeSummedEnvelop = new double[TheoreticalEnvelope.Size];

            AbundanceDistributionAcrossCharge = new double[2];
            BestCorrelationScoreAcrossCharge = new double[2];
            BestDistanceScoreAcrossCharge = new double[2];
            BestIntensityScoreAcrossCharge = new double[2];

            EnvelopeDistanceScoreAcrossCharge = new double[2];
            EnvelopeCorrelationScoreAcrossCharge = new double[2];
            EnvelopeIntensityScoreAcrossCharge = new double[2];
            BestCharge = new int[2];

            XicCorrelationBetweenBestCharges = new double[2];
            _initScore = false;
        }

        public void AddEnvelopes(int minCharge, int maxCharge, int minScanNum, int maxScanNum,
            IList<ObservedIsotopeEnvelope> envelopes)
        {
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            var minCol = ms1ScanNumToIndex[minScanNum];
            var maxCol = ms1ScanNumToIndex[maxScanNum];

            var nRows = maxCharge - minCharge + 1;
            var nCols = maxCol - minCol + 1;

            MinCharge = minCharge;
            MaxCharge = maxCharge;
            MinScanNum = minScanNum;
            MaxScanNum = maxScanNum;

            Envelopes = new ObservedIsotopeEnvelope[nRows][];
            for(var i = 0; i < nRows; i++) Envelopes[i] = new ObservedIsotopeEnvelope[nCols];
            
            foreach (var envelope in envelopes)
            {
                var i = envelope.Charge - MinCharge;
                var j = ms1ScanNumToIndex[envelope.ScanNum] - minCol;

                if (i < 0 || i >= nRows || j < 0 || j >= nCols) continue;
                Envelopes[i][j] = envelope;
            }
        }

        public void UpdateWithDecoyScore(List<Ms1Spectrum> ms1Spectra)
        {
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            var ms1ScanNums = Run.GetMs1ScanVector();
            var minCol = ms1ScanNumToIndex[MinScanNum];
            var maxCol = ms1ScanNumToIndex[MaxScanNum];
            MinCharge = DetectableMinCharge;
            MaxCharge = DetectableMaxCharge;
            
            var rnd = new Random();
            var comparer = new MzComparerWithBinning(28);
            var mostAbuInternalIndex = TheoreticalEnvelope.IndexOrderByRanking[0];

            var nRows = MaxCharge - MinCharge + 1;
            var nCols = maxCol - minCol + 1;

            Envelopes = new ObservedIsotopeEnvelope[nRows][];
            for (var i = 0; i < nRows; i++) Envelopes[i] = new ObservedIsotopeEnvelope[nCols];
            
            for (var charge = DetectableMinCharge; charge <= DetectableMaxCharge; charge++)
            {
                var mostAbuMz = TheoreticalEnvelope.GetIsotopeMz(charge, mostAbuInternalIndex);
                if (Run.MaxMs1Mz < mostAbuMz || mostAbuMz < Run.MinMs1Mz) continue;

                for (var col = minCol; col <= maxCol; col++)
                {
                    var localWin = ms1Spectra[col].GetLocalMzWindow(mostAbuMz);

                    var numMzBins = comparer.GetBinNumber(localWin.MaxMz) - comparer.GetBinNumber(localWin.MinMz) + 1;
                    var peakSet = new Ms1Peak[TheoreticalEnvelope.Size];

                    for (var k = 0; k < peakSet.Length; k++)
                    {
                        var r = rnd.Next(0, numMzBins);
                        if (r < localWin.PeakCount)
                            peakSet[k] = (Ms1Peak) ms1Spectra[col].Peaks[r + localWin.PeakStartIndex];
                    }

                    var env = new ObservedIsotopeEnvelope(Mass, charge, ms1ScanNums[col], peakSet, TheoreticalEnvelope);
                    //AddObservedEnvelope(env);
                    Envelopes[charge - MinCharge][col - minCol] = env;
                }
            }
            UpdateScore(ms1Spectra, false);
        }

        private void ClearScore()
        {
            Array.Clear(AbundanceDistributionAcrossCharge, 0, 2);

            Array.Clear(BestCorrelationScoreAcrossCharge, 0, 2);
            Array.Clear(BestIntensityScoreAcrossCharge, 0, 2);
            BestDistanceScoreAcrossCharge[0] = 1;
            BestDistanceScoreAcrossCharge[1] = 1;

            Array.Clear(EnvelopeCorrelationScoreAcrossCharge, 0, 2);
            Array.Clear(EnvelopeIntensityScoreAcrossCharge, 0, 2);
            EnvelopeDistanceScoreAcrossCharge[0] = 1;
            EnvelopeDistanceScoreAcrossCharge[1] = 1;

            Array.Clear(BestCharge, 0, 2);

            Array.Clear(XicCorrelationBetweenBestCharges, 0, 2);
        }

        public void UpdateScore(List<Ms1Spectrum> ms1Spectra, bool pValueCheck = true)
        {
            var nRows = MaxCharge - MinCharge + 1;
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            var minCol = ms1ScanNumToIndex[MinScanNum];
            var maxCol = ms1ScanNumToIndex[MaxScanNum];
            var nCols = maxCol - minCol + 1;
            var mostAbuIdx = TheoreticalEnvelope.IndexOrderByRanking[0];
            
            ClearScore();

            var bestChargeDist = new double[]{10.0d, 10.0d};
            // sum envelopes at each charge 
            var summedIntensity = new double[TheoreticalEnvelope.Size];

            var xicLen = nCols;
            var xicStartIdx = 0;
            if (nCols < 13)
            {
                xicLen = 13;
                xicStartIdx = (int) Math.Floor((xicLen - nCols)*0.5);
            }
            
            var xic2 = new double[2][];
            xic2[0] = new double[xicLen];
            xic2[1] = new double[xicLen];
            var chargeXic = new double[nRows][];

            var tempBestBcDist = 10.0d;
            var repEnvelopeBcDist = 10.0d;
            ObservedIsotopeEnvelope repEnvelope = null;

            var repEnvelopeBcDist2 = 10.0d;
            ObservedIsotopeEnvelope repEnvelope2 = null;


            var tempBestDistanceScoreAcrossCharge = new double[2]{10, 10};
            var tempBestIntensityScoreAcrossCharge = new double[2];
            var tempBestCorrelationScoreAcrossCharge = new double[2];

            for (var i = 0; i < nRows; i++)
            {
                var charge = i + MinCharge;
                var mostAbuMz = TheoreticalEnvelope.GetIsotopeMz(charge, mostAbuIdx);
                Array.Clear(summedIntensity, 0, summedIntensity.Length);

                chargeXic[i] = new double[xicLen];

                var chargeIdx = (charge % 2 == 0) ? EvenCharge : OddCharge;
                var summedMostAbuIsotopeIntensity = 0d;
                var summedReferenceIntensity = 0d;

                for (var j = 0; j < nCols; j++)
                {
                    var envelope = Envelopes[i][j];
                    var col = minCol + j;
                    
                    var localWin = ms1Spectra[col].GetLocalMzWindow(mostAbuMz);

                    if (envelope == null) continue;
                    
                    envelope.Peaks.SumEnvelopeTo(summedIntensity);
                    var mostAbuPeak = envelope.Peaks[mostAbuIdx];
                        
                    if (mostAbuPeak != null && mostAbuPeak.Active)
                    {
                        summedMostAbuIsotopeIntensity += mostAbuPeak.Intensity;
                        summedReferenceIntensity += localWin.HighestIntensity;
                    }
                    AbundanceDistributionAcrossCharge[chargeIdx] += envelope.Abundance;

                    var newBcDist = TheoreticalEnvelope.GetBhattacharyyaDistance(envelope.Peaks);
                    var newCorr = TheoreticalEnvelope.GetPearsonCorrelation(envelope.Peaks);

                    if (newBcDist < 0.2 || newCorr > 0.6) xic2[chargeIdx][xicStartIdx+j] += envelope.Abundance;
                    if (newBcDist < 0.2 || newCorr > 0.6) chargeXic[i][xicStartIdx + j] = envelope.Abundance;

                    var levelOneEnvelope = true;
                    var levelTwoEnvelope = true;

                    if (pValueCheck)
                    {
                        var poissonPvalue = localWin.GetPoissonTestPvalue(envelope.Peaks, TheoreticalEnvelope.Size);
                        var rankSumPvalue = localWin.GetRankSumTestPvalue(envelope.Peaks, TheoreticalEnvelope.Size);
                        levelOneEnvelope = (rankSumPvalue < 0.01 && poissonPvalue < 0.01);
                        levelTwoEnvelope = (rankSumPvalue < 0.01 || poissonPvalue < 0.01);
                    }

                    if (levelOneEnvelope)
                    {
                        if (newBcDist < BestDistanceScoreAcrossCharge[chargeIdx])
                        {
                            BestDistanceScoreAcrossCharge[chargeIdx] = newBcDist;
                            if (localWin.MedianIntensity > 0)
                                BestIntensityScoreAcrossCharge[chargeIdx] = envelope.HighestIntensity / localWin.HighestIntensity;
                            else BestIntensityScoreAcrossCharge[chargeIdx] = 1.0d;
                        }
                        
                        BestCorrelationScoreAcrossCharge[chargeIdx] = Math.Max(BestCorrelationScoreAcrossCharge[chargeIdx], newCorr);

                        if (newBcDist < repEnvelopeBcDist)
                        {
                            repEnvelopeBcDist = newBcDist;
                            repEnvelope = envelope;
                        }

                        // in the initial scoring, classify major and minor envelopes
                        if (!_initScore)
                        {
                            if (newBcDist < 0.07 || newCorr > 0.7) envelope.GoodEnough = true;
                        }
                    }

                    if (levelTwoEnvelope)
                    {
                        // without considering p-values
                        if (newBcDist < tempBestDistanceScoreAcrossCharge[chargeIdx])
                        {
                            tempBestDistanceScoreAcrossCharge[chargeIdx] = newBcDist;
                            if (localWin.MedianIntensity > 0)
                                tempBestIntensityScoreAcrossCharge[chargeIdx] = envelope.HighestIntensity / localWin.HighestIntensity;
                            else tempBestIntensityScoreAcrossCharge[chargeIdx] = 1.0d;
                        }
                        tempBestCorrelationScoreAcrossCharge[chargeIdx] = Math.Max(tempBestCorrelationScoreAcrossCharge[chargeIdx], newCorr);

                        if (newBcDist < repEnvelopeBcDist2)
                        {
                            repEnvelopeBcDist2 = newBcDist;
                            repEnvelope2 = envelope;
                        }
                    }
                }

                var bcDist = TheoreticalEnvelope.GetBhattacharyyaDistance(summedIntensity);
                EnvelopeDistanceScoreAcrossCharge[chargeIdx] = Math.Min(bcDist, EnvelopeDistanceScoreAcrossCharge[chargeIdx]);
                EnvelopeCorrelationScoreAcrossCharge[chargeIdx] = Math.Max(TheoreticalEnvelope.GetPearsonCorrelation(summedIntensity), EnvelopeCorrelationScoreAcrossCharge[chargeIdx]);

                if (BestCharge[chargeIdx] < 1 || bcDist < bestChargeDist[chargeIdx])
                {
                    BestCharge[chargeIdx] = charge;
                    bestChargeDist[chargeIdx] = bcDist;
                    if (summedReferenceIntensity > 0)
                        EnvelopeIntensityScoreAcrossCharge[chargeIdx] = summedMostAbuIsotopeIntensity/summedReferenceIntensity;
                    //if (summedMedianIntensity > 0) EnvelopeIntensityScoreAcrossCharge[chargeIdx] = Math.Min(1.0, 0.1*(summedMostAbuIsotopeIntensity / summedMedianIntensity));
                }

                if (bcDist < tempBestBcDist)
                {
                    tempBestBcDist = bcDist;
                    Array.Copy(summedIntensity, RepresentativeSummedEnvelop, RepresentativeSummedEnvelop.Length);
                }
            }

            // when good envellope is observed at only either even or odd charge...
            if (BestCorrelationScoreAcrossCharge[0] > 0.7 && BestCorrelationScoreAcrossCharge[1] < 0.1)
            {
                var i = 1;    
                BestCorrelationScoreAcrossCharge[i] = tempBestCorrelationScoreAcrossCharge[i];
                BestIntensityScoreAcrossCharge[i] = tempBestIntensityScoreAcrossCharge[i];
                BestDistanceScoreAcrossCharge[i] = tempBestDistanceScoreAcrossCharge[i];
            }
            if (BestCorrelationScoreAcrossCharge[1] > 0.7 && BestCorrelationScoreAcrossCharge[0] < 0.1)
            {
                var i = 0;
                BestCorrelationScoreAcrossCharge[i] = tempBestCorrelationScoreAcrossCharge[i];
                BestIntensityScoreAcrossCharge[i] = tempBestIntensityScoreAcrossCharge[i];
                BestDistanceScoreAcrossCharge[i] = tempBestDistanceScoreAcrossCharge[i];
            }

            
            // normalize abudnace across charges
            var s = AbundanceDistributionAcrossCharge[0] + AbundanceDistributionAcrossCharge[1];
            if (s > 0)
            {
                for (var chargeIdx = 0; chargeIdx < 2; chargeIdx++) AbundanceDistributionAcrossCharge[chargeIdx] = AbundanceDistributionAcrossCharge[chargeIdx] / s;    
            }

            if (nCols > 1)
            {
                var evenChargeIdx = BestCharge[EvenCharge] - MinCharge;
                var oddChargeIdx = BestCharge[OddCharge] - MinCharge;
                XicCorrelationBetweenBestCharges[0] = FitScoreCalculator.GetPearsonCorrelation(Smoother.Smooth(chargeXic[evenChargeIdx]), Smoother.Smooth(chargeXic[oddChargeIdx]));
                XicCorrelationBetweenBestCharges[1] = FitScoreCalculator.GetPearsonCorrelation(Smoother.Smooth(xic2[EvenCharge]), Smoother.Smooth(xic2[OddCharge]));
            }

            if (repEnvelope == null && repEnvelope2 != null) repEnvelope = repEnvelope2;

            if (repEnvelope != null)
            {
                // set representative charge, mz and scanNum
                RepresentativeCharge = repEnvelope.Charge;
                RepresentativeMz = repEnvelope.RepresentativePeak.Mz;
                RepresentativeScanNum = repEnvelope.ScanNum;                
            }

            _initScore = true;
        }

        public void SetChargeRange(int minCharge, int maxCharge)
        {
            MinCharge = minCharge;
            MaxCharge = MaxCharge;
        }

        public void SetAbundance(double abu)
        {
            Abundance = abu;
        }

        internal void Expand(ObservedIsotopeEnvelope envelope)
        {
            if (MaxScanNum < 0 || envelope.ScanNum > MaxScanNum)
            {
                MaxScanNum = envelope.ScanNum;
            }
            if (MinScanNum < 0 || envelope.ScanNum < MinScanNum)
            {
                MinScanNum = envelope.ScanNum;
            }
            if (MaxCharge < 0 || envelope.Charge > MaxCharge)
            {
                MaxCharge = envelope.Charge;
            }
            if (MinCharge < 0 || envelope.Charge < MinCharge)
            {
                MinCharge = envelope.Charge;
            }            
        }
        
        public double CoElutionLength(LcMsPeakCluster other)
        {
            if (other.MaxScanNum >= MinScanNum && other.MaxScanNum <= MaxScanNum)
            {
                return (other.MaxElutionTime - Math.Max(MinElutionTime, other.MinElutionTime));
            }

            if (MaxScanNum >= other.MinScanNum && MaxScanNum <= other.MaxScanNum)
            {
                return (MaxElutionTime - Math.Max(other.MinElutionTime, MinElutionTime));
            }
            
            return 0d;
        }
        
        public void ExpandElutionRange()
        {
            // considering DDA instrument
            if (Run.MaxMsLevel > 1 && NetLength < 0.01)
            {
                var ms1ScanNums = Run.GetMs1ScanVector();
                var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

                var minCol = ms1ScanNumToIndex[MinScanNum];
                var maxCol = ms1ScanNumToIndex[MaxScanNum];

                for (var i = minCol - 1; i >= 0; i--)
                {
                    if (MinScanNum - ms1ScanNums[i] > (minCol - i))
                    {
                        MinScanNum = ms1ScanNums[i];
                        break;
                    }
                }
                for (var i = maxCol + 1; i < ms1ScanNums.Length; i++)
                {
                    if (ms1ScanNums[i] - MaxScanNum > (i - maxCol))
                    {
                        MaxScanNum = ms1ScanNums[i];
                        break;
                    }
                }
            }
        }
        
        public IEnumerable<ObservedIsotopeEnvelope> EnumerateEnvelopes()
        {
            for (var i = 0; i < Envelopes.Length; i++)
            {
                for (var j = 0; j < Envelopes[0].Length; j++)
                {
                    var envelope = Envelopes[i][j];
                    if (envelope != null) yield return envelope;
                }
            }
        }

        public IEnumerable<Ms1Peak> GetMajorPeaks()
        {
            foreach (var envelope in EnumerateEnvelopes())
            {
                if (!envelope.GoodEnough) continue;

                for (var i = 0; i < TheoreticalEnvelope.Size; i++)
                {
                    if (TheoreticalEnvelope.Isotopes[i].Ratio > 0.3)
                    {
                        var peak = envelope.Peaks[i];
                        if (peak != null && peak.Active) yield return peak;                        
                    }
                }
            }
        }

        public IEnumerable<Ms1Peak> GetMinorPeaks()
        {
            foreach (var envelope in EnumerateEnvelopes())
            {
                if (envelope.GoodEnough) continue;

                for (var i = 0; i < TheoreticalEnvelope.Size; i++)
                {
                    if (TheoreticalEnvelope.Isotopes[i].Ratio <= 0.3)
                    {
                        var peak = envelope.Peaks[i];
                        if (peak != null && peak.Active) yield return peak;
                    }
                }
            }            
        }

        public void InActivateMajorPeaks()
        {
            foreach (var peak in GetMajorPeaks()) peak.InActivate();
        }

        public HashSet<LcMsPeakCluster> OverlappedFeatures
        {
            // definition
            // Ovelapped features = features whose scores are affected by inactivating major peaks of this feature;
            //                    = features that share any major peaks in this features as their peaks
            //                     + features that share any minor peaks in this  feature as their major peaks
            get
            {
                var ret = new HashSet<LcMsPeakCluster>();
                foreach (var peak in GetMajorPeaks())
                {
                    foreach (var f in peak.GetAllTaggedFeatures()) ret.Add(f);
                }
                foreach (var peak in GetMinorPeaks())
                {
                    foreach (var f in peak.GetMajorTaggedFeatures()) ret.Add(f);
                }

                return ret;
            }
        }

        public bool GoodEnougth
        {
            get
            {
                if (Mass > 35000) return BestCorrelationScore > 0.85;
                if (Mass > 25000) return BestCorrelationScore > 0.8;
                if (Mass > 15000) return BestCorrelationScore > 0.75;
                return BestCorrelationScore > 0.7;                
            }
        }

        public readonly int DetectableMaxCharge;
        public readonly int DetectableMinCharge; 
        public ObservedIsotopeEnvelope[][] Envelopes;
        public double Score { get; internal set; }

        public readonly int[] BestCharge;

        public readonly double[] RepresentativeSummedEnvelop;
        
        public readonly double[] EnvelopeDistanceScoreAcrossCharge;
        public readonly double[] EnvelopeCorrelationScoreAcrossCharge;
        public readonly double[] EnvelopeIntensityScoreAcrossCharge;
        
        public readonly double[] AbundanceDistributionAcrossCharge;

        public readonly double[] BestCorrelationScoreAcrossCharge;
        public readonly double[] BestDistanceScoreAcrossCharge;
        public readonly double[] BestIntensityScoreAcrossCharge;

        public readonly double[] XicCorrelationBetweenBestCharges;

        public const int EvenCharge = 0;
        public const int OddCharge = 1;

        public double BestCorrelationScore { get { return Math.Max(BestCorrelationScoreAcrossCharge.Max(), EnvelopeCorrelationScoreAcrossCharge.Max());  } }
        public readonly TheoreticalIsotopeEnvelope TheoreticalEnvelope;
        private static readonly SavitzkyGolaySmoother Smoother = new SavitzkyGolaySmoother(9, 2);
        internal byte Flag;

        private bool _initScore;
    }
}
