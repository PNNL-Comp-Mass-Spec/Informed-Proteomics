using System;
using System.Collections.Generic;
using System.Linq;
using System.Linq.Expressions;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsPeakCluster : LcMsFeature
    {
        public LcMsPeakCluster(LcMsRun run, ObservedIsotopeEnvelope observedEnvelope)
            : base(observedEnvelope.MonoMass, observedEnvelope.Charge, observedEnvelope.RepresentativePeak.Mz, observedEnvelope.ScanNum, observedEnvelope.Abundance)
        {
            Envelopes = new List<ObservedIsotopeEnvelope>();
            Run = run;
            AddObservedEnvelope(observedEnvelope);
            _theoreticalIsotopeEnvelope = observedEnvelope.TheoreticalEnvelope;
            Flag = 0;
        }

        public void UpdateWithDecoyScore(List<Ms1Spectrum> ms1Spectra, int minScanCharge = 2, int maxScanCharge = 60)
        {
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();
            var ms1ScanNums = Run.GetMs1ScanVector();
            var minCol = ms1ScanNumToIndex[MinScanNum];
            var maxCol = ms1ScanNumToIndex[MaxScanNum];

            var minCharge = (int)Math.Max(Math.Floor(Mass / Run.MaxMs1Mz), minScanCharge);
            var maxCharge = (int)Math.Min(Math.Ceiling(Mass / Run.MinMs1Mz), maxScanCharge);
            

            //int n = _ms1PeakList.Count;
            var rnd = new Random();
            var comparer = new MzComparerWithBinning(28);
            //var decoyEnvelopes = new List<ObservedIsotopeEnvelope>();
            var mostAbuInternalIndex = _theoreticalIsotopeEnvelope.IndexOrderByRanking[0];
            
            Clear();
            for (var charge = minCharge; charge <= maxCharge; charge++)
            {
                var mostAbuMz = _theoreticalIsotopeEnvelope.GetIsotopeMz(charge, mostAbuInternalIndex);
                if (Run.MaxMs1Mz < mostAbuMz || mostAbuMz < Run.MinMs1Mz) continue;

                for (var col = minCol; col <= maxCol; col++)
                {
                    var localWin = ms1Spectra[col].GetLocalMzWindow(mostAbuMz, Mass);

                    if (localWin == null)
                    {
                        continue;
                    }

                    var numMzBins = comparer.GetBinNumber(localWin.MaxMz) - comparer.GetBinNumber(localWin.MinMz) + 1;
                    var peakSet = new Ms1Peak[_theoreticalIsotopeEnvelope.Size];

                    for (var k = 0; k < peakSet.Length; k++)
                    {
                        var r = rnd.Next(0, numMzBins);
                        if (r < localWin.PeakCount)
                            peakSet[k] = (Ms1Peak) ms1Spectra[col].Peaks[r + localWin.PeakStartIndex];
                        //else peakSet[k] = null;
                    }
                    
                    var env = new ObservedIsotopeEnvelope(Mass, charge, ms1ScanNums[col], peakSet, _theoreticalIsotopeEnvelope);
                    //decoyEnvelopes.Add(env);
                    AddObservedEnvelope(env);
                }
            }
            
            MinScanNum = ms1ScanNums[minCol];
            MaxScanNum = ms1ScanNums[maxCol];

            //MinCharge = minCharge;
            //MaxCharge = maxCharge;
            //Envelopes = decoyEnvelopes;
            UpdateScore(ms1Spectra);
        }

        public void UpdateScore(List<Ms1Spectrum> ms1Spectra)
        {
            var nRows = MaxCharge - MinCharge + 1;
            var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

            var summedPeakIntensity = new double[nRows];
            var summedMedianIntensity = new double[nRows];
            
            var summedIntensity = new double[nRows][]; //, _theoreticalIsotopeEnvelope.Size];
            for(var i = 0; i < nRows; i++) summedIntensity[i] = new double[_theoreticalIsotopeEnvelope.Size];

            var mostAbuIdx = Envelopes[0].TheoreticalEnvelope.IndexOrderByRanking[0];
            
            //RepresentativeSummedEnvelop = new double[_theoreticalIsotopeEnvelope.Size];
            
            AbundanceDistributionAcrossCharge = new double[nRows]; 

            BestCorrelationScoreAcrossCharge = new double[nRows];
            BestDistanceScoreAcrossCharge = new double[nRows];
            BestIntensityScoreAcrossCharge = new double[nRows];

            var bestEnvelopeAcrossCharge = new ObservedIsotopeEnvelope[nRows];
            
            // sum envelopes at each charge 
            foreach (var envelope in Envelopes)
            {
                var charge = envelope.Charge;
                var row = charge - MinCharge;
                var col = ms1ScanNumToIndex[envelope.ScanNum];
                envelope.Peaks.SumEnvelopeTo(summedIntensity[row]);
                var mostAbuPeak = envelope.Peaks[mostAbuIdx];
                var mostAbuMz = _theoreticalIsotopeEnvelope.GetIsotopeMz(charge, mostAbuIdx);
                //var repPeak = envelope.RepresentativePeak;

                var newBcDist = _theoreticalIsotopeEnvelope.GetBhattacharyyaDistance(envelope.Peaks);
                if (!(BestDistanceScoreAcrossCharge[row] > 0) || newBcDist < BestDistanceScoreAcrossCharge[row])
                {
                    BestDistanceScoreAcrossCharge[row] = newBcDist;
                    bestEnvelopeAcrossCharge[row] = envelope;
                }

                var newCorr = _theoreticalIsotopeEnvelope.GetPearsonCorrelation(envelope.Peaks);
                BestCorrelationScoreAcrossCharge[charge - MinCharge] =
                    Math.Max(BestCorrelationScoreAcrossCharge[charge - MinCharge], newCorr);

                var localWin = ms1Spectra[col].GetLocalMzWindow(mostAbuMz, Mass);

                if (mostAbuPeak != null && mostAbuPeak.Active)
                {
                    summedPeakIntensity[row] += mostAbuPeak.Intensity;
                    BestIntensityScoreAcrossCharge[row] = Math.Max(BestIntensityScoreAcrossCharge[row], mostAbuPeak.Intensity / localWin.MedianIntensity);
                }
                
                summedMedianIntensity[row] += localWin.MedianIntensity;
                /*
                if (mostAbuPeak != null && mostAbuPeak.Active)
                {
                    var intensity = mostAbuPeak.Intensity;
                    var refIntensity = ms1Spectra[mostAbuPeak.Ms1SpecIndex].GetLocalMedianIntensity(mostAbuPeak, Mass);
                    summedPeakIntensity[row] += intensity;
                    summedMedianIntensity[row] += refIntensity;

                    
                }
                else
                {
                    var repPeak = envelope.RepresentativePeak;
                    if (repPeak != null)
                    {
                        var refIntensity = ms1Spectra[repPeak.Ms1SpecIndex].GetLocalMedianIntensity(repPeak, Mass);
                        summedMedianIntensity[row] += refIntensity;
                    }
                }
                */
                AbundanceDistributionAcrossCharge[row] += envelope.Abundance;
            }
            
            EnvelopeDistanceScoreAcrossCharge = new double[nRows];
            EnvelopeCorrelationScoreAcrossCharge = new double[nRows];
            EnvelopeIntensityScoreAcrossCharge  = new double[nRows];
            
            var tempBestBcDist = 10.0d;
            var tempBestBcDistCharge = 0;
            for (var row = 0; row < nRows; row++)
            {
                var bcDist = _theoreticalIsotopeEnvelope.GetBhattacharyyaDistance(summedIntensity[row]);
                if (!(EnvelopeDistanceScoreAcrossCharge[row] > 0) || bcDist < EnvelopeDistanceScoreAcrossCharge[row])
                    EnvelopeDistanceScoreAcrossCharge[row] = bcDist;

                var corr = _theoreticalIsotopeEnvelope.GetPearsonCorrelation(summedIntensity[row]);
                EnvelopeCorrelationScoreAcrossCharge[row] = Math.Max(EnvelopeCorrelationScoreAcrossCharge[row], corr);

                if (bcDist < tempBestBcDist)
                {
                    tempBestBcDist = bcDist;
                    tempBestBcDistCharge = row + MinCharge;
                }

                if (summedMedianIntensity[row] > 0)
                    EnvelopeIntensityScoreAcrossCharge[row] = summedPeakIntensity[row] / summedMedianIntensity[row];
            }
            RepresentativeSummedEnvelop = summedIntensity[tempBestBcDistCharge - MinCharge];

            // normalize abudnace across charges
            var s = AbundanceDistributionAcrossCharge.Sum();
            for (var i = 0; i < AbundanceDistributionAcrossCharge.Length; i++)
            {
                AbundanceDistributionAcrossCharge[i] = AbundanceDistributionAcrossCharge[i] / s;
            }

            // set representative charge, mz and scanNum
            var repEnvelope = bestEnvelopeAcrossCharge[tempBestBcDistCharge - MinCharge];
            RepresentativeCharge = repEnvelope.Charge;
            RepresentativeMz = repEnvelope.RepresentativePeak.Mz;
            RepresentativeScanNum = repEnvelope.ScanNum;
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

        public void AddObservedEnvelope(ObservedIsotopeEnvelope envelope)
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
            Envelopes.Add(envelope);
        }

        public void Clear()
        {
            MaxScanNum = -1;
            MinScanNum = -1;
            MinCharge = -1;
            MaxCharge = -1;
            Envelopes.Clear();
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
        
        public IEnumerable<Ms1Peak> GetMajorPeaks()
        {
            foreach (var envelope in Envelopes)
            {
                foreach (var k in _theoreticalIsotopeEnvelope.IndexOrderByRanking)
                {
                    if (_theoreticalIsotopeEnvelope.Isotopes[k].Ratio < 0.4) yield break;
                    var peak = envelope.Peaks[k];
                    if (peak != null && peak.Active) yield return envelope.Peaks[k];
                }
            }
        }

        public IEnumerable<Ms1Peak> GetMinorPeaks()
        {
            foreach (var envelope in Envelopes)
            {
                for(var j = _theoreticalIsotopeEnvelope.Size - 1; j >= 0; j--)
                {
                    var k = _theoreticalIsotopeEnvelope.IndexOrderByRanking[j];
                    if (_theoreticalIsotopeEnvelope.Isotopes[k].Ratio >= 0.4) yield break;
                    
                    var peak = envelope.Peaks[k];
                    if (peak != null && peak.Active) yield return envelope.Peaks[k];
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
        
        public List<ObservedIsotopeEnvelope> Envelopes { get; private set; }
        public double Score { get; internal set; }
        
        public double[] RepresentativeSummedEnvelop { get; internal set; }
        
        public double[] EnvelopeDistanceScoreAcrossCharge { get; internal set; }
        public double[] EnvelopeCorrelationScoreAcrossCharge { get; internal set; }
        public double[] EnvelopeIntensityScoreAcrossCharge { get; internal set; }
        public double[] AbundanceDistributionAcrossCharge { get; internal set; }
        public double[] BestCorrelationScoreAcrossCharge { get; private set; }
        public double[] BestDistanceScoreAcrossCharge { get; private set; }
        public double[] BestIntensityScoreAcrossCharge { get; private set; }

        public double BestCorrelationScore { get { return BestCorrelationScoreAcrossCharge.Max();  } }
        public double BestDistance { get { return BestDistanceScoreAcrossCharge.Where(d => d > 0).Min(); } }
        public double BestIntensityScore { get { return BestIntensityScoreAcrossCharge.Max(); } }

        /*
        public double BestSummedEnvelopeDistanceEvenCharge { get; private set; }
        public double BestSummedEnvelopeDistanceOddCharge { get; private set; }
        public double BestEnvelopeDistance { get; private set; }
        public double BestEnvelopeCorrelation { get; private set; }
       */

        internal double SummedEnvelopeBcDistCutoff = 0.1;

        public double tempInitialCorr = 0;
        public double tempInitialDist = 0;

        private readonly TheoreticalIsotopeEnvelope _theoreticalIsotopeEnvelope;

        internal byte Flag;
    }
}
