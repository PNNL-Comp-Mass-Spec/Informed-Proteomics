using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsFeatureScorer
    {
        public LcMsFeatureScorer(List<Ms1Spectrum> spectra)
        {
            _ms1Spectra = spectra;
            _peakScorer = new LcMsPeakScorer[_ms1Spectra.Count];
            _peakScorerForLargeMassFeature = new LcMsPeakScorer[_ms1Spectra.Count];

            for (var i = 0; i < _ms1Spectra.Count; i++)
            {
                _peakScorer[i] = new LcMsPeakScorer(_ms1Spectra[i], 20);
                _peakScorerForLargeMassFeature[i] = new LcMsPeakScorer(_ms1Spectra[i], 19);
            }
        }
        private static readonly double LogP2 = -Math.Log(0.02, 2);
        public LcMsFeatureScore ComputeScores(LcMsPeakCluster feature, bool initialEvaluation = false)
        {
            if (feature.Envelopes.Count < 1) return null;

            var bestEnvelopeCorrelation = 0.0d;
            var bestBhattacharyyaDistance = 1.0d;
            var bestRankSumScore = 1.0d;
            var bestPoissonScore = 1.0d;

            double goodEnvCorrTh;
            double goodEnvBcTh;

            if (feature.RepresentativeMass < 15000)
            {
                goodEnvCorrTh = 0.6;
                goodEnvBcTh = 0.25;
            }
            else if (feature.RepresentativeMass < 25000)
            {
                goodEnvCorrTh = 0.4;
                goodEnvBcTh = 0.3;
            }
            else
            {
                goodEnvCorrTh = 0.3;
                goodEnvBcTh = 0.3;
            }

            //ObservedIsotopeEnvelope repEnvelope = null;
            //var envelopePerTime = new double[feature.ColumnLength][];
            //var envelopePerCharge = new double[feature.ChargeLength][];
            //for (var i = 0; i < ColumnLength; i++) envelopePerTime[i] = new double[feature.Envelopes[0].Size];
            //for (var i = 0; i < ChargeLength; i++) envelopePerCharge[i] = new double[feature.Envelopes[0].Size];
            var abuPerCharge = new double[feature.ChargeLength];
            var GoodEnvelopeCount = 0;
            var memberCorr = new double[feature.Envelopes.Count];
            var memberBc = new double[feature.Envelopes.Count];
            var theoreticalEnvelope = feature.Envelopes[0].TheoreticalEnvelope;

            for (var e = 0; e < feature.Envelopes.Count; e++)
            {
                var envelope = feature.Envelopes[e];
                var spectrum = _ms1Spectra[envelope.RepresentativePeak.Ms1SpecIndex];

                abuPerCharge[envelope.Charge - feature.MinCharge] += envelope.Abundance;

                var peakScorer = (feature.Mass < 15000) ? _peakScorer[envelope.RepresentativePeak.Ms1SpecIndex] : _peakScorerForLargeMassFeature[envelope.RepresentativePeak.Ms1SpecIndex];
                //var statSigTestResult = spectrum.TestStatisticalSignificance(TheoreticalEnvelope, envelope);

                var statSigTestResult = peakScorer.PreformStatisticalSignificanceTest(envelope);

                var envCorr = envelope.GetPearsonCorrelation(theoreticalEnvelope);
                var bcDistance = envelope.GetBhattacharyyaDistance(theoreticalEnvelope);
                /*
                for (var i = 0; i < theoreticalEnvelope.Size; i++)
                {
                    if (envelope.Peaks[i] == null || !envelope.Peaks[i].Active) continue;
                    envelopePerTime[envelope.Col - MinCol][i] += envelope.Peaks[i].Intensity;
                    envelopePerCharge[envelope.Row - MinRow][i] += envelope.Peaks[i].Intensity;
                }
                */

                if (envCorr > goodEnvCorrTh && bcDistance < goodEnvBcTh && statSigTestResult.PoissonScore > LogP2 && statSigTestResult.RankSumScore > LogP2)
                {
                    if (initialEvaluation)
                    {
                        if (peakScorer.CheckChargeState(envelope))
                        {
                            if (envCorr > 0.6 && bcDistance < 0.2) envelope.GoodEnough = true;
                            GoodEnvelopeCount++;
                        }
                    }
                    else
                    {
                        GoodEnvelopeCount++;
                    }

                    if (bestPoissonScore < statSigTestResult.PoissonScore) bestPoissonScore = statSigTestResult.PoissonScore;
                    if (bestRankSumScore < statSigTestResult.RankSumScore) bestRankSumScore = statSigTestResult.RankSumScore;
                    if (bestEnvelopeCorrelation < envCorr) bestEnvelopeCorrelation = envCorr;
                    if (bcDistance < bestBhattacharyyaDistance)
                    {
                        bestBhattacharyyaDistance = bcDistance;
                        //repEnvelope = envelope;
                    }
                }

                memberCorr[e] = envCorr;
                memberBc[e] = bcDistance;
            }

            //if (GoodEnvelopeCount < 1 || repEnvelope == null) return null;
            if (GoodEnvelopeCount < 1) return null;

            /*var repPeak = repEnvelope.Peaks[repEnvelope.RefIsotopeInternalIndex];
            var repPeak = repEnvelope.RepresentativePeak;
            var representativeCharge = repEnvelope.Charge;
            var representativeMass = Ion.GetMonoIsotopicMass(repPeak.Mz, RepresentativeCharge, TheoreticalEnvelope[repEnvelope.RefIsotopeInternalIndex].Index);
            var representativeMz = repPeak.Mz;
            var representativeScanNum = Run.GetMs1ScanVector()[repEnvelope.Col];*/

            var bestBcPerTime = feature.EnvelopeDistanceScoreAcrossTime.Min();
            /*
            var bestBcPerTime = 10d;
            for (var col = MinCol; col <= MaxCol; col++)
            {
                var bc = TheoreticalEnvelope.GetBhattacharyyaDistance(envelopePerTime[col - MinCol]);
                if (bc < bestBcPerTime) bestBcPerTime = bc;
            }*/

            var bestBcPerCharge = 10d;
            var bestBcEvenCharge = 10d;
            var bestBcOddCharge = 10d;
            var maxAbuRow = 0;

            for (var c = feature.MinCharge; c <= feature.MaxCharge; c++)
            {
                var row = c - feature.MinCharge;
                var bc = feature.EnvelopeDistanceScoreAcrossCharge[row];
                if (bc < bestBcPerCharge) bestBcPerCharge = bc;

                if (c%2 == 0)
                {
                    if (bc < bestBcEvenCharge) bestBcEvenCharge = bc;
                }
                else
                {
                    if (bc < bestBcOddCharge) bestBcOddCharge = bc;
                }
                if (abuPerCharge[maxAbuRow] < abuPerCharge[row]) maxAbuRow = row;
            }

            var score = new LcMsFeatureScore();
            if (feature.ChargeLength > 1)
            {
                var meanStd = abuPerCharge.MeanStandardDeviation();
                score.AbundanceChangesOverCharges = meanStd.Item2/meanStd.Item1;
                score.BhattacharyyaDistanceSummedOverEvenCharges = bestBcEvenCharge;
                score.BhattacharyyaDistanceSummedOverOddCharges = bestBcOddCharge;
            }
            else
            {
                score.AbundanceChangesOverCharges = 0;
                score.BhattacharyyaDistanceSummedOverEvenCharges = bestBcPerCharge;
                score.BhattacharyyaDistanceSummedOverOddCharges = bestBcPerCharge;
            }


            score.EnvelopeCorrelation = bestEnvelopeCorrelation;
            score.RankSum = bestRankSumScore;
            score.Poisson = bestPoissonScore;
            score.BhattacharyyaDistance = bestBhattacharyyaDistance;
            score.BhattacharyyaDistanceSummedOverCharges = bestBcPerCharge;
            score.BhattacharyyaDistanceSummedOverTimes = bestBcPerTime;
            
            //SetSummedEnvelope(memberCorr, memberBc);
            return score;
        }
        
        private readonly List<Ms1Spectrum> _ms1Spectra;
        private readonly LcMsPeakScorer[] _peakScorerForLargeMassFeature;
        private readonly LcMsPeakScorer[] _peakScorer;

    }
}
