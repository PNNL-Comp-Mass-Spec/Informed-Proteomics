using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;

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
        
        /*
        public void UpdateScores(IList<Ms1Spectrum> spectra)
        {
            if (Envelopes.Count < 1) return;

            var bestEnvelopeCorrelation = 0.0d;
            var bestBhattacharyyaDistance = 1.0d;
            var bestRankSumScore = 1.0d;
            var bestPoissonScore = 1.0d;

            double goodEnvCorrTh;
            double goodEnvBcTh;

            if (RepresentativeMass < 15000)
            {
                goodEnvCorrTh = 0.6;
                goodEnvBcTh = 0.25;
            }
            else if (RepresentativeMass < 25000)
            {
                goodEnvCorrTh = 0.4;
                goodEnvBcTh = 0.3;
            }
            else
            {
                goodEnvCorrTh = 0.3;
                goodEnvBcTh = 0.3;
            }

            ObservedEnvelope repEnvelope = null;
            var envelopePerTime = new double[ColumnLength][];
            var envelopePerCharge = new double[ChargeLength][];
            for (var i = 0; i < ColumnLength; i++) envelopePerTime[i] = new double[TheoreticalEnvelope.Count];
            for (var i = 0; i < ChargeLength; i++) envelopePerCharge[i] = new double[TheoreticalEnvelope.Count];
            GoodEnvelopeCount = 0;

            var memberCorr = new double[Envelopes.Count];
            var memberBc = new double[Envelopes.Count];

            for (var e = 0; e < Envelopes.Count; e++)
            {
                var envelope = Envelopes[e];
                var spectrum = spectra[envelope.Col];
                var statSigTestResult = spectrum.TestStatisticalSignificance(TheoreticalEnvelope, envelope);
                var envCorr = envelope.GetPearsonCorrelation(TheoreticalEnvelope.Envelope);
                var bcDistance = envelope.GetBhattacharyyaDistance(TheoreticalEnvelope.EnvelopePdf);

                for (var i = 0; i < TheoreticalEnvelope.Count; i++)
                {
                    if (envelope.Peaks[i] == null || !envelope.Peaks[i].Active) continue;
                    envelopePerTime[envelope.Col - MinCol][i] += envelope.Peaks[i].Intensity;
                    envelopePerCharge[envelope.Row - MinRow][i] += envelope.Peaks[i].Intensity;
                }

                if (envCorr > goodEnvCorrTh && bcDistance < goodEnvBcTh && statSigTestResult.PoissonScore > LogP2 && statSigTestResult.RankSumScore > LogP2)
                {
                    if (!_scoreInit)
                    {
                        if (spectrum.CorrectChargeState(envelope, _minCharge))
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
                        repEnvelope = envelope;
                    }
                }

                memberCorr[e] = envCorr;
                memberBc[e] = bcDistance;
            }

            if (GoodEnvelopeCount < 1 || repEnvelope == null) return;

            var repPeak = repEnvelope.Peaks[repEnvelope.RefIsotopeInternalIndex];
            RepresentativeCharge = _minCharge + repEnvelope.Row;
            RepresentativeMass = Ion.GetMonoIsotopicMass(repPeak.Mz, RepresentativeCharge, TheoreticalEnvelope[repEnvelope.RefIsotopeInternalIndex].Index);
            RepresentativeMz = repPeak.Mz;
            RepresentativeScanNum = Run.GetMs1ScanVector()[repEnvelope.Col];

            var bestBcPerTime = 10d;
            for (var col = MinCol; col <= MaxCol; col++)
            {
                var bc = TheoreticalEnvelope.GetBhattacharyyaDistance(envelopePerTime[col - MinCol]);
                if (bc < bestBcPerTime) bestBcPerTime = bc;
            }

            var bestBcPerCharge = 10d;
            var bestBcEvenCharge = 10d;
            var bestBcOddCharge = 10d;
            var abuPerCharge = new double[ChargeLength];
            var maxAbuRow = MinRow;

            //if (BcDistPerCharge == null || BcDistPerCharge.Length != ChargeLength) BcDistPerCharge = new double[ChargeLength];
            //else Array.Clear(BcDistPerCharge, 0, BcDistPerCharge.Length);

            for (var row = MinRow; row <= MaxRow; row++)
            {
                var bc = TheoreticalEnvelope.GetBhattacharyyaDistance(envelopePerCharge[row - MinRow]);
                //BcDistPerCharge[row - MinRow] = bc;
                if (bc < bestBcPerCharge) bestBcPerCharge = bc;

                if ((row + _minCharge) % 2 == 0)
                {
                    if (bc < bestBcEvenCharge) bestBcEvenCharge = bc;
                }
                else
                {
                    if (bc < bestBcOddCharge) bestBcOddCharge = bc;
                }
                abuPerCharge[row - MinRow] = envelopePerCharge[row - MinRow].Sum();

                if (abuPerCharge[maxAbuRow - MinRow] < abuPerCharge[row - MinRow]) maxAbuRow = row;
            }

            if (ChargeLength > 1)
            {
                var meanStd = abuPerCharge.MeanStandardDeviation();
                SetScore(Ms1FeatureScore.AbundanceChangesOverCharges, meanStd.Item2 / meanStd.Item1);
                SetScore(Ms1FeatureScore.BhattacharyyaDistanceSummedOverEvenCharges, bestBcEvenCharge);
                SetScore(Ms1FeatureScore.BhattacharyyaDistanceSummedOverOddCharges, bestBcOddCharge);
            }
            else
            {
                SetScore(Ms1FeatureScore.AbundanceChangesOverCharges, 0);
                SetScore(Ms1FeatureScore.BhattacharyyaDistanceSummedOverEvenCharges, bestBcPerCharge);
                SetScore(Ms1FeatureScore.BhattacharyyaDistanceSummedOverOddCharges, bestBcPerCharge);
            }

            var bestMzErrorPpm = 10d;
            var totalMzErrorPpm = 0d;
            var totalMzPairCount = 0;
            for (var i = 0; i < Envelopes.Count; i++)
            {
                if (memberBc[i] > 0.3 && memberCorr[i] < 0.5) continue;
                var envelope = Envelopes[i];
                var mzErrorPpm = 0d;
                var n = 0;
                var charge = _minCharge + envelope.Row;
                for (var j = 0; j < TheoreticalEnvelope.Count; j++)
                {
                    if (envelope.Peaks[j] == null || !envelope.Peaks[j].Active) continue;
                    var theoreticalMz = Ion.GetIsotopeMz(RepresentativeMass, charge, TheoreticalEnvelope[j].Index);
                    mzErrorPpm += (Math.Abs(envelope.Peaks[j].Mz - theoreticalMz) * 1e6) / theoreticalMz;
                    n++;
                }

                totalMzErrorPpm += mzErrorPpm;
                totalMzPairCount += n;
                mzErrorPpm /= n;
                if (mzErrorPpm < bestMzErrorPpm) bestMzErrorPpm = mzErrorPpm;
            }

            SetScore(Ms1FeatureScore.TotalMzError, totalMzErrorPpm / totalMzPairCount);
            SetScore(Ms1FeatureScore.EnvelopeCorrelation, bestEnvelopeCorrelation);
            SetScore(Ms1FeatureScore.RankSum, bestRankSumScore);
            SetScore(Ms1FeatureScore.Poisson, bestPoissonScore);
            SetScore(Ms1FeatureScore.BhattacharyyaDistance, bestBhattacharyyaDistance);
            SetScore(Ms1FeatureScore.BhattacharyyaDistanceSummedOverCharges, bestBcPerCharge);
            SetScore(Ms1FeatureScore.BhattacharyyaDistanceSummedOverTimes, bestBcPerTime);
            SetScore(Ms1FeatureScore.MzError, bestMzErrorPpm);
            SetSummedEnvelope(memberCorr, memberBc);
            Probability = GetProbabilityByLogisticRegression();

            _scoreInit = true;
        }
        */
        private readonly List<Ms1Spectrum> _ms1Spectra;
        private readonly LcMsPeakScorer[] _peakScorerForLargeMassFeature;
        private readonly LcMsPeakScorer[] _peakScorer;

    }
}
