using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MathAndStats;
using InformedProteomics.FeatureFinding.Data;
using InformedProteomics.FeatureFinding.IsotopicEnvelope;
using MathNet.Numerics.Distributions;

namespace InformedProteomics.FeatureFinding.Scoring
{
    public class LcMsPeakScorer
    {
        public LcMsPeakScorer(Ms1Spectrum spec, int numBits4WinSize = 19) // 19 bits -> 4096 ppm,  20 bits -> 2048
        {
            Spectrum = spec;
            var peaks = spec.Peaks;
            _windowComparer = new MzComparerWithBinning(numBits4WinSize);
            _minBinNum = _windowComparer.GetBinNumber(peaks[0].Mz);
            _maxBinNum = _windowComparer.GetBinNumber(peaks[peaks.Length - 1].Mz);
            var numberOfBins = _maxBinNum - _minBinNum + 1;

            _peakStartIndex = new int[2][];
            _peakRanking = new int[2][][];
            var intensities = new List<double>[2][];

            for (var i = 0; i < 2; i++)
            {
                _peakStartIndex[i] = new int[numberOfBins];
                _peakRanking[i] = new int[numberOfBins][];
                intensities[i] = new List<double>[numberOfBins];

                for (var j = 0; j < numberOfBins; j++)
                {
                    _peakStartIndex[i][j] = peaks.Length - 1;
                    intensities[i][j] = new List<double>();
                }
            }

            for (var i = 0; i < peaks.Length; i++)
            {
                var binNum = _windowComparer.GetBinNumber(peaks[i].Mz);
                var binMzAverage = _windowComparer.GetMzAverage(binNum);
                var binIdx = binNum - _minBinNum;

                intensities[0][binIdx].Add(peaks[i].Intensity);
                if (i < _peakStartIndex[0][binIdx])
                {
                    _peakStartIndex[0][binIdx] = i;
                }

                if (peaks[i].Mz < binMzAverage)
                {
                    intensities[1][binIdx].Add(peaks[i].Intensity);
                    if (i < _peakStartIndex[1][binIdx])
                    {
                        _peakStartIndex[1][binIdx] = i;
                    }
                }
                else if (binNum < _maxBinNum) // skip this at the rightmost bin
                {
                    intensities[1][binIdx + 1].Add(peaks[i].Intensity);
                    if (i < _peakStartIndex[1][binIdx + 1])
                    {
                        _peakStartIndex[1][binIdx + 1] = i;
                    }
                }
            }

            for (var i = 0; i < 2; i++)
            {
                for (var binIdx = 0; binIdx < numberOfBins; binIdx++)
                {
                    if (intensities[i][binIdx].Count < 1)
                    {
                        continue;
                    }

                    _peakRanking[i][binIdx] = GetRankings(intensities[i][binIdx].ToArray(), out _);
                }
            }
        }

        public bool CheckChargeState(ObservedIsotopeEnvelope envelope)
        {
            var checkCharge = envelope.Charge;
            if (checkCharge > 20)
            {
                return true; //high charge (> +20), just pass
            }

            var peakStartIndex = envelope.MinMzPeak.IndexInSpectrum;
            var peakEndIndex = envelope.MaxMzPeak.IndexInSpectrum;
            var nPeaks = peakEndIndex - peakStartIndex + 1;

            if (nPeaks < 10)
            {
                return false;
            }

            if (envelope.NumberOfPeaks > nPeaks * 0.7)
            {
                return true;
            }

            var tolerance = new Tolerance(5);
            var threshold = nPeaks * 0.5;
            var mzTol = tolerance.GetToleranceAsMz(Spectrum.Peaks[peakStartIndex].Mz);

            var minCheckCharge = Math.Max(checkCharge * 2 - 1, 4);
            var maxCheckCharge = Math.Min(checkCharge * 5 + 1, 60);
            var maxDeltaMz = Constants.C13MinusC12 / minCheckCharge + mzTol;
            var nChargeGaps = new int[maxCheckCharge - minCheckCharge + 1];

            for (var i = peakStartIndex; i <= peakEndIndex; i++)
            {
                for (var j = i + 1; j <= peakEndIndex; j++)
                {
                    var deltaMz = Spectrum.Peaks[j].Mz - Spectrum.Peaks[i].Mz;

                    if (deltaMz > maxDeltaMz)
                    {
                        break;
                    }

                    for (var c = Math.Round(1 / (deltaMz + mzTol)); c <= Math.Round(1 / (deltaMz - mzTol)); c++)
                    {
                        if (c < minCheckCharge || c > maxCheckCharge)
                        {
                            continue;
                        }

                        var k = (int)c - minCheckCharge;
                        nChargeGaps[k]++;

                        if (nChargeGaps[k] + 1 > threshold && nChargeGaps[k] + 1 > 1.25 * envelope.NumberOfPeaks)
                        {
                            return false;
                        }
                    }
                }
            }

            return true;
        }

        public IsotopeEnvelopeStatisticalInfo PreformStatisticalSignificanceTest(ObservedIsotopeEnvelope envelope)
        {
            //var refPeak = envelope.Peaks[envelope.RefIsotopeInternalIndex];

            double mostAbuMz;
            var mostAbutPeakInternalIndex = envelope.TheoreticalEnvelope.IndexOrderByRanking[0];
            if (envelope.Peaks[mostAbutPeakInternalIndex] != null)
            {
                mostAbuMz = envelope.Peaks[mostAbutPeakInternalIndex].Mz;
            }
            else
            {
                mostAbuMz = envelope.TheoreticalEnvelope.GetIsotopeMz(envelope.Charge, mostAbutPeakInternalIndex);
            }

            var rankings = GetLocalRankings(mostAbuMz, out var peakStartIndex, out var mzBoundary);

            // smallest delta_mz = 0.01 (th) ?
            var ret = new IsotopeEnvelopeStatisticalInfo
            {
                LocalMzStart = mzBoundary.Item1,
                LocalMzEnd = mzBoundary.Item2,
                NumberOfLocalPeaks = rankings.Length,
                NumberOfPossiblePeaks = (int)Math.Ceiling(100 * (mzBoundary.Item2 - mzBoundary.Item1)),
                NumberOfIsotopePeaks = envelope.Size,
            };

            // calculate rankSum test score
            var rankSum = 0;
            var nRankSum = 0;
            for (var i = 0; i < envelope.Size; i++)
            {
                if (envelope.Peaks[i] == null || !envelope.Peaks[i].Active)
                {
                    continue;
                }

                ret.NumberOfMatchedIsotopePeaks++;

                //if (isotopeList[i].Ratio > RelativeIntensityThresholdForRankSum)
                //{
                var localIndex = envelope.Peaks[i].IndexInSpectrum - peakStartIndex;
                if (localIndex >= rankings.Length || localIndex < 0)
                {
                    continue;
                }

                rankSum += rankings[localIndex];
                nRankSum++;
                //}
            }

            var pValue = FitScoreCalculator.GetRankSumPValue(ret.NumberOfLocalPeaks, nRankSum, rankSum);
            ret.RankSumScore = (pValue > 0) ? -Math.Log(pValue, 2) : 50;

            // calculate poisson test score
            var n = ret.NumberOfPossiblePeaks;
            var k = ret.NumberOfIsotopePeaks; // # of theoretical isotope ions of the mass within the local window
            var n1 = ret.NumberOfLocalPeaks; // # of detected ions within the local window
            var k1 = ret.NumberOfMatchedIsotopePeaks; // # of matched ions generating isotope envelope profile

            var lambda = (n1 / (double)n) * k;
            pValue = 1 - Poisson.CDF(lambda, k1);
            ret.PoissonScore = (pValue > 0) ? -Math.Log(pValue, 2) : 50;
            return ret;
        }

        private int[] GetLocalRankings(double mz, out int peakStartIndex, out Tuple<double, double> mzBoundary)
        {
            var binNum = _windowComparer.GetBinNumber(mz);
            var binIndex = binNum - _minBinNum;
            byte binShift = 0;

            var d0 = Math.Abs(_windowComparer.GetMzAverage(binNum) - mz);
            var d1 = Math.Abs(_windowComparer.GetMzStart(binNum) - mz);
            var d2 = Math.Abs(_windowComparer.GetMzEnd(binNum) - mz);

            if (d1 < d2 && d1 < d0)
            {
                binShift = 1;
            }
            else if (d2 < d1 && d2 < d0 && binNum < _maxBinNum)
            {
                binShift = 1;
                binIndex++;
            }

            peakStartIndex = _peakStartIndex[binShift][binIndex];
            mzBoundary = GetMzBoundary(binShift, binNum);

            return _peakRanking[binShift][binIndex];
        }

        private Tuple<double, double> GetMzBoundary(byte binShift, int binNum)
        {
            var minMz = (binShift == 1) ? _windowComparer.GetMzAverage(binNum - 1) : _windowComparer.GetMzStart(binNum);
            var maxMz = (binShift == 1) ? _windowComparer.GetMzAverage(binNum) : _windowComparer.GetMzEnd(binNum);
            return new Tuple<double, double>(minMz, maxMz);
        }

        private int[] GetRankings(double[] values, out double medianValue)
        {
            var index = Enumerable.Range(0, values.Length).ToArray();
            Array.Sort(values, index);
            medianValue = values[values.Length / 2];

            var ranking = 1;
            var rankingList = new int[index.Length];
            for (var i = index.Length - 1; i >= 0; i--)
            {
                rankingList[index[i]] = ranking++;
            }
            return rankingList;
        }

        public readonly Ms1Spectrum Spectrum;
        private readonly MzComparerWithBinning _windowComparer;
        private readonly int[][] _peakStartIndex;
        private readonly int[][][] _peakRanking;
        private readonly int _minBinNum;
        private readonly int _maxBinNum;
    }

    public class IsotopeEnvelopeStatisticalInfo
    {
        public int NumberOfPossiblePeaks { get; internal set; } // n
        public int NumberOfIsotopePeaks { get; internal set; } // k
        public int NumberOfLocalPeaks { get; internal set; } // n1
        public int NumberOfMatchedIsotopePeaks { get; internal set; } // k1

        public uint IsotopePeakRankSum { get; internal set; }

        public double RankSumScore { get; internal set; }
        public double PoissonScore { get; internal set; }

        public double LocalMzStart { get; internal set; }
        public double LocalMzEnd { get; internal set; }

        public IsotopeEnvelopeStatisticalInfo()
        {
            NumberOfPossiblePeaks = 0;
            NumberOfLocalPeaks = 0;
            NumberOfIsotopePeaks = 0;
            NumberOfMatchedIsotopePeaks = 0;

            LocalMzStart = 0d;
            LocalMzEnd = 0d;
            IsotopePeakRankSum = 0;
            RankSumScore = 0d;
            PoissonScore = 0d;
        }
    }
}
