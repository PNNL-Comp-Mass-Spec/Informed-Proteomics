using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using Microsoft.Win32;

namespace InformedProteomics.Backend.Data.Spectrometry
{
   
    public class Ms1Spectrum : Spectrum
    {
        public readonly double MedianIntensity;
        public double MinMz { get { return Peaks[0].Mz; } }
        public double MaxMz { get { return Peaks[Peaks.Length - 1].Mz; } }
        public static readonly double RelativeSignificantIntesnityThreshold = 0.7d;
        public static readonly double RelativeIntesnityThresholdForRankSum = 0.1d;

        public Ms1Spectrum(int scanNum, int index, Ms1Peak[] peaks) : base(scanNum)
        {
            Index = index;
            Peaks = peaks;
            MedianIntensity = Peaks.Select(p => p.Intensity).Median();
            
            _peakRanker = new PeakRanker(Peaks, 19);
            //_peakRanker.SetLocalRanking((Ms1Peak[]) Peaks);
            
            _peakRanker2 = new PeakRanker(Peaks, 20);
            //peakRanker2.SetLocalRanking((Ms1Peak[])Peaks);
            
        }
        
        public Ms1Spectrum(int scanNum, ushort index, Peak[] peaks) : base(scanNum)
        {
            Index = index;
            MsLevel = 1;
            Peaks = new Ms1Peak[peaks.Length];
            for (var i = 0; i < Peaks.Length; i++) Peaks[i] = new Ms1Peak(peaks[i].Mz, peaks[i].Intensity, i) { Ms1SpecIndex = index };
            MedianIntensity = Peaks.Select(p => p.Intensity).Median();
            
            _peakRanker = new PeakRanker(Peaks, 19);
            //_peakRanker.SetLocalRanking((Ms1Peak[])Peaks);

            _peakRanker2 = new PeakRanker(Peaks, 20);
            //peakRanker2.SetLocalRanking((Ms1Peak[])Peaks);
        }

        public bool CorrectChargeState(ObservedEnvelope envelope, int minScanCharge)
        {
            var checkCharge = envelope.Row + minScanCharge;
            if (checkCharge > 20) return true; //high charge (> +20), just pass

            var peakStartIndex = envelope.MinMzPeak.IndexInSpectrum;
            var peakEndIndex = envelope.MaxMzPeak.IndexInSpectrum;
            var nPeaks = peakEndIndex - peakStartIndex + 1;
            
            if (nPeaks < 10) return false;
            if (envelope.NumberOfPeaks > nPeaks*0.7) return true;

            var tolerance = new Tolerance(5);
            var threshold = nPeaks*0.5;
            var mzTol = tolerance.GetToleranceAsTh(Peaks[peakStartIndex].Mz);

            var minCheckCharge = Math.Max(checkCharge * 2 - 1, 4);
            var maxCheckCharge = Math.Min(checkCharge * 5 + 1, 60);
            var maxDeltaMz = Constants.C13MinusC12 / minCheckCharge + mzTol;
            var nChargeGaps = new int[maxCheckCharge - minCheckCharge + 1];

            for (var i = peakStartIndex; i <= peakEndIndex; i++)
            {
                for (var j = i + 1; j <= peakEndIndex; j++)
                {
                    var deltaMz = Peaks[j].Mz - Peaks[i].Mz;

                    if (deltaMz > maxDeltaMz) break;
                    for (var c = Math.Round(1 / (deltaMz + mzTol)); c <= Math.Round(1 / (deltaMz - mzTol)); c++)
                    {
                        if (c < minCheckCharge || c > maxCheckCharge) continue;
                        var k = (int) c - minCheckCharge;
                        nChargeGaps[k]++;

                        if (nChargeGaps[k] + 1 > threshold && nChargeGaps[k] + 1 > 1.25*envelope.NumberOfPeaks) return false;
                    }
                }
            }

            return true;
        }

        public Ms1Peak[] GetAllIsotopePeaks(double monoIsotopeMass, int charge, IsotopeList isotopeList, Tolerance tolerance)
        {
            var observedPeaks = new Ms1Peak[isotopeList.Count];
            var mz = Ion.GetIsotopeMz(monoIsotopeMass, charge, isotopeList[0].Index);
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;

            var index = Array.BinarySearch(Peaks, new Ms1Peak(minMz, 0, 0));
            if (index < 0) index = ~index;

            var bestPeakIndex = -1;
            var bestIntensity = 0.0;
            // go up
            var i = index;
            while (i >= 0 && i < Peaks.Length)
            {
                if (Peaks[i].Mz >= maxMz) break;
                if (Peaks[i].Intensity > bestIntensity)
                {
                    bestIntensity = Peaks[i].Intensity;
                    bestPeakIndex = i;
                    observedPeaks[0] = (Ms1Peak)Peaks[bestPeakIndex];
                }
                ++i;
            }

            var peakIndex = (bestPeakIndex >= 0) ? bestPeakIndex + 1 : index;
            // go up
            for (var j = 1; j < isotopeList.Count; j++)
            {
                var isotopeIndex = isotopeList[j].Index;
                
                var isotopeMz = Ion.GetIsotopeMz(monoIsotopeMass, charge, isotopeIndex);
                tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                minMz = isotopeMz - tolTh;
                maxMz = isotopeMz + tolTh;
                for (i = peakIndex; i < Peaks.Length; i++)
                {
                    var peakMz = Peaks[i].Mz;
                    if (peakMz > maxMz)
                    {
                        peakIndex = i;
                        break;
                    }
                    if (peakMz >= minMz)    // find match, move to prev isotope
                    {
                        var peak = Peaks[i];
                        if (observedPeaks[j] == null ||
                            peak.Intensity > observedPeaks[j].Intensity)
                        {
                            observedPeaks[j] = (Ms1Peak) peak;
                        }
                    }
                }
            }
            return observedPeaks;            
        }
        /*
        public double GetLocalMedianIntensity(Ms1Peak peak, double targetMass)
        {
            if (targetMass < 15000)
            {
                return _peakRanker2.GetLocalMedianIntensity(peak.Mz);
            }
            return _peakRanker.GetLocalMedianIntensity(peak.Mz);
        }*/

        public LocalMzWindow GetLocalMzWindow(double mz, double targetMass)
        {
            if (targetMass < 15000)
            {
                return _peakRanker2.GetLocalMzWindow(mz);
            }
            return _peakRanker.GetLocalMzWindow(mz);            
        }
        
        public StatSigTestResult TestStatisticalSignificance(IsotopeList isotopeList, ObservedEnvelope envelope)
        {
            int peakStartIndex;
            Tuple<double, double> mzBoundary;
            //var mostAbuPeak = isotopePeaks[isotopeList.SortedIndexByIntensity[0]];
            //if (mostAbuPeak == null) return null;
            var refPeak = envelope.Peaks[envelope.RefIsotopeInternalIndex];
            var rankings = _peakRanker.GetLocalRankings(refPeak.Mz, out peakStartIndex, out mzBoundary);
            // smallest delta_mz = 0.01 (th) ?
            var ret = new StatSigTestResult
            {
                LocalMzStart = mzBoundary.Item1,
                LocalMzEnd = mzBoundary.Item2,
                NumberOfLocalPeaks = rankings.Length,
                NumberOfPossiblePeaks = (int)Math.Ceiling(100 * (mzBoundary.Item2 - mzBoundary.Item1)),
                NumberOfIsotopePeaks = isotopeList.Count
            };

            // calculate ranksum test score
            var ranksum = 0;
            var nRankSum = 0;
            for (var i = 0; i < isotopeList.Count; i++)
            {
                if (envelope.Peaks[i] == null || !envelope.Peaks[i].Active) continue;
                ret.NumberOfMatchedIsotopePeaks++;
                if (isotopeList[i].Ratio > RelativeIntesnityThresholdForRankSum)
                {
                    var localIndex = envelope.Peaks[i].IndexInSpectrum - peakStartIndex;
                    if (localIndex >= rankings.Length || localIndex < 0) continue;
                    ranksum += rankings[localIndex];
                    nRankSum++;
                }
            }

            var pvalue = FitScoreCalculator.GetRankSumPvalue(ret.NumberOfLocalPeaks, nRankSum, ranksum);
            ret.RankSumScore = (pvalue > 0) ? -Math.Log(pvalue, 2) : 50;

            // calculate poisson test score
            var n = ret.NumberOfPossiblePeaks;
            var k = ret.NumberOfIsotopePeaks; // # of theretical isotope ions of the mass within the local window
            var n1 = ret.NumberOfLocalPeaks; // # of detected ions within the local window
            var k1 = ret.NumberOfMatchedIsotopePeaks; // # of matched ions generating isotope envelope profile

            var lambda = ((double)n1 / (double)n) * k;
            pvalue = 1 - Poisson.CDF(lambda, k1);
            ret.PoissonScore = (pvalue > 0) ? -Math.Log(pvalue, 2) : 50;
            return ret;
        }
      
        internal int Index { private set; get; }
        
        private readonly PeakRanker _peakRanker;
        private readonly PeakRanker _peakRanker2;

    }

    public static class Ms1PeakExtension
    {
        public static void Shuffle<T>(this IList<T> list)
        {
            int n = list.Count;
            var rnd = new Random();
            while (n > 1)
            {
                int k = (rnd.Next(0, n) % n);
                n--;
                T value = list[k];
                list[k] = list[n];
                list[n] = value;
            }
        }

        public static void SumEnvelopeTo(this Ms1Peak[] peaks, double[] targetEnvelope)
        {
            for (var i = 0; i < targetEnvelope.Length; i++)
            {
                if (peaks[i] != null && peaks[i].Active) targetEnvelope[i] += peaks[i].Intensity;
            }
        }

        public static double Likelihood(this Ms1Peak[] peaks, IsotopeList isotopeList)
        {
            var ret = 0d;
            for (var i = 0; i < isotopeList.Count; i++)
            {
                if (peaks[i] != null && peaks[i].Active)
                {
                    ret += isotopeList.EnvelopePdf[i];
                }
            }
            return ret;
        }

        public static double GetChiSquareSignificanceScore(this Ms1Peak[] isotopePeaks, double[] theoreticalEnvelopePdf)
        {
            var k = theoreticalEnvelopePdf.Length - 1;
            var x = 0d;
            var s2 = 0d;

            for (var i = 0; i < theoreticalEnvelopePdf.Length; i++)
            {
                if (isotopePeaks[i] != null && isotopePeaks[i].Active)
                    s2 += isotopePeaks[i].Intensity;
            }
            for (var i = 0; i < theoreticalEnvelopePdf.Length; i++)
            {
                var p = theoreticalEnvelopePdf[i];
                var q = (isotopePeaks[i] != null && isotopePeaks[i].Active) ? isotopePeaks[i].Intensity / s2 : 0;
                x += ((p - q) * (p - q)) / (p + q);
            }

            x *= 0.5;
            var pvalue = ChiSquared.CDF(k, x);

            return (pvalue > 0) ? -Math.Log(pvalue, 2) : 50;
        }

        public static double GetPearsonCorrelation(this Ms1Peak[] isotopePeaks, double[] theoreticalEnvelope)
        {
            var m1 = 0.0;
            var m2 = 0.0;

            for (var i = 0; i < theoreticalEnvelope.Length; i++)
            {
                m1 += theoreticalEnvelope[i];
                if (isotopePeaks[i] != null && isotopePeaks[i].Active) m2 += isotopePeaks[i].Intensity;
            }

            m1 /= theoreticalEnvelope.Length;
            m2 /= theoreticalEnvelope.Length;

            // compute Pearson correlation
            var cov = 0.0;
            var s1 = 0.0;
            var s2 = 0.0;

            for (var i = 0; i < theoreticalEnvelope.Length; i++)
            {
                var d1 = theoreticalEnvelope[i] - m1;
                var d2 = (isotopePeaks[i] != null && isotopePeaks[i].Active) ? isotopePeaks[i].Intensity - m2 : -m2;
                cov += d1 * d2;
                s1 += d1 * d1;
                s2 += d2 * d2;
            }

            if (s1 <= 0 || s2 <= 0) return 0;

            return cov < 0 ? 0d : cov / Math.Sqrt(s1 * s2);
        }

        public static double GetBhattacharyyaDistance(this Ms1Peak[] isotopePeaks, double[] theoreticalEnvelopePdf)
        {
            var s2 = 0d;

            for (var i = 0; i < theoreticalEnvelopePdf.Length; i++)
            {
                if (isotopePeaks[i] != null && isotopePeaks[i].Active)
                    s2 += isotopePeaks[i].Intensity;
            }

            if (!(s2 > 0)) return IsotopeEnvelope.MaxBhattacharyyaDistance;

            var bc = 0d;
            for (var i = 0; i < theoreticalEnvelopePdf.Length; i++)
            {
                var p = theoreticalEnvelopePdf[i];
                var q = (isotopePeaks[i] != null && isotopePeaks[i].Active) ? isotopePeaks[i].Intensity / s2 : 0;
                bc += Math.Sqrt(p * q);
            }

            return -Math.Log(bc);
        }
    }

    public class StatSigTestResult
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
        
        public StatSigTestResult()
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

    public class PeakRanker
    {
        public PeakRanker(Peak[] peaks, int numBits4WinSize = 19) // 19 bits -> 4096 ppm,  20 bits -> 2048
        {
            _windowComparer = new MzComparerWithBinning(numBits4WinSize);
            _minBinNum = _windowComparer.GetBinNumber(peaks[0].Mz);
            _maxBinNum = _windowComparer.GetBinNumber(peaks[peaks.Length - 1].Mz);
            var numberOfbins = _maxBinNum - _minBinNum + 1;

            _peakStartIndex     = new int[2][];
            _medianIntensity = new double[2][];
            _peakRanking        = new int[2][][];
            var intensities     = new List<double>[2][];

            for (var i = 0; i < 2; i++)
            {
                _peakStartIndex[i]  = new int[numberOfbins];
                _medianIntensity[i] = new double[numberOfbins];
                _peakRanking[i]     = new int[numberOfbins][];
                intensities[i]      = new List<double>[numberOfbins];

                for (var j = 0; j < numberOfbins; j++)
                {
                    _peakStartIndex[i][j]   = peaks.Length - 1;
                    intensities[i][j]       = new List<double>();
                }
            }

            for (var i = 0; i < peaks.Length; i++)
            {
                var binNum = _windowComparer.GetBinNumber(peaks[i].Mz);
                var binMzAverage = _windowComparer.GetMzAverage(binNum);
                var binIdx = binNum - _minBinNum;

                intensities[0][binIdx].Add(peaks[i].Intensity);
                if (i < _peakStartIndex[0][binIdx]) _peakStartIndex[0][binIdx] = i;

                if (peaks[i].Mz < binMzAverage)
                {
                    intensities[1][binIdx].Add(peaks[i].Intensity);
                    if (i < _peakStartIndex[1][binIdx]) _peakStartIndex[1][binIdx] = i;
                }
                else if (binNum < _maxBinNum) // skip this at the rightmost bin
                {
                    intensities[1][binIdx + 1].Add(peaks[i].Intensity);
                    if (i < _peakStartIndex[1][binIdx + 1]) _peakStartIndex[1][binIdx + 1] = i;
                }
            }

            for (var i = 0; i < 2; i++)
            {
                for (var binIdx = 0; binIdx < numberOfbins; binIdx++)
                {
                    if (intensities[i][binIdx].Count < 1) continue;

                    double medianIntensity;
                    _peakRanking[i][binIdx] = GetRankings(intensities[i][binIdx].ToArray(), out medianIntensity);
                    _medianIntensity[i][binIdx] = medianIntensity;

                }
            }
        }
        /*
        public double GetLocalMedianIntensity(double mz) 
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

            return _medianIntensity[binShift][binIndex];

            //var peakStartIndex = _peakStartIndex[binShift][binIndex];
            //var peakEndIndex = peakStartIndex + _peakRanking[binShift][binIndex].Length - 1;
            //return new Tuple<int, int>(peakStartIndex, peakEndIndex);
        }
        */
        public int[] GetLocalRankings(double mz, out int peakStartIndex, out Tuple<double, double> mzBoundary)
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

        public LocalMzWindow GetLocalMzWindow(double mz)
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

            if (binIndex < 0 || binIndex >= _maxBinNum - _minBinNum + 1) return null;
            
            //var mzBoundary = GetMzBoundary(binShift, binIndex);
            var minMz = (binShift == 1) ? _windowComparer.GetMzAverage(binNum - 1) : _windowComparer.GetMzStart(binNum);
            var maxMz = (binShift == 1) ? _windowComparer.GetMzAverage(binNum) : _windowComparer.GetMzEnd(binNum);

            
            if (_peakRanking[binShift][binIndex] == null)
            {
                var emptyWin = new LocalMzWindow()
                {
                    MinMz = minMz,
                    MaxMz = maxMz,
                    PeakStartIndex = -1,
                    PeakCount = 0,
                    MedianIntensity = 0,
                };
                return emptyWin;
            }

            var peakStartIndex = _peakStartIndex[binShift][binIndex];
            var numOfPeaks = _peakRanking[binShift][binIndex].Length;

            var window = new LocalMzWindow()
            {
                MinMz = minMz,
                MaxMz = maxMz,
                PeakStartIndex = peakStartIndex,
                PeakCount = numOfPeaks,
                MedianIntensity = _medianIntensity[binShift][binIndex],
            };

            return window;
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
            
            
            
            medianValue = values[values.Length/2];

            var ranking = 1;
            var rankingList = new int[index.Length];
            for (var i = index.Length - 1; i >= 0; i--)
            {
                rankingList[index[i]] = ranking++;
            }
            return rankingList;
        }
      
        private readonly MzComparerWithBinning _windowComparer;
        private readonly int[][] _peakStartIndex;
        private readonly int[][][] _peakRanking;
        
        private readonly double[][] _medianIntensity;

        private readonly int _minBinNum;
        private readonly int _maxBinNum;
    }

    public class LocalMzWindow
    {
        public double MinMz { get; internal set; }
        public double MaxMz { get; internal set; }
        public int PeakStartIndex { get; internal set; }
        public int PeakCount { get; internal set; }
        public int PeakEndIndex { get { return PeakStartIndex + PeakCount - 1; }}

        public double MedianIntensity { get; internal set; }
    }

}
