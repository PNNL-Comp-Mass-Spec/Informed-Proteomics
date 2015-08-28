﻿using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.Backend.Data.Spectrometry
{

    public class Ms1Spectrum : Spectrum
    {
        public readonly double MedianIntensity;

        public double MinMz
        {
            get { return Peaks[0].Mz; }
        }

        public double MaxMz
        {
            get { return Peaks[Peaks.Length - 1].Mz; }
        }

        public Ms1Spectrum(int scanNum, int index, Ms1Peak[] peaks) : base(scanNum)
        {
            Index = index;
            Peaks = new Peak[peaks.Length];
            peaks.CopyTo(Peaks, 0);
            MedianIntensity = Peaks.Select(p => p.Intensity).Median();
            PreArrangeLocalMzWindows();
        }

        public Ms1Spectrum(int scanNum, int index, Peak[] peaks) : base(scanNum)
        {
            Index = index;
            MsLevel = 1;
            Peaks = new Peak[peaks.Length];
            var sIndex = (ushort) index;
            for (var i = 0; i < Peaks.Length; i++)
                Peaks[i] = new Ms1Peak(peaks[i].Mz, peaks[i].Intensity, i) {Ms1SpecIndex = sIndex};
            MedianIntensity = Peaks.Select(p => p.Intensity).Median();
            PreArrangeLocalMzWindows();
        }

        public Ms1Peak[] GetAllIsotopePeaks(double monoIsotopeMass, int charge, TheoreticalIsotopeEnvelope isotopeList, Tolerance tolerance)
        {
            var observedPeaks = new Ms1Peak[isotopeList.Size];
            var mz = isotopeList.GetIsotopeMz(charge, 0);

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
                    observedPeaks[0] = (Ms1Peak) Peaks[bestPeakIndex];
                }
                ++i;
            }

            var peakIndex = (bestPeakIndex >= 0) ? bestPeakIndex + 1 : index;
            // go up
            for (var j = 1; j < isotopeList.Size; j++)
            {
                var isotopeMz = isotopeList.GetIsotopeMz(charge, j);
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
                    if (peakMz >= minMz) // find match, move to prev isotope
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
      
        internal int Index { private set; get; }

        private int[][] _peakStartIndex;
        private int[][][] _peakRanking;
        
        private double[][] _medianIntensity;
        private double[][] _highestIntensity;
        private int[][] _intensePeakCount;

        private const double MzWindowSize = 6;
        
        private void PreArrangeLocalMzWindows()
        {
            var numberOfbins = (int)Math.Round((MaxMz - MinMz) / MzWindowSize) + 1;

            _peakStartIndex     = new int[2][];
            _medianIntensity = new double[2][];
            _highestIntensity = new double[2][];
            _peakRanking        = new int[2][][];
            _intensePeakCount = new int[2][];
            var intensities     = new List<double>[2][];

            for (var i = 0; i < 2; i++)
            {
                _peakStartIndex[i]      = new int[numberOfbins];
                _medianIntensity[i]     = new double[numberOfbins];
                _highestIntensity[i]    = new double[numberOfbins];
                _peakRanking[i]         = new int[numberOfbins][];
                _intensePeakCount[i] = new int[numberOfbins];

                intensities[i]          = new List<double>[numberOfbins];

                for (var j = 0; j < numberOfbins; j++)
                {
                    _peakStartIndex[i][j]   = Peaks.Length - 1;
                    intensities[i][j]       = new List<double>();
                }
            }
            
            for (var i = 0; i < Peaks.Length; i++)
            {
                var binIdx = (int) Math.Round((Peaks[i].Mz - MinMz)/MzWindowSize);
                var binCenterMz = MinMz + MzWindowSize*binIdx;

                intensities[0][binIdx].Add(Peaks[i].Intensity);
                if (i < _peakStartIndex[0][binIdx]) _peakStartIndex[0][binIdx] = i;

                if (Peaks[i].Mz < binCenterMz)
                {
                    intensities[1][binIdx].Add(Peaks[i].Intensity);
                    if (i < _peakStartIndex[1][binIdx]) _peakStartIndex[1][binIdx] = i;
                }
                else if (binIdx < numberOfbins - 1) // skip this at the rightmost bin
                {
                    intensities[1][binIdx + 1].Add(Peaks[i].Intensity);
                    if (i < _peakStartIndex[1][binIdx + 1]) _peakStartIndex[1][binIdx + 1] = i;
                }
            }

            for (var i = 0; i < 2; i++)
            {
                for (var binIdx = 0; binIdx < numberOfbins; binIdx++)
                {
                    if (intensities[i][binIdx].Count < 1) continue;

                    double medianIntensity, highestIntensity;
                    _peakRanking[i][binIdx] = GetRankings(intensities[i][binIdx].ToArray(), out medianIntensity, out highestIntensity);
                    _medianIntensity[i][binIdx] = medianIntensity;
                    _highestIntensity[i][binIdx] = highestIntensity;

                    var intensePeakThreshold = highestIntensity * 0.1;

                    _intensePeakCount[i][binIdx] = intensities[i][binIdx].Count(x => x > intensePeakThreshold);
                }
            }
        }

        public LocalMzWindow GetLocalMzWindow(double mz)
        {
            var binIndex = (int)Math.Round((mz - MinMz) / MzWindowSize);

            var binCenterMz = MinMz + MzWindowSize * binIndex;
            var binStartMz = binCenterMz - MzWindowSize*0.5;
            var binEndMz = binCenterMz + MzWindowSize * 0.5;

            var numberOfbins = (int)Math.Round((MaxMz - MinMz) / MzWindowSize) + 1;
            byte binShift = 0;

            var d0 = Math.Abs(binCenterMz - mz);
            var d1 = Math.Abs(binStartMz - mz);
            var d2 = Math.Abs(binEndMz - mz);

            if (d1 < d2 && d1 < d0)
            {
                binShift = 1;
            }
            else if (d2 < d1 && d2 < d0 && binIndex < numberOfbins - 1)
            {
                binShift = 1;
                binIndex++;
            }

            if (binShift == 1)
            {
                binStartMz = binCenterMz - MzWindowSize;
                binEndMz = binCenterMz;
            }

            if (binIndex < 0 || binIndex >= numberOfbins || _peakRanking[binShift][binIndex] == null)
            {
                var emptyWin = new LocalMzWindow()
                {
                    MinMz = binStartMz,
                    MaxMz = binEndMz,
                    PeakStartIndex = -1,
                    PeakCount = 0,
                    MedianIntensity = 0,
                    HighestIntensity =  0,
                    IntensePeakCount = 0,
                    PeakRanking = null,
                };
                return emptyWin;
            }
            
            var peakStartIndex = _peakStartIndex[binShift][binIndex];
            var numOfPeaks = _peakRanking[binShift][binIndex].Length;
            
            var window = new LocalMzWindow()
            {
                MinMz = binStartMz,
                MaxMz = binEndMz,
                PeakStartIndex = peakStartIndex,
                PeakCount = numOfPeaks,
                MedianIntensity = _medianIntensity[binShift][binIndex],
                HighestIntensity = _highestIntensity[binShift][binIndex],
                IntensePeakCount = _intensePeakCount[binShift][binIndex],
                PeakRanking = _peakRanking[binShift][binIndex],
            };

            return window;
        }
     
        private int[] GetRankings(double[] values, out double medianValue, out double highestValue)
        {
            var index = Enumerable.Range(0, values.Length).ToArray();
            Array.Sort(values, index);
            
            medianValue = values[values.Length/2];
            highestValue = values[index.Length - 1];

            var ranking = 1;
            var rankingList = new int[index.Length];
            for (var i = index.Length - 1; i >= 0; i--)
            {
                rankingList[index[i]] = ranking++;
            }
            return rankingList;
        }

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


    public class LocalMzWindow
    {
        public double MinMz { get; internal set; }
        public double MaxMz { get; internal set; }
        public int PeakStartIndex { get; internal set; }
        public int PeakCount { get; internal set; }
        public int PeakEndIndex { get { return PeakStartIndex + PeakCount - 1; }}

        public double MedianIntensity { get; internal set; }
        public double HighestIntensity { get; internal set; }

        public int[] PeakRanking { get; internal set; }

        public int IntensePeakCount { get; internal set; }

        public double GetRankSumTestPvalue(Ms1Peak[] peaks, int envelopeSize)
        {
            if (PeakRanking == null) return 1.0d;
            
            // calculate ranksum test score
            var ranksum = 0;
            var nRankSum = 0;
            for (var i = 0; i < envelopeSize; i++)
            {
                if (peaks[i] == null || !peaks[i].Active) continue;

                var localIndex = peaks[i].IndexInSpectrum - PeakStartIndex;
                if (localIndex >= PeakCount || localIndex < 0) continue;
                ranksum += PeakRanking[localIndex];
                nRankSum++;
            }

            var pvalue = FitScoreCalculator.GetRankSumPvalue(PeakCount, nRankSum, ranksum);
            return pvalue;
        }

        public double GetPoissonTestPvalue(Ms1Peak[] peaks, int envelopeSize)
        {
            if (PeakRanking == null) return 1.0d;

            var intensePeakThreshold = HighestIntensity*0.1;
            var numberOfMatchedIsotopePeaks = peaks.Count(p => p != null && p.Active && p.Intensity > intensePeakThreshold);
            var numberOfPossiblePeaks = (int) Math.Ceiling(100*(MaxMz - MinMz));

            // calculate poisson test score
            var n = numberOfPossiblePeaks;
            var k = envelopeSize; // # of theretical isotope ions of the mass within the local window
            //var n1 = PeakCount; // # of detected ions within the local window
            var n1 = IntensePeakCount;
            var k1 = numberOfMatchedIsotopePeaks; // # of matched ions generating isotope envelope profile

            var lambda = ((double)n1 / (double)n) * k;
            var pvalue = 1 - Poisson.CDF(lambda, k1);

            return pvalue;
        }
    }

}
