using System;
using System.Collections.Generic;
using System.Drawing.Drawing2D;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Interpolation.Algorithms;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Spectrum
    {
        public Spectrum(IList<double> mzArr, IList<double> intensityArr, int scanNum)
        {
            Peaks = new Peak[mzArr.Count];
            for(var i=0; i<mzArr.Count; i++) Peaks[i] = new Peak(mzArr[i], intensityArr[i]);
            ScanNum = scanNum;
        }

        public Spectrum(ICollection<Peak> peaks, int scanNum)
        {
            Peaks = new Peak[peaks.Count];
            peaks.CopyTo(Peaks, 0);
            ScanNum = scanNum;
        }

        public int ScanNum { get; private set; }
        public int MsLevel
        {
            get { return _msLevel; }
            set { _msLevel = value; }
        }

        // Peaks are assumed to be sorted according to m/z
        public Peak[] Peaks { get; private set; }

        /// <summary>
        /// Finds the maximum intensity peak within the specified range
        /// </summary>
        /// <param name="mz">m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <returns>maximum intensity peak</returns>
        public Peak FindPeak(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;

            return FindPeak(minMz, maxMz);
        }

        /// <summary>
        /// Checks whether this spectrum contains all isotope peaks whose relative intensity is equal or larter than the threshold
        /// </summary>
        /// <param name="ion">ion</param>
        /// <param name="tolerance">tolerance</param>
        /// <param name="relativeIntensityThreshold">relative intensity threshold of the theoretical isotope profile</param>
        /// <returns>true if spectrum contains all ions; false otherwise.</returns>
        public bool ContainsIon(Ion ion, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var baseIsotopeIndex = ion.Composition.GetMostAbundantIsotopeZeroBasedIndex();
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelop();
            var baseIsotopMz = ion.GetIsotopeMz(baseIsotopeIndex);
            var baseIsotopePeakIndex = FindPeakIndex(baseIsotopMz, tolerance);
            if (baseIsotopePeakIndex < 0) return false;

            // go down
            var peakIndex = baseIsotopePeakIndex;
            for (var isotopeIndex = baseIsotopeIndex - 1; isotopeIndex >= 0; isotopeIndex--)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex - 1; i >= 0; i--)
                {
                    var peakMz = Peaks[i].Mz;
                    if (peakMz < minMz) return false;
                    if (peakMz <= maxMz)    // find match, move to prev isotope
                    {
                        peakIndex = i;
                        break;
                    }
                }
            }

            // go up
            peakIndex = baseIsotopePeakIndex;
            for (var isotopeIndex = baseIsotopeIndex + 1; isotopeIndex < isotopomerEnvelope.Length; isotopeIndex++)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex + 1; i < Peaks.Length; i++)
                {
                    var peakMz = Peaks[i].Mz;
                    if (peakMz > maxMz) return false;
                    if (peakMz >= minMz)    // find match, move to prev isotope
                    {
                        peakIndex = i;
                        break;
                    }
                }
            }

            return true;
        }

        /// <summary>
        /// Finds all isotope peaks corresponding to theoretical profiles with relative intensity higher than the threshold
        /// </summary>
        /// <param name="ion">ion</param>
        /// <param name="tolerance">tolerance</param>
        /// <param name="relativeIntensityThreshold">relative intensity threshold of the theoretical isotope profile</param>
        /// <returns>all isotope peaks found in the spectrum</returns>
        public IEnumerable<Peak> GetAllIonPeaks(Ion ion, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelop();
            var observedPeaks = new List<Peak>();

            var downHill = false;
            var prevTheoIntensity = 0f;
            // Assume that the isotopomer envelop is unimodal
            for (var isotopeIndex = 0; isotopeIndex <= isotopomerEnvelope.Length; isotopeIndex--)
            {
                var theoIntensity = isotopomerEnvelope[isotopeIndex];
                if (!downHill && theoIntensity < prevTheoIntensity) downHill = true;
                if (theoIntensity < relativeIntensityThreshold)
                {
                    if (downHill) break;
                }
                else
                {
                    var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                    var peak = FindPeak(isotopeMz, tolerance);
                    if (peak != null) observedPeaks.Add(peak);
                }
                prevTheoIntensity = theoIntensity;
            }

            return observedPeaks;
        }

        /// <summary>
        /// Computes the fit score between the ion and corresponding peaks in the spectrum
        /// </summary>
        /// <param name="ion">ion</param>
        /// <param name="tolerance">tolerance</param>
        /// <param name="relativeIntensityThreshold">relative intensity threshold of the theoretical isotope profile</param>
        /// <returns>fit score</returns>
        public double GetFitScore(Ion ion, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelop();
            var observedPeaks = new float[isotopomerEnvelope.Length];

            var downHill = false;
            var prevTheoIntensity = 0f;
            var maxObservedIntensity = float.NegativeInfinity;
            // Assume that the isotopomer envelop is unimodal
            for (var isotopeIndex = 0; isotopeIndex <= isotopomerEnvelope.Length; isotopeIndex--)
            {
                var theoIntensity = isotopomerEnvelope[isotopeIndex];
                if (!downHill && theoIntensity < prevTheoIntensity) downHill = true;
                if (theoIntensity < relativeIntensityThreshold)
                {
                    if (downHill) break;
                }
                else
                {
                    var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                    var peak = FindPeak(isotopeMz, tolerance);
                    if (peak != null)
                    {
                        var intensity = (float)peak.Intensity;
                        observedPeaks[isotopeIndex] = intensity;
                        if (intensity > maxObservedIntensity) maxObservedIntensity = intensity;
                    }
                    else observedPeaks[isotopeIndex] = 0f;
                }
                prevTheoIntensity = theoIntensity;
            }

            return FitScoreCalculator.GetFitOfNormalizedVectors(isotopomerEnvelope, observedPeaks);
        }

        // Added by Chris
        public IEnumerable<Peak> GetIonPeaks(Ion ion, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var baseIsotopeIndex = ion.Composition.GetMostAbundantIsotopeZeroBasedIndex();
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelop();
            var baseIsotopMz = ion.GetIsotopeMz(baseIsotopeIndex);
            var baseIsotopePeakIndex = FindPeakIndex(baseIsotopMz, tolerance);
            var finalPeaks = new List<Peak>();
            if (baseIsotopePeakIndex < 0) return null;

            // go down
            var peakIndex = baseIsotopePeakIndex;
            for (var isotopeIndex = baseIsotopeIndex - 1; isotopeIndex >= 0; isotopeIndex--)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex - 1; i >= 0; i--)
                {
                    var peakMz = Peaks[i].Mz;
                    if (peakMz < minMz) break;
                    if (peakMz <= maxMz)    // find match, move to prev isotope
                    {
                        finalPeaks.Add(Peaks[i]);
                        peakIndex = i;
                        break;
                    }
                }
            }

            // go up
            peakIndex = baseIsotopePeakIndex;
            for (var isotopeIndex = baseIsotopeIndex + 1; isotopeIndex < isotopomerEnvelope.Length; isotopeIndex++)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex + 1; i < Peaks.Length; i++)
                {
                    var peakMz = Peaks[i].Mz;
                    if (peakMz > maxMz) return null;
                    if (peakMz >= minMz)    // find match, move to prev isotope
                    {
                        finalPeaks.Add(Peaks[i]);
                        peakIndex = i;
                        break;
                    }
                }
            }

            return finalPeaks;
        }


        /// <summary>
        /// Finds the maximum intensity peak within the specified range
        /// </summary>
        /// <param name="minMz">minimum m/z</param>
        /// <param name="maxMz">maximum m/z</param>
        /// <returns>maximum intensity peak within the range</returns>
        public Peak FindPeak(double minMz, double maxMz)
        {
            var peakIndex = FindPeakIndex(minMz, maxMz);
            if (peakIndex < 0) return null;
            return Peaks[peakIndex];
        }

        public void Display()
        {
            var sb = new StringBuilder();
            sb.Append("--------- Spectrum -----------------\n");
            foreach (var peak in Peaks)
            {
                sb.Append(peak.Mz);
                sb.Append("\t");
                sb.Append(peak.Intensity);
                sb.Append("\n");
            }
            sb.Append("--------------------------- end ---------------------------------------\n");

            Console.Write(sb.ToString());
        }

        public void FilterNoise(double signalToNoiseRatio = 1.4826)
        {
            if (Peaks.Length < 2) return;
            Array.Sort(Peaks, new IntensityComparer());
            var medianIntPeak = Peaks[Peaks.Length / 2];
            var noiseLevel = medianIntPeak.Intensity;

            var filteredPeaks = Peaks.TakeWhile(peak => !(peak.Intensity < noiseLevel * signalToNoiseRatio)).ToList();

            filteredPeaks.Sort();
            Peaks = filteredPeaks.ToArray();
        }

        private int _msLevel = 1;

        private int FindPeakIndex(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return FindPeakIndex(minMz, maxMz);
        }

        private int FindPeakIndex(double minMz, double maxMz)
        {
            var index = Array.BinarySearch(Peaks, new Peak((minMz + maxMz) / 2, 0));
            if (index < 0) index = ~index;

            var bestPeakIndex = -1;
            var bestIntensity = 0.0;

            // go down
            var i = index - 1;
            while (i >= 0 && i < Peaks.Length)
            {
                if (Peaks[i].Mz <= minMz) break;
                if (Peaks[i].Intensity > bestIntensity)
                {
                    bestIntensity = Peaks[i].Intensity;
                    bestPeakIndex = i;
                }
                --i;
            }

            // go up
            i = index;
            while (i >= 0 && i < Peaks.Length)
            {
                if (Peaks[i].Mz >= maxMz) break;
                if (Peaks[i].Intensity > bestIntensity)
                {
                    bestIntensity = Peaks[i].Intensity;
                    bestPeakIndex = i;
                }
                ++i;
            }
            return bestPeakIndex;
        }
    }
}
