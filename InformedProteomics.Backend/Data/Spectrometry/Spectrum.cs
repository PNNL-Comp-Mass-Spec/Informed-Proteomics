using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Utils;

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

        public double ElutionTime { get; set; }

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
        /// Gets a list of peaks within [minMz, maxMz] 
        /// </summary>
        /// <param name="minMz">minimum m/z</param>
        /// <param name="maxMz">maximum m/z</param>
        /// <returns>list of peaks within [minMz, maxMz]</returns>
        public List<Peak> GetPeakListWithin(double minMz, double maxMz)
        {
            var peakList = new List<Peak>();

            GetPeakListWithin(minMz, maxMz, ref peakList);

            return peakList;
        }

        /// <summary>
        /// Gets a list of peaks within [minMz, maxMz] and add to peakList
        /// </summary>
        /// <param name="minMz">minimum m/z</param>
        /// <param name="maxMz">maximum m/z</param>
        /// <param name="peakList">list of peaks where the peaks will be added</param>
        /// <returns>list of peaks within [minMz, maxMz]</returns>
        public void GetPeakListWithin(double minMz, double maxMz, ref List<Peak> peakList)
        {
            var index = Array.BinarySearch(Peaks, new Peak(minMz - float.Epsilon, 0));
            if (index < 0) index = ~index;

            // go up
            var i = index;
            while (i < Peaks.Length)
            {
                if (Peaks[i].Mz > maxMz) break;
                peakList.Add(Peaks[i]);
                ++i;
            }
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
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelope();
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
        /// <returns>array of observed isotope peaks in the spectrum. null if no peak found.</returns>
        public Peak[] GetAllIsotopePeaks(Ion ion, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var mostAbundantIsotopeIndex = ion.Composition.GetMostAbundantIsotopeZeroBasedIndex();
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelope();
            var baseIsotopMz = ion.GetIsotopeMz(mostAbundantIsotopeIndex);
            var mostAbundantIsotopeMatchedPeakIndex = FindPeakIndex(baseIsotopMz, tolerance);
            if (mostAbundantIsotopeMatchedPeakIndex < 0) return null;

            var observedPeaks = new Peak[isotopomerEnvelope.Length];
            observedPeaks[mostAbundantIsotopeIndex] = Peaks[mostAbundantIsotopeMatchedPeakIndex];

            // go down
            var peakIndex = mostAbundantIsotopeMatchedPeakIndex - 1;
            for (var isotopeIndex = mostAbundantIsotopeIndex - 1; isotopeIndex >= 0; isotopeIndex--)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex; i >= 0; i--)
                {
                    var peakMz = Peaks[i].Mz;
                    if (peakMz < minMz)
                    {
                        peakIndex = i;
                        break;
                    }
                    if (peakMz <= maxMz)    // find match, move to prev isotope
                    {
                        var peak = Peaks[i];
                        if (observedPeaks[isotopeIndex] == null ||
                            peak.Intensity > observedPeaks[isotopeIndex].Intensity)
                        {
                            observedPeaks[isotopeIndex] = peak;
                        }
                    }
                }
            }

            // go up
            peakIndex = mostAbundantIsotopeMatchedPeakIndex + 1;
            for (var isotopeIndex = mostAbundantIsotopeIndex + 1; isotopeIndex < isotopomerEnvelope.Length; isotopeIndex++)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex; i < Peaks.Length; i++)
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
                        if (observedPeaks[isotopeIndex] == null ||
                            peak.Intensity > observedPeaks[isotopeIndex].Intensity)
                        {
                            observedPeaks[isotopeIndex] = peak;
                        }
                    }
                }
            }

            return observedPeaks;
        }

        public Peak[] GetAllIsotopePeaks(double monoIsotopeMass, int charge, IsotopomerEnvelope envelope, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var mostAbundantIsotopeIndex = envelope.MostAbundantIsotopeIndex;
            var isotopomerEnvelope = envelope.Envolope;
            var mostAbundantIsotopeMz = Ion.GetIsotopeMz(monoIsotopeMass, charge, mostAbundantIsotopeIndex);
            var mostAbundantIsotopePeakIndex = FindPeakIndex(mostAbundantIsotopeMz, tolerance);
            if (mostAbundantIsotopePeakIndex < 0) return null;

            var observedPeaks = new Peak[isotopomerEnvelope.Length];
            observedPeaks[mostAbundantIsotopeIndex] = Peaks[mostAbundantIsotopePeakIndex];

            // go down
            var peakIndex = mostAbundantIsotopePeakIndex - 1;
            for (var isotopeIndex = mostAbundantIsotopeIndex - 1; isotopeIndex >= 0; isotopeIndex--)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = Ion.GetIsotopeMz(monoIsotopeMass, charge, isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex; i >= 0; i--)
                {
                    var peakMz = Peaks[i].Mz;
                    if (peakMz < minMz)
                    {
                        peakIndex = i;
                        break;
                    }
                    if (peakMz <= maxMz)    // find match, move to prev isotope
                    {
                        var peak = Peaks[i];
                        if (observedPeaks[isotopeIndex] == null ||
                            peak.Intensity > observedPeaks[isotopeIndex].Intensity)
                        {
                            observedPeaks[isotopeIndex] = peak;
                        }
                    }
                }
            }

            // go up
            peakIndex = mostAbundantIsotopePeakIndex + 1;
            for (var isotopeIndex = mostAbundantIsotopeIndex + 1; isotopeIndex < isotopomerEnvelope.Length; isotopeIndex++)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = Ion.GetIsotopeMz(monoIsotopeMass, charge, isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex; i < Peaks.Length; i++)
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
                        if (observedPeaks[isotopeIndex] == null ||
                            peak.Intensity > observedPeaks[isotopeIndex].Intensity)
                        {
                            observedPeaks[isotopeIndex] = peak;
                        }
                    }
                }
            }
            return observedPeaks;
        }

        /// <summary>
        /// Computes the Pearson correlation between the ion and corresponding peaks in the spectrum
        /// </summary>
        /// <param name="ion">ion</param>
        /// <param name="tolerance">tolerance</param>
        /// <param name="relativeIntensityThreshold">relative intensity threshold of the theoretical isotope profile</param>
        /// <returns>Pearson correlation</returns>
        public double GetCorrScore(Ion ion, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var observedPeaks = GetAllIsotopePeaks(ion, tolerance, relativeIntensityThreshold);
            if (observedPeaks == null) return 0;

            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelope();
            var observedIntensities = new double[observedPeaks.Length];

            for (var i = 0; i < observedPeaks.Length; i++)
            {
                var observedPeak = observedPeaks[i];
                observedIntensities[i] = observedPeak != null ? (float)observedPeak.Intensity : 0.0;
            }
            return FitScoreCalculator.GetPearsonCorrelation(isotopomerEnvelope, observedIntensities);
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
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelope();
            var observedPeaks = GetAllIsotopePeaks(ion, tolerance, relativeIntensityThreshold);

            if (observedPeaks == null) return 1;
            var theoIntensities = new double[observedPeaks.Length];
            Array.Copy(isotopomerEnvelope, theoIntensities, theoIntensities.Length);

            var maxObservedIntensity = observedPeaks.Select(p => p != null ? p.Intensity : 0).Max();
            var normalizedObs = observedPeaks.Select(p => p != null ? p.Intensity / maxObservedIntensity : 0).ToArray();
            return FitScoreCalculator.GetFitOfNormalizedVectors(isotopomerEnvelope, normalizedObs);
        }

        /// <summary>
        /// Computes the cosine between the ion and corresponding peaks in the spectrum
        /// </summary>
        /// <param name="ion">ion</param>
        /// <param name="tolerance">tolerance</param>
        /// <param name="relativeIntensityThreshold">relative intensity threshold of the theoretical isotope profile</param>
        /// <returns>cosine value</returns>
        public double GetConsineScore(Ion ion, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var observedPeaks = GetAllIsotopePeaks(ion, tolerance, relativeIntensityThreshold);
            if (observedPeaks == null) return 0;

            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelope();
            var theoIntensities = new double[observedPeaks.Length];
            var observedIntensities = new double[observedPeaks.Length];

            for (var i = 0; i < observedPeaks.Length; i++)
            {
                theoIntensities[i] = isotopomerEnvelope[i];
                var observedPeak = observedPeaks[i];
                observedIntensities[i] = observedPeak != null ? (float)observedPeak.Intensity : 0.0;
            }
            return FitScoreCalculator.GetCosine(isotopomerEnvelope, observedIntensities);
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

        public int FindPeakIndex(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return FindPeakIndex(minMz, maxMz);
        }

        public int FindPeakIndex(double minMz, double maxMz)
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
