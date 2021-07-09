using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MathAndStats;
using MathNet.Numerics.Statistics;

// ReSharper disable UnusedMember.Global

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Peak List utilities
    /// </summary>
    public static class PeakListUtils
    {
        /// <summary>
        /// Find the highest-intensity peak that matches <paramref name="mz"/> within a tolerance
        /// </summary>
        /// <param name="peakList"></param>
        /// <param name="mz"></param>
        /// <param name="tolerance"></param>
        /// <returns>Peak in the given m/z range with the highest intensity</returns>
        public static Peak FindPeak(List<Peak> peakList, double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsMz(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return FindPeak(peakList, minMz, maxMz);
        }

        /// <summary>
        /// Find the highest-intensity peak with the m/z range specified
        /// </summary>
        /// <param name="peakList"></param>
        /// <param name="minMz"></param>
        /// <param name="maxMz"></param>
        /// <returns>Peak in the given m/z range with the highest intensity</returns>
        public static Peak FindPeak(List<Peak> peakList, double minMz, double maxMz)
        {
            //var index = peakList.BinarySearch(new Peak(mz, 0.0), comparer);
            //return index < 0 ? null : peakList[index];
            var index = peakList.BinarySearch(new Peak((minMz + maxMz) / 2, 0));
            if (index < 0)
            {
                index = ~index;
            }

            var bestPeakIndex = -1;
            var bestIntensity = 0.0;

            // go down
            var i = index - 1;
            while (i >= 0 && i < peakList.Count)
            {
                if (peakList[i].Mz <= minMz)
                {
                    break;
                }

                if (peakList[i].Intensity > bestIntensity)
                {
                    bestIntensity = peakList[i].Intensity;
                    bestPeakIndex = i;
                }
                --i;
            }

            // go up
            i = index;
            while (i >= 0 && i < peakList.Count)
            {
                if (peakList[i].Mz >= maxMz)
                {
                    break;
                }

                if (peakList[i].Intensity > bestIntensity)
                {
                    bestIntensity = peakList[i].Intensity;
                    bestPeakIndex = i;
                }
                ++i;
            }
            return bestPeakIndex == -1 ? null : peakList[bestPeakIndex];
        }

        /// <summary>
        /// Find all peaks that are in the range specified by <paramref name="mz"/> and a tolerance
        /// </summary>
        /// <param name="peakList"></param>
        /// <param name="mz"></param>
        /// <param name="tolerance"></param>
        /// <returns>List of matched peaks</returns>
        public static IList<Peak> FindAllPeaks(List<Peak> peakList, double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsMz(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return FindAllPeaks(peakList, minMz, maxMz);
        }

        /// <summary>
        /// Find all peaks that are in the range provided
        /// </summary>
        /// <param name="peakList"></param>
        /// <param name="minMz"></param>
        /// <param name="maxMz"></param>
        /// <returns>List of matched peaks</returns>
        public static IList<Peak> FindAllPeaks(List<Peak> peakList, double minMz, double maxMz)
        {
            //var index = peakList.BinarySearch(new Peak(mz, 0.0), comparer);
            //return index < 0 ? null : peakList[index];
            var index = peakList.BinarySearch(new Peak((minMz + maxMz) / 2, 0));
            if (index < 0)
            {
                index = ~index;
            }

            var matchedPeakList = new List<Peak>();
            // go down
            var i = index - 1;
            while (i >= 0 && i < peakList.Count)
            {
                if (peakList[i].Mz <= minMz)
                {
                    break;
                }

                matchedPeakList.Add(peakList[i]);
                --i;
            }

            // go up
            i = index;
            while (i >= 0 && i < peakList.Count)
            {
                if (peakList[i].Mz >= maxMz)
                {
                    break;
                }

                matchedPeakList.Add(peakList[i]);
                ++i;
            }
            matchedPeakList.Sort();
            return matchedPeakList;
        }

        /// <summary>
        /// Get the Pearson correlation for the provided experimental and theoretical peaks
        /// </summary>
        /// <param name="spec"></param>
        /// <param name="isoProfile"></param>
        /// <param name="comparer"></param>
        /// <returns>Pearson correlation</returns>
        public static double GetPearsonCorrelation(IList<Peak> spec, IList<Peak> isoProfile,
            IComparer<Peak> comparer)
        {
            var numPeaksSpec = spec.Count;
            var numTheoProfilePeaks = isoProfile.Count;
            var index1 = 0;
            var index2 = 0;

            var theoIntensities = isoProfile.Select(p => p.Intensity).ToArray();
            var observedIntensities = new double[numTheoProfilePeaks];

            while (index1 < numPeaksSpec && index2 < numTheoProfilePeaks)
            {
                var comp = comparer.Compare(spec[index1], isoProfile[index2]);
                if (comp < 0)
                {
                    ++index1;
                }
                else if (comp > 0)
                {
                    ++index2;
                }
                else
                {
                    if (spec[index1].Intensity > observedIntensities[index2])
                    {
                        observedIntensities[index2] = spec[index1].Intensity;
                    }
                    ++index1;
                    //++index2;
                }
            }

            return FitScoreCalculator.GetPearsonCorrelation(theoIntensities, observedIntensities);
        }

        /// <summary>
        /// Get the cosine score for the provided experimental and theoretical peaks
        /// </summary>
        /// <param name="spec"></param>
        /// <param name="isoProfile"></param>
        /// <param name="comparer"></param>
        /// <returns>Cosine score</returns>
        public static double GetCosine(IList<Peak> spec, IList<Peak> isoProfile, IComparer<Peak> comparer)
        {
            var numPeaksSpec = spec.Count;
            var numTheoProfilePeaks = isoProfile.Count;
            var index1 = 0;
            var index2 = 0;

            var theoIntensities = isoProfile.Select(p => p.Intensity).ToArray();
            var observedIntensities = new double[numTheoProfilePeaks];

            while (index1 < numPeaksSpec && index2 < numTheoProfilePeaks)
            {
                var comp = comparer.Compare(spec[index1], isoProfile[index2]);
                if (comp < 0)
                {
                    ++index1;
                }
                else if (comp > 0)
                {
                    ++index2;
                }
                else
                {
                    if (spec[index1].Intensity > observedIntensities[index2])
                    {
                        observedIntensities[index2] = spec[index1].Intensity;
                    }
                    ++index1;
                    //++index2;
                }
            }

            return FitScoreCalculator.GetCosine(theoIntensities, observedIntensities);
        }

        /// <summary>
        /// Get the intersection of the provided peak lists
        /// </summary>
        /// <param name="peakList1"></param>
        /// <param name="peakList2"></param>
        /// <param name="comparer"></param>
        /// <returns>List of matching peaks</returns>
        public static List<Peak> GetIntersection(List<Peak> peakList1, List<Peak> peakList2,
            IComparer<Peak> comparer)
        {
            var count1 = peakList1.Count;
            var count2 = peakList2.Count;
            var index1 = 0;
            var index2 = 0;

            var intersection = new List<Peak>();
            while (index1 < count1 && index2 < count2)
            {
                var comp = comparer.Compare(peakList1[index1], peakList2[index2]);
                if (comp < 0)
                {
                    ++index1;
                }
                else if (comp > 0)
                {
                    ++index2;
                }
                else
                {
                    intersection.Add(peakList1[index1]);
                    ++index1;
                    ++index2;
                }
            }

            return intersection;
        }

        /// <summary>
        /// Compare two peak lists and create a list of peaks that are in peakList1 but not in peakList2
        /// </summary>
        /// <param name="peakList1"></param>
        /// <param name="peakList2"></param>
        /// <param name="comparer"></param>
        /// <returns>List of peaks only in peakList1</returns>
        public static List<Peak> GetExceptWith(List<Peak> peakList1, List<Peak> peakList2,
            IComparer<Peak> comparer)
        {
            var count1 = peakList1.Count;
            var count2 = peakList2.Count;
            var index1 = 0;
            var index2 = 0;

            var exceptWith = new List<Peak>();
            while (index1 < count1 && index2 < count2)
            {
                var comp = comparer.Compare(peakList1[index1], peakList2[index2]);
                if (comp < 0)
                {
                    exceptWith.Add(peakList1[index1]);
                    ++index1;
                }
                else if (comp > 0)
                {
                    ++index2;
                }
                else
                {
                    ++index1;
                    //++index2;
                }
            }

            return exceptWith;
        }

        /// <summary>
        /// Get the sum of the two peak lists
        /// </summary>
        /// <param name="peakList1"></param>
        /// <param name="peakList2"></param>
        /// <param name="comparer"></param>
        /// <returns>Combined list of peaks (without any duplicates)</returns>
        public static List<Peak> Sum(IList<Peak> peakList1, IList<Peak> peakList2,
            IComparer<Peak> comparer)
        {
            var count1 = peakList1.Count;
            var count2 = peakList2.Count;
            var index1 = 0;
            var index2 = 0;

            var sum = new List<Peak>();

            while (index1 < count1 && index2 < count2)
            {
                var p1 = peakList1[index1];
                var p2 = peakList2[index2];

                var comp = comparer.Compare(peakList1[index1], peakList2[index2]);
                if (comp < 0)
                {
                    sum.Add(p1);
                    ++index1;
                }
                else if (comp > 0)
                {
                    sum.Add(p2);
                    ++index2;
                }
                else
                {
                    peakList1[index1] = new Peak(p1.Mz, p1.Intensity + p2.Intensity);
                    ++index2;
                }
            }

            while (index1 < count1)
            {
                sum.Add(peakList1[index1++]);
            }

            while (index2 < count2)
            {
                sum.Add(peakList2[index2++]);
            }

            return sum;
        }

        /// <summary>
        /// Get a peak list where peaks before a threshold are removed
        /// </summary>
        /// <param name="peakList"></param>
        /// <param name="filteredPeakList"></param>
        /// <param name="signalToMedianRatio"></param>
        public static void FilterNoise(IList<Peak> peakList,
            ref List<Peak> filteredPeakList, double signalToMedianRatio = 1.4826)
        {
            //var medianIntensity = peakList.OrderByDescending(p => p.Intensity).Select(p => p.Intensity).Median();
            var medianIntensity = peakList.Select(p => p.Intensity).Median();
            filteredPeakList.AddRange(peakList.OrderByDescending(p => p.Intensity).TakeWhile(p => !(p.Intensity < medianIntensity * signalToMedianRatio)).OrderBy(p => p.Mz));
        }
    }
}
