using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Utils
{
    public class PeakListUtils
    {
        public static Peak FindPeak(List<Peak> peakList, double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return FindPeak(peakList, minMz, maxMz);
        }

        public static Peak FindPeak(List<Peak> peakList, double minMz, double maxMz)
        {
            //var index = peakList.BinarySearch(new Peak(mz, 0.0), comparer);
            //return index < 0 ? null : peakList[index];
            var index = peakList.BinarySearch(new Peak((minMz + maxMz) / 2, 0));
            if (index < 0) index = ~index;

            var bestPeakIndex = -1;
            var bestIntensity = 0.0;

            // go down
            var i = index - 1;
            while (i >= 0 && i < peakList.Count)
            {
                if (peakList[i].Mz <= minMz) break;
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
                if (peakList[i].Mz >= maxMz) break;
                if (peakList[i].Intensity > bestIntensity)
                {
                    bestIntensity = peakList[i].Intensity;
                    bestPeakIndex = i;
                }
                ++i;
            }
            return bestPeakIndex == -1 ? null : peakList[bestPeakIndex];
        }

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
                if (comp < 0) ++index1;
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

        public static double GetCosine(IList<Peak> spec, IList<Peak> isoProfile,IComparer<Peak> comparer)
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
                if (comp < 0) ++index1;
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
                var comp = peakList1[index1].CompareTo(peakList2[index2]);
                if (comp < 0) ++index1;
                else if (comp > 0) ++index2;
                else
                {
                    intersection.Add(peakList1[index1]);
                    ++index1;
                    ++index2;
                }
            }

            return intersection;
        }

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
                else if (comp > 0) ++index2;
                else
                {
                    ++index1;
                    //++index2;
                }
            }

            return exceptWith;
        }

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

                var comp = p1.CompareTo(p2);
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

            while (index1 < count1) sum.Add(peakList1[index1++]);
            while (index2 < count2) sum.Add(peakList2[index2++]);

            return sum;
        }
    }
}
