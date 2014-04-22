using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class SortedPeakSet
    {


        //public static int FindPeakIndex(List<Peak> peakList, double mz, Tolerance tolerance)
        //{
        //    var tolTh = tolerance.GetToleranceAsTh(mz);
        //    var minMz = mz - tolTh;
        //    var maxMz = mz + tolTh;
        //    return FindPeakIndex(peakList, minMz, maxMz);
        //}

        //public static int FindPeakIndex(List<Peak> peakList, double minMz, double maxMz)
        //{
        //    var index = peakList.BinarySearch(new Peak((minMz + maxMz) / 2, 0));
        //    if (index < 0) index = ~index;

        //    var bestPeakIndex = -1;
        //    var bestIntensity = 0.0;

        //    // go down
        //    var i = index - 1;
        //    while (i >= 0 && i < peakList.Count)
        //    {
        //        if (peakList[i].Mz <= minMz) break;
        //        if (peakList[i].Intensity > bestIntensity)
        //        {
        //            bestIntensity = peakList[i].Intensity;
        //            bestPeakIndex = i;
        //        }
        //        --i;
        //    }

        //    // go up
        //    i = index;
        //    while (i >= 0 && i < peakList.Count)
        //    {
        //        if (peakList[i].Mz >= maxMz) break;
        //        if (peakList[i].Intensity > bestIntensity)
        //        {
        //            bestIntensity = peakList[i].Intensity;
        //            bestPeakIndex = i;
        //        }
        //        ++i;
        //    }
        //    return bestPeakIndex;
        //}
    }
}
