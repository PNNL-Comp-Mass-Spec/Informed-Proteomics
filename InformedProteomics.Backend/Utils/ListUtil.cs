using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Utils
{
    public static class ListUtil
    {
        public static IList<int> GetIntersection(List<int> intList1, List<int> intList2,
             IComparer<Peak> comparer)
        {
            var count1 = intList1.Count;
            var count2 = intList2.Count;
            var index1 = 0;
            var index2 = 0;

            var intersection = new List<int>();
            while (index1 < count1 && index2 < count2)
            {
                var comp = intList1[index1] - intList2[index2];
                if (comp < 0) ++index1;
                else if (comp > 0) ++index2;
                else intersection.Add(intList1[index1]);
            }

            return intersection;
        }

        public static void AddSorted<T>(this List<T> @this, T item) where T : IComparable<T>
        {
            if (@this.Count == 0)
            {
                @this.Add(item);
                return;
            }
            if (@this[@this.Count - 1].CompareTo(item) <= 0)
            {
                @this.Add(item);
                return;
            }
            if (@this[0].CompareTo(item) >= 0)
            {
                @this.Insert(0, item);
                return;
            }
            int index = @this.BinarySearch(item);
            if (index < 0)
                index = ~index;
            @this.Insert(index, item);
        }
    }
}
