using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Xic: List<XicPeak>
    {
        public double GetCorrelation(Xic other)
        {
            if (this.Count == 0 || other == null || other.Count == 0) return 0;

            var minScanNum = Math.Min(this[0].ScanNum, other[0].ScanNum);
            var maxScanNum = Math.Max(this[Count-1].ScanNum, other[other.Count-1].ScanNum);

            var intArr1 = new double[maxScanNum-minScanNum+1];
            foreach (var p in this) intArr1[p.ScanNum - minScanNum] = p.Intensity;

            var intArr2 = new double[maxScanNum - minScanNum + 1];
            foreach (var p in other) intArr2[p.ScanNum - minScanNum] = p.Intensity;

            var correlation = SimpleMath.GetCorrelation(intArr1, intArr2);
            return correlation;
        }

        public double GetSumIntensities()
        {
            return this.Sum(p => p.Intensity);
        }

        public bool ContainsScanNum(int scanNum)
        {
            return this.Any(xicPeak => xicPeak.ScanNum == scanNum);
        }
    }
}
