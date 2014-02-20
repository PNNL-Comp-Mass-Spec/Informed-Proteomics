using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Utils;
using MultiDimensionalPeakFinding;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Xic: List<XicPoint>
    {
        //public int MinScanNum { get; private set; }
        //public int MaxScanNum { get; private set; }

    	private static SavitzkyGolaySmoother _smoother;

		static Xic()
		{
			_smoother = new SavitzkyGolaySmoother(9, 2);
		}

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

        public int GetApexScanNum()
        {
            var maxIntensity = double.MinValue;
            var apexScan = -1;
            foreach (var p in this)
            {
                if (p.Intensity > maxIntensity)
                {
                    apexScan = p.ScanNum;
                    maxIntensity = p.Intensity;
                }
            }
            return apexScan;
        }

		public int GetNearestApexScanNum(int scanNumber, bool performSmoothing = true)
		{
			// If there are not very many points, just return the global apex
			if (this.Count < 6) return GetApexScanNum();

			List<XicPoint> xicPointList = new List<XicPoint>();

			if(performSmoothing)
			{
				double[] intensityValues = this.Select(x => x.Intensity).ToArray();
				intensityValues = _smoother.Smooth(intensityValues);

				for(int i = 0; i < this.Count; i++)
				{
					xicPointList.Add(new XicPoint(this[i].ScanNum, intensityValues[i]));
				}
			}
			else
			{
				xicPointList = this;
			}

			// Find the XIC Point that is closest to the input scan number
			XicPoint searchPoint = new XicPoint(scanNumber, 0);
			int indexOfClosestScan = xicPointList.BinarySearch(searchPoint, new AnonymousComparer<XicPoint>((x, y) => x.ScanNum.CompareTo(y.ScanNum)));
			indexOfClosestScan = indexOfClosestScan < 0 ? ~indexOfClosestScan : indexOfClosestScan;
			XicPoint closestXicPoint = xicPointList[indexOfClosestScan];

			// Figure out if we want to search for an apex by moving left or right
			bool moveRight;
			if (indexOfClosestScan <= 1) moveRight = true;
			else if (indexOfClosestScan >= xicPointList.Count - 2) moveRight = false;
			else if (xicPointList[indexOfClosestScan + 1].Intensity > closestXicPoint.Intensity) moveRight = true;
			else moveRight = false;

			// Check to the right
			if(moveRight)
			{
				if (indexOfClosestScan + 1 >= xicPointList.Count) return GetApexScanNum();
				double previousIntensity = xicPointList[indexOfClosestScan + 1].Intensity;

				for (int i = indexOfClosestScan + 2; i < xicPointList.Count; i++)
				{
					double currentIntensity = xicPointList[i].Intensity;
					if (currentIntensity < previousIntensity) return xicPointList[i-1].ScanNum;
					previousIntensity = currentIntensity;
				}
			}
			// Check to the left
			else
			{
				if (indexOfClosestScan - 1 < 0) return GetApexScanNum();
				double previousIntensity = this[indexOfClosestScan - 1].Intensity;

				for (int i = indexOfClosestScan - 2; i >= 0; i--)
				{
					double currentIntensity = this[i].Intensity;
					if (currentIntensity < previousIntensity) return this[i+1].ScanNum;
					previousIntensity = currentIntensity;
				}
			}

			// I don't think it is possible, but if we make it this far, then we should just return the apex of the whole XIC because a single peak was not discovered
			return GetApexScanNum();
		}
    }
}
