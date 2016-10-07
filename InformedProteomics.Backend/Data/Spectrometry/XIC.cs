using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Xic: List<XicPoint>
    {
        //public int Min { get; private set; }
        //public int Max { get; private set; }

        private static readonly SavitzkyGolaySmoother Smoother;
        static Xic()
        {
            Smoother = new SavitzkyGolaySmoother(9, 2);
        }

        public double GetCorrelation(Xic other)
        {
            if (Count == 0 || other == null || other.Count == 0) return 0;

            var count1 = Count;
            var count2 = other.Count;
            var index1 = 0;
            var index2 = 0;

            var intList1 = new List<double>();
            var intList2 = new List<double>();
            while (index1 < count1 && index2 < count2)
            {
                var comp = this[index1].ScanNum - other[index2].ScanNum;
                if (comp < 0) ++index1;
                else if (comp > 0) ++index2;
                else
                {
                    intList1.Add(this[index1].Intensity);
                    intList2.Add(other[index2].Intensity);
                    ++index1;
                    ++index2;
                }
            }

            var correlation = FitScoreCalculator.GetPearsonCorrelation(intList1.ToArray(), intList2.ToArray());
            return correlation;
        }

        public double GetCosine(Xic other)
        {
            if (Count == 0 || other == null || other.Count == 0) return 0;

            //var minScanNum = Math.Min(this[0].ScanNum, other[0].ScanNum);
            //var maxScanNum = Math.Max(this[Count-1].ScanNum, other[other.Count-1].ScanNum);

            var minScanNum = this[0].ScanNum;
            var maxScanNum = this[Count - 1].ScanNum;

            var intArr1 = new double[maxScanNum - minScanNum + 1];
            foreach (var p in this)
            {
                intArr1[p.ScanNum - minScanNum] = p.Intensity;
            }

            var intArr2 = new double[maxScanNum - minScanNum + 1];
            foreach (var p in other)
            {
                var index = p.ScanNum - minScanNum;
                if (index < 0 || index >= intArr2.Length) continue;
                intArr2[p.ScanNum - minScanNum] = p.Intensity;
            }

            var correlation = FitScoreCalculator.GetCosine(intArr1, intArr2);
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
            if (Count < 6) return GetApexScanNum();

            var xicPointList = new List<XicPoint>();

            if(performSmoothing)
            {
                double[] intensityValues = this.Select(x => x.Intensity).ToArray();
                intensityValues = Smoother.Smooth(intensityValues);

                for(var i = 0; i < Count; i++)
                {
                    xicPointList.Add(new XicPoint(this[i].ScanNum, this[i].Mz, intensityValues[i]));
                }
            }
            else
            {
                xicPointList = this;
            }

            // Find the XIC Point that is closest to the input scan number
            var searchPoint = new XicPoint(scanNumber, 0, 0);
            int indexOfClosestScan = xicPointList.BinarySearch(searchPoint, new AnonymousComparer<XicPoint>((x, y) => x.ScanNum.CompareTo(y.ScanNum)));
            indexOfClosestScan = indexOfClosestScan < 0 ? ~indexOfClosestScan : indexOfClosestScan;
            if (indexOfClosestScan >= xicPointList.Count) indexOfClosestScan = xicPointList.Count - 1;
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

        /// <summary>
        /// Display the chromatogram
        /// </summary>
        /// <param name="maxPointsToShow">Maximum number of data points to show</param>
        /// <remarks>Set maxPoints to 0 to see all of the data points</remarks>
        public void Display(int maxPointsToShow = 50)
        {
            var threshold = 1;
            if (maxPointsToShow > 0)
            {
                threshold = (int)Math.Round(this.Count / (double)maxPointsToShow);
                if (threshold < 1)
                    threshold = 1;
            }

            var index = 0;
            var pointsShown = 0;

            foreach (var p in this)
            {
                if (index == 0 || index % threshold == 0)
                {
                    Console.WriteLine(p.ScanNum + "\t" + p.Intensity);
                    pointsShown++;
                }

                index++;
            }

            Console.WriteLine("Displayed {0} out of {1} data points", pointsShown, this.Count);
        }

        // sort XicPoints and select one peak per scan
        public static Xic GetSelectedXic(Xic xic)
        {
            if (xic.Count == 0) return xic;
            xic.Sort(); // Need to guarantee a sorted Xic

            // select one best peak for each scan
            var newXic = new Xic();

            var prevScanNum = xic[0].ScanNum;
            var bestPeak = xic[0];
            for (var i = 1; i < xic.Count; i++)
            {
                var xicPeak = xic[i];
                if (xicPeak.ScanNum > prevScanNum)
                {
                    newXic.Add(bestPeak);
                    bestPeak = xicPeak;
                }
                else
                {
                    if (xicPeak.Intensity > bestPeak.Intensity)
                    {
                        bestPeak = xicPeak;
                    }
                }
                prevScanNum = xicPeak.ScanNum;
            }
            newXic.Add(bestPeak);
            return newXic;
        }

        public override bool Equals(object obj)
        {
            return obj as Xic != null && Equals((Xic) obj);
        }

        protected bool Equals(Xic other)
        {
            for (var i = 0; i < Count; i++)
            {
                if (!this[i].Equals(other[i])) return false;
            }
            return true;
        }

        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
    }
}
