using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using MultiDimensionalPeakFinding.PeakDetection;
using Point = System.Drawing.Point;

namespace InformedProteomics.Backend.IMS
{
    public class Feature
    {
        public ushort ScanLcStart { get; private set; }// 100
        public byte ScanLcLength { get; private set; }// 5 -> 100-104
        public byte ScanLcRepOffset { get; private set; }// to highest point 
        public ushort ScanImsStart { get; private set; }
        public byte ScanImsLength { get; private set; }
        public byte ScanImsRepOffset { get; private set; }
        public float IntensityMax { get; private set; }
        public float SumIntensities { get; private set; }
        public ushort NumPoints { get; private set; }
        public float[] LcApexPeakProfile { get; private set; }//ScanLcLength
        public float[] ImsApexPeakProfile { get; private set; }//ScanImsLength

        public Feature()
        {
        }

        public Feature(FeatureBlobStatistics featureBlobStatistics)
        {
            ScanLcStart = featureBlobStatistics.ScanLcStart;
            ScanLcLength = featureBlobStatistics.ScanLcLength;
            ScanLcRepOffset = featureBlobStatistics.ScanLcRepOffset;
            ScanImsStart = featureBlobStatistics.ScanImsStart;
            ScanImsLength = featureBlobStatistics.ScanImsLength;
            ScanImsRepOffset = featureBlobStatistics.ScanImsRepOffset;
            IntensityMax = featureBlobStatistics.IntensityMax;
            SumIntensities = featureBlobStatistics.SumIntensities;
            NumPoints = featureBlobStatistics.NumPoints;
            LcApexPeakProfile = featureBlobStatistics.LcApexPeakProfile;
            ImsApexPeakProfile = featureBlobStatistics.ImsApexPeakProfile;
        }

        public Rectangle GetBoundary()
        {
            return new Rectangle(ScanLcStart, ScanImsStart, ScanLcLength, ScanImsLength);
        }

        public Point GetHighestPoint()
        {
            return new Point(ScanLcStart+ScanLcRepOffset, ScanImsStart+ScanImsRepOffset);
        }

        public override string ToString()
        {
            return string.Format(
                "LC: [{0},{1}], IMS: [{2},{3}], Apex: [{4},{5}] SumIntensities: {6}, NumPoints: {7}",
                ScanLcStart,
                (ScanLcStart + ScanLcLength - 1),
                ScanImsStart,
                ScanImsStart + ScanImsLength - 1,
                ScanLcStart + ScanLcRepOffset,
                ScanImsStart + ScanImsRepOffset,
                SumIntensities,
                NumPoints
                );
        }
    }
}
