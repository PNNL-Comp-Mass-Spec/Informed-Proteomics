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
            ScanLcStart = (ushort) featureBlobStatistics.ScanLcMin;

            int scanLength = featureBlobStatistics.ScanLcMax - featureBlobStatistics.ScanLcMin + 1;
            if (scanLength > byte.MaxValue) ScanLcLength = byte.MaxValue;
            else ScanLcLength = (byte) (scanLength);

            int scanRepOffset = featureBlobStatistics.ScanLcRep - featureBlobStatistics.ScanLcMin;
            if (scanRepOffset > byte.MaxValue) ScanLcRepOffset = byte.MaxValue;
            else ScanLcRepOffset = (byte)(scanRepOffset);

            ScanImsStart = (ushort)featureBlobStatistics.ScanImsMin;
            int imsLength = featureBlobStatistics.ScanImsMax - featureBlobStatistics.ScanImsMin + 1;
            if (imsLength > byte.MaxValue) ScanImsLength = byte.MaxValue;
            else ScanImsLength = (byte)(imsLength);

            int imsRepOffset = featureBlobStatistics.ScanImsRep - featureBlobStatistics.ScanImsMin;
            if (imsRepOffset > byte.MaxValue) ScanImsRepOffset = byte.MaxValue;
            else ScanImsRepOffset = (byte)(imsRepOffset);

            IntensityMax = (float) featureBlobStatistics.IntensityMax;
            SumIntensities = (float) featureBlobStatistics.SumIntensities;
            NumPoints = (ushort) featureBlobStatistics.NumPoints;

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

        public bool IsEmpty() // added by kyowon
        {
            return NumPoints == 0;
        }
    }
}
