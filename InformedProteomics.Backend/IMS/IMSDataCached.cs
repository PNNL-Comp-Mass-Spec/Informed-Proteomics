using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using UIMFLibrary;

namespace InformedProteomics.Backend.IMS
{
    public class ImsDataCached : ImsData
    {
        private Dictionary<int, FeatureSet> _precursorFeatureSetMap;
        private int _numPrecursorFeatures;

        private Dictionary<int, FeatureSet> _fragmentFeatureSetMap;
        private int _numFragmentFeatures;

        public ImsDataCached(string filePath) : this(filePath, 400.0, 2000.0, 10.0, 2500.0, 
            new Tolerance(25, DataReader.ToleranceType.PPM), new Tolerance(25, DataReader.ToleranceType.PPM))
        {
        }

        public ImsDataCached(string filePath, double minPrecursorMz, double maxPrecursorMz, 
            double minFragmentMz, double maxFragmentMz, Tolerance precursorTolerance, Tolerance fragmentTolerance)
            : base(filePath)
        {
            MinPrecursorMz = minPrecursorMz;
            MaxPrecursorMz = maxPrecursorMz;
            MinFragmentMz = minFragmentMz;
            MaxFragmentMz = maxFragmentMz;
            PrecursorTolerance = precursorTolerance;
            FragmentTolerance = fragmentTolerance;
        }

        public double MinPrecursorMz { get; private set; }
        public double MaxPrecursorMz { get; private set; }
        public double MinFragmentMz { get; private set; }
        public double MaxFragmentMz { get; private set; }
        public Tolerance PrecursorTolerance { get; private set; }
        public Tolerance FragmentTolerance { get; private set; }

        public int GetNumberOfPrecursorFeatures()
        {
            return _numPrecursorFeatures;
        }

        public int GetNumberOfFragmentFeatures()
        {
            return _numFragmentFeatures;
        }

        public FeatureSet GetPrecursorFeatures(double precursorMz)
        {
            return GetFeatures(precursorMz, true);
        }

        /// <summary>
        /// Gets the precursor feature within the boundary of the precursorFeature
        /// </summary>
        /// <param name="precursorMz">m/z of the precursor</param>
        /// <param name="precursorFeature">precursorFeature defining the retention and drift time boundary</param>
        /// <returns>the precursor feature within the boundary of the precursorFeature</returns>
        public Feature GetPrecursorFeature(double precursorMz, Feature precursorFeature)
        {
            return GetFeature(precursorMz, precursorFeature, true);
        }

        public FeatureSet GetFragmentFeatures(double fragmentMz)
        {
            return GetFeatures(fragmentMz, false);
        }

        /// <summary>
        /// Extracts the fragment feature within the boundary of the precursorFeature 
        /// </summary>
        /// <param name="fragmentMz"></param>
        /// <param name="precursorFeature"></param>
        /// <returns>the fragment feature within the boundary of the precursorFeature</returns>
        public Feature GetFramentFeature(double fragmentMz, Feature precursorFeature)
        {
            return GetFeature(fragmentMz, precursorFeature, false);
        }

        public void CreatePrecursorFeatures()
        {
            // Precursor
            _precursorFeatureSetMap = new Dictionary<int, FeatureSet>();

            int totalNumPrecursorFeatures = 0;

            int minPrecursorTargetBin = GetBinFromMz(MinPrecursorMz);
            int maxPrecursorTargetBin = GetBinFromMz(MaxPrecursorMz);
            if (maxPrecursorTargetBin >= GetNumberOfBins())
                maxPrecursorTargetBin = GetNumberOfBins() - 1;

            for (int targetBin = minPrecursorTargetBin; targetBin <= maxPrecursorTargetBin; targetBin++)
            {
                double mz = GetMzFromBin(targetBin);
                FeatureSet featureSet = GetFeatures(mz, PrecursorTolerance, DataReader.FrameType.MS1);
                if (featureSet.GetFeatures().Any())
                    _precursorFeatureSetMap.Add(targetBin, featureSet);
                totalNumPrecursorFeatures += featureSet.GetFeatures().Count();
            }

            _numPrecursorFeatures = totalNumPrecursorFeatures;
        }

        public void CreateFragmentFeatures()
        {
            // Fragment
            _fragmentFeatureSetMap = new Dictionary<int, FeatureSet>();

            int totalNumFragmentFeatures = 0;

            const double minFragmentMz = 0.0;
            const double maxFragmentMz = 2500.0;

            int minFragmentTargetBin = GetBinFromMz(minFragmentMz);
            int maxFragmentTargetBin = GetBinFromMz(maxFragmentMz);
            if (maxFragmentTargetBin >= GetNumberOfBins())
                maxFragmentTargetBin = GetNumberOfBins() - 1;

            Console.WriteLine("Generating fragment features (MinMz: " + minFragmentMz + " MaxMz: " + maxFragmentMz + ")");
            for (int targetBin = minFragmentTargetBin; targetBin <= maxFragmentTargetBin; targetBin++)
            {
                double mz = GetMzFromBin(targetBin);
                FeatureSet featureSet = GetFeatures(mz, FragmentTolerance, DataReader.FrameType.MS2);
                if (featureSet.GetFeatures().Any())
                    _fragmentFeatureSetMap.Add(targetBin, featureSet);
                totalNumFragmentFeatures += featureSet.GetFeatures().Count();
            }

            _numFragmentFeatures = totalNumFragmentFeatures;
        }

        private FeatureSet GetFeatures(double mz, bool isPrecursor)
        {
            int mzBin = GetBinFromMz(mz);
            if (isPrecursor)
            {
                return _precursorFeatureSetMap[mzBin];
            }
            return _fragmentFeatureSetMap[mzBin];
        }

        private Feature GetFeature(double mz, Feature precursorFeature, bool isPrecursor)
        {
            Feature bestFeature = null;
            int bestIntersectionArea = 0;
            Rectangle precursorBoundary = precursorFeature.GetBoundary();

            // TODO: this is not optimal
            FeatureSet features = GetFeatures(mz, isPrecursor);
            foreach (Feature feature in features)
            {
                if (precursorBoundary.Contains(feature.GetHighestPoint()))
                {
                    Rectangle boundary = feature.GetBoundary();
                    Rectangle intersection = Rectangle.Intersect(precursorBoundary, boundary);
                    int intersectionArea = intersection.Width * intersection.Height;
                    double portionIntersection = (double)intersectionArea / (boundary.Width * boundary.Height);
                    if (portionIntersection < 0.9f)
                        continue;
                    if (intersectionArea > bestIntersectionArea)
                    {
                        bestFeature = feature;
                        bestIntersectionArea = intersectionArea;
                    }
                }
            }

            return bestFeature;
        }
    }
}
