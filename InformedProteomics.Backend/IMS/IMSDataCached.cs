using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using UIMFLibrary;

namespace InformedProteomics.Backend.IMS
{
    public class ImsDataCached : ImsData
    {
        private readonly Dictionary<int, FeatureSet> _precursorFeatureSetMap;
        private readonly bool[] _isPrecursorCached;

        private readonly Dictionary<int, FeatureSet> _fragmentFeatureSetMap;
        private readonly bool[] _isFragmentCached;

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

            _precursorFeatureSetMap = new Dictionary<int, FeatureSet>();
            _fragmentFeatureSetMap = new Dictionary<int, FeatureSet>();

            var numBins = GetNumberOfBins();
            _isFragmentCached = new bool[numBins];
            _isPrecursorCached = new bool[numBins];
        }

        public double MinPrecursorMz { get; private set; }
        public double MaxPrecursorMz { get; private set; }
        public double MinFragmentMz { get; private set; }
        public double MaxFragmentMz { get; private set; }
        public Tolerance PrecursorTolerance { get; private set; }
        public Tolerance FragmentTolerance { get; private set; }

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

        public int CreatePrecursorFeatures()
        {
            return CreateFeatures(true);
        }

        public int CreateFragmentFeatures()
        {
            return CreateFeatures(false);
        }

        private int CreateFeatures(bool isPrecursor)
        {
            var frameType = isPrecursor ? DataReader.FrameType.MS1 : DataReader.FrameType.MS2;
            var featureSetMap = isPrecursor ? _precursorFeatureSetMap : _fragmentFeatureSetMap;
            var isCached = isPrecursor ? _isPrecursorCached : _isFragmentCached;
            var tolerance = isPrecursor ? PrecursorTolerance : FragmentTolerance;
            var minTargetBin = isPrecursor ? GetBinFromMz(MinPrecursorMz) : GetBinFromMz(MinFragmentMz);
            var maxTargetBin = isPrecursor ? GetBinFromMz(MaxPrecursorMz) : GetBinFromMz(MaxFragmentMz);
            if (maxTargetBin >= GetNumberOfBins())
                maxTargetBin = GetNumberOfBins() - 1;

            int totalNumFeatures = 0;

            for (int targetBin = minTargetBin; targetBin <= maxTargetBin; targetBin++)
            {
                double mz = GetMzFromBin(targetBin);
                //if(targetBin == 81635)
                //    Console.WriteLine("***MzBin: " + targetBin + "\tMz: " + mz);
                FeatureSet featureSet = GetFeatures(mz, tolerance, frameType);
                if (featureSet.GetFeatures().Any())
                    featureSetMap.Add(targetBin, featureSet);
                totalNumFeatures += featureSet.GetFeatures().Count();
                isCached[targetBin] = true;
            }

            return totalNumFeatures;
        }

        private FeatureSet GetFeatures(double mz, bool isPrecursor)
        {
            int mzBin = GetBinFromMz(mz);
            //Console.WriteLine("*****MzBin: " + mzBin + "\tMz: " + mz);
            var featureSetMap = isPrecursor ? _precursorFeatureSetMap : _fragmentFeatureSetMap;
            var isCached = isPrecursor ? _isPrecursorCached : _isFragmentCached;

            if (isCached[mzBin])
                return featureSetMap[mzBin];
            double recoveredMz = GetMzFromBin(mzBin);
            var featureSet = isPrecursor ?
                                 GetFeatures(recoveredMz, PrecursorTolerance, DataReader.FrameType.MS1)
                                 : GetFeatures(recoveredMz, FragmentTolerance, DataReader.FrameType.MS2);
            featureSetMap.Add(mzBin, featureSet);
            isCached[mzBin] = true;
            return featureSet;
        }

        private Feature GetFeature(double mz, Feature precursorFeature, bool isPrecursor)
        {
            Feature bestFeature = null;
            const float portionIntersectioinThreshold = 0.8f; // added by Kyowon - testing
            int bestIntersectionArea = 0;
            Rectangle precursorBoundary = precursorFeature.GetBoundary();

            // TODO: this may not be optimal
            FeatureSet features = GetFeatures(mz, isPrecursor);
            foreach (Feature feature in features)
            {
                if (precursorBoundary.Contains(feature.GetHighestPoint()))
                {
                    Rectangle boundary = feature.GetBoundary();
                    Rectangle intersection = Rectangle.Intersect(precursorBoundary, boundary);
                    int intersectionArea = intersection.Width * intersection.Height;
                    double portionIntersection = (double)intersectionArea / (boundary.Width * boundary.Height);
                    if (portionIntersection < portionIntersectioinThreshold)
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
