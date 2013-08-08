using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using UIMFLibrary;

namespace InformedProteomics.Backend.IMS
{
    public class ImsDataCached : ImsData
    {
        private readonly ConcurrentDictionary<int, FeatureSet> _precursorFeatureSetMap;
        //private readonly bool[] _isPrecursorCached;
        private bool _isAllPrecursorFeaturesCached = false;

        private readonly ConcurrentDictionary<int, FeatureSet> _fragmentFeatureSetMap;
        //private readonly bool[] _isFragmentCached;
        private bool _isAllFragmentFeaturesCached = false;

        public ImsDataCached(string filePath) : this(filePath, 400.0, 1000.0, 10.0, 2000.0, 
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

            _precursorFeatureSetMap = new ConcurrentDictionary<int, FeatureSet>();
            _fragmentFeatureSetMap = new ConcurrentDictionary<int, FeatureSet>();
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
            _isAllPrecursorFeaturesCached = true;
            return CreateFeatures(true);
        }

        public int CreateFragmentFeatures()
        {
            _isAllFragmentFeaturesCached = true;
            return CreateFeatures(false);
        }

        public int CreatePrecursorFeaturesMultiThreads()
        {
            _isAllPrecursorFeaturesCached = true;
            return CreateFeaturesMultiThreading(true);
        }

        public int CreateFragmentFeaturesMultiThreads()
        {
            _isAllFragmentFeaturesCached = true;
            return CreateFeaturesMultiThreading(false);
        }

        private int CreateFeatures(bool isPrecursor)
        {
            var frameType = isPrecursor ? DataReader.FrameType.MS1 : DataReader.FrameType.MS2;
            var featureSetMap = isPrecursor ? _precursorFeatureSetMap : _fragmentFeatureSetMap;
            var tolerance = isPrecursor ? PrecursorTolerance : FragmentTolerance;
            var minTargetBin = isPrecursor ? GetBinFromMz(MinPrecursorMz) : GetBinFromMz(MinFragmentMz);
            var maxTargetBin = isPrecursor ? GetBinFromMz(MaxPrecursorMz) : GetBinFromMz(MaxFragmentMz);
            if (maxTargetBin >= GetNumberOfBins())
                maxTargetBin = GetNumberOfBins() - 1;

            var totalNumFeatures = 0;

            for (var targetBin = minTargetBin; targetBin <= maxTargetBin; targetBin++)
            {
                var mz = GetMzFromBin(targetBin);
                var featureSet = GetFeatures(mz, tolerance, frameType);
                featureSetMap.TryAdd(targetBin, featureSet);
                totalNumFeatures += featureSet.GetFeatures().Count();
            }

            return totalNumFeatures;
        }

        private void WriteFeatures(bool isPrecursor)
        {
            
        }

        private int CreateFeaturesMultiThreading(bool isPrecursor)
        {
            var frameType = isPrecursor ? DataReader.FrameType.MS1 : DataReader.FrameType.MS2;
            var featureSetMap = isPrecursor ? _precursorFeatureSetMap : _fragmentFeatureSetMap;
            var tolerance = isPrecursor ? PrecursorTolerance : FragmentTolerance;
            var minTargetBin = isPrecursor ? GetBinFromMz(MinPrecursorMz) : GetBinFromMz(MinFragmentMz);
            var maxTargetBin = isPrecursor ? GetBinFromMz(MaxPrecursorMz) : GetBinFromMz(MaxFragmentMz);
            if (maxTargetBin >= GetNumberOfBins())
                maxTargetBin = GetNumberOfBins() - 1;

            var featureSetDictionary = GetFeaturesMultiThreading(minTargetBin, maxTargetBin, tolerance, frameType);
            var totalNumFeatures = 0;
            foreach (var entry in featureSetDictionary)
            {
                var targetBin = entry.Key;
                var featureSet = entry.Value;
                featureSetMap.TryAdd(targetBin, featureSet);
                totalNumFeatures += featureSet.Count();
            }

            return totalNumFeatures;
        }

        private FeatureSet GetFeatures(double mz, bool isPrecursor)
        {
            var mzBin = GetBinFromMz(mz);
            var featureSetMap = isPrecursor ? _precursorFeatureSetMap : _fragmentFeatureSetMap;

            var isAllFeaturesCached = isPrecursor ? _isAllPrecursorFeaturesCached : _isAllFragmentFeaturesCached;

            FeatureSet curFeatureSet;
            if (featureSetMap.TryGetValue(mzBin, out curFeatureSet))
            {
                return curFeatureSet;
            }

            var recoveredMz = GetMzFromBin(mzBin);
            var featureSet = isPrecursor ?
                                 GetFeatures(recoveredMz, PrecursorTolerance, DataReader.FrameType.MS1)
                                 : GetFeatures(recoveredMz, FragmentTolerance, DataReader.FrameType.MS2);
            featureSetMap.TryAdd(mzBin, featureSet);
            return featureSet;
        }

        private Feature GetFeature(double mz, Feature precursorFeature, bool isPrecursor)
        {
            Feature bestFeature = null;
            const float portionIntersectioinThreshold = 0.1f; // added by Kyowon - testing
            int bestIntersectionArea = 0;
            var precursorBoundary = precursorFeature.GetBoundary();
            
            // TODO: this may not be optimal
            var features = GetFeatures(mz, isPrecursor);
            foreach (var feature in features)
            {
                if (precursorBoundary.Contains(feature.GetHighestPoint()))
                {
                    var boundary = feature.GetBoundary();
                    var intersection = Rectangle.Intersect(precursorBoundary, boundary);
                    var intersectionArea = intersection.Width * intersection.Height;
                    var portionIntersection = (double)intersectionArea / (precursorBoundary.Width * precursorBoundary.Height);
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
