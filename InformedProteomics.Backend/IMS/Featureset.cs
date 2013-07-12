using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using MultiDimensionalPeakFinding;
using MultiDimensionalPeakFinding.PeakDetection;
using UIMFLibrary;

namespace InformedProteomics.Backend.IMS
{
    [Serializable]
    public class FeatureSet : IEnumerable<Feature>
    {
        private List<Feature> _featureList;

        public FeatureSet(List<IntensityPoint> intensityPointList)
        {
            var pointList = WaterShedMapUtil.BuildWatershedMap(intensityPointList);
            Smoother.Smooth(ref pointList);
            FindFeatures(pointList);
        }

        public FeatureSet(IEnumerable<FeatureBlobStatistics> featureStatistics)
        {
            foreach (var featureStat in featureStatistics)
            {
                _featureList.Add(new Feature(featureStat));
            }
        }

        public FeatureSet(double[,] intensityBlock)
        {
            Smoother.Smooth(ref intensityBlock); 
            FindFeatures(intensityBlock);
        }

        public IEnumerable<Feature> GetFeatures()
        {
            return _featureList;
        }

        public IEnumerable<FeatureBlobStatistics> GetFeatures(int frameIndexMin, int frameIndexMax, int scanMin, int scanMax)
        {
            return null;
        }

        #region Private Methods

        private void FindFeatures(IEnumerable<Point> pointList)
        {
            var featureBlobList = FeatureDetection.DoWatershedAlgorithm(pointList);
            var featureBlobStatList = featureBlobList.Select(featureBlob => featureBlob.Statistics).ToList();

            _featureList = featureBlobStatList.Select(featureBlobStatistics => new Feature(featureBlobStatistics)).ToList();
        }

        private void FindFeatures(double[,] intensityBlock)
        {
            var pointList = WaterShedMapUtil.BuildWatershedMap(intensityBlock);
            var featureBlobList = FeatureDetection.DoWatershedAlgorithm(pointList);
            var featureBlobStatList = featureBlobList.Select(featureBlob => featureBlob.Statistics).ToList();

            _featureList = featureBlobStatList.Select(featureBlobStatistics => new Feature(featureBlobStatistics)).ToList();
        }

        #endregion

        private static readonly SavitzkyGolaySmoother Smoother = new SavitzkyGolaySmoother(9, 2);

        public IEnumerator<Feature> GetEnumerator()
        {
            return _featureList.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }
}
