using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using MultiDimensionalPeakFinding;
using UIMFLibrary;

namespace InformedProteomics.Backend.IMS
{
    public class ImsData
    {
        private readonly UimfUtil _uimfUtil;
        private readonly FeatureDetectionUtil _featureDetectionUtil;
        private readonly string _filePath;

        public ImsData(string filePath)
        {
            _filePath = filePath;
            _uimfUtil = new UimfUtil(filePath);
            _featureDetectionUtil = new FeatureDetectionUtil(filePath);
        }

        public FeatureSet GetFeatures(double mz, Tolerance tolerance, DataReader.FrameType frameType)
        {
            var intensityBlock = _uimfUtil.GetXic(mz, tolerance.GetValue(), frameType, tolerance.GetUnit());
            var features = new FeatureSet(intensityBlock);

            return features;
        }

        public Dictionary<int,FeatureSet> GetFeaturesMultiThreading(int minTargetBin, int maxTargetBin, Tolerance tolerance, DataReader.FrameType frameType)
        {
            var featureDictionary = new Dictionary<int, FeatureSet>();

            var targetMzList = Enumerable.Range(minTargetBin, maxTargetBin - minTargetBin + 1).Select(targetBin => _uimfUtil.GetMzFromBin(targetBin)).ToList();
            var featureStatisticsDictionary = _featureDetectionUtil.GetFeatureStatistics(targetMzList, tolerance.GetValue(), frameType, tolerance.GetUnit());
            foreach (var entry in featureStatisticsDictionary)
            {
                var targetBin = _uimfUtil.GetBinFromMz(entry.Key);
                var featureStats = entry.Value;
                featureDictionary.Add(targetBin, new FeatureSet(featureStats));
            }
            return featureDictionary;
        }

        public FeatureSet GetFeatures(int targetBin, DataReader.FrameType frameType)
        {
            List<IntensityPoint> intensityPointList = _uimfUtil.GetXic(targetBin, frameType);
            var features = new FeatureSet(intensityPointList);

            return features;
        }

        public int GetNumberOfBins()
        {
            return _uimfUtil.GetNumberOfBins();
        }

        public double GetMzFromBin(int bin)
        {
            return _uimfUtil.GetMzFromBin(bin);
        }

        public int GetBinFromMz(double mz)
        {
            return _uimfUtil.GetBinFromMz(mz);
        }

        public string GetFilePath()
        {
            return _filePath;
        }
    }
}
