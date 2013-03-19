using System.Collections.Generic;
using DeconTools.Backend;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Utils;
using MultiDimensionalPeakFinding;
using UIMFLibrary;

namespace InformedProteomics.Backend.IMS
{
    public class ImsData
    {
        private readonly UimfUtil _uimfUtil;
        private readonly string _filePath;

        public ImsData(string filePath)
        {
            _filePath = filePath;
            _uimfUtil = new UimfUtil(filePath);
        }

        public FeatureSet GetFeatures(double mz, Tolerance tolerance, DataReader.FrameType frameType)
        {
            List<IntensityPoint> intensityBlock = _uimfUtil.GetXic(mz, tolerance.GetValue(), frameType, tolerance.GetUnit());
            var xic = new FeatureSet(intensityBlock);

            return xic;
        }

        public FeatureSet GetFeatures(int targetBin, DataReader.FrameType frameType)
        {
            List<IntensityPoint> intensityPointList = _uimfUtil.GetXic(targetBin, frameType);
            var xic = new FeatureSet(intensityPointList);

            return xic;
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
