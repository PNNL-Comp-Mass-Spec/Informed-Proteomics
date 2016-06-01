using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.FeatureFinding.MassFeature
{
    public class Ms1Peak : Peak
    {
        public Ms1Peak(double mz, double intensity, int indexInSpec)
            : base(mz, intensity)
        {
            IndexInSpectrum = indexInSpec;
            Active = true;
            _countMinorTaggedFeatures = 0;
            _countMajorTaggedFeatures = 0;
        }

        public void InActivate()
        {
            Active = false;
        }

        public void Activate()
        {
            Active = true;
        }

        public void TagMajorPeakOf(LcMsPeakCluster feature)
        {
            //if (_majorTaggedFeatures == null) _majorTaggedFeatures = new List<LcMsPeakCluster>();
            //_majorTaggedFeatures.Add(feature);
            if (_majorTaggedFeatures == null)
            {
                _majorTaggedFeatures = new LcMsPeakCluster[2];
                _majorTaggedFeatures[_countMajorTaggedFeatures++] = feature;
            }
            else
            {
                if (_countMajorTaggedFeatures >= _majorTaggedFeatures.Length)
                {
                    Array.Resize(ref _majorTaggedFeatures, _majorTaggedFeatures.Length * 2);
                }
                _majorTaggedFeatures[_countMajorTaggedFeatures++] = feature;
            }

        }
        public void TagMinorPeakOf(LcMsPeakCluster feature)
        {
            //if (_minorTaggedFeatures == null) _minorTaggedFeatures = new List<LcMsPeakCluster>();
            //_minorTaggedFeatures.Add(feature);
            if (_minorTaggedFeatures == null)
            {
                _minorTaggedFeatures = new LcMsPeakCluster[2];
                _minorTaggedFeatures[_countMinorTaggedFeatures++] = feature;
            }
            else
            {
                if (_countMinorTaggedFeatures >= _minorTaggedFeatures.Length)
                {
                    Array.Resize(ref _minorTaggedFeatures, _minorTaggedFeatures.Length * 2);
                }
                _minorTaggedFeatures[_countMinorTaggedFeatures++] = feature;
            }
        }

        public IEnumerable<LcMsPeakCluster> GetAllTaggedFeatures()
        {
            //if (_minorTaggedFeatures != null) foreach (var f in _minorTaggedFeatures) yield return f;
            //if (_majorTaggedFeatures != null) foreach (var f in _majorTaggedFeatures) yield return f;
            for (var i = 0; i < _countMajorTaggedFeatures; i++) yield return _majorTaggedFeatures[i];
            for (var i = 0; i < _countMinorTaggedFeatures; i++) yield return _minorTaggedFeatures[i];
        }
        public IEnumerable<LcMsPeakCluster> GetMajorTaggedFeatures()
        {
            //if (_majorTaggedFeatures != null) foreach (var f in _majorTaggedFeatures) yield return f;
            for (var i = 0; i < _countMajorTaggedFeatures; i++) yield return _majorTaggedFeatures[i];
        }

        public int IndexInSpectrum { get; private set; }
        public bool Active { get; private set; }
        public int Ms1SpecIndex { get; set; }
        //private List<LcMsPeakCluster> _minorTaggedFeatures;
        //private List<LcMsPeakCluster> _majorTaggedFeatures;

        private LcMsPeakCluster[] _minorTaggedFeatures;
        private int _countMinorTaggedFeatures;

        private LcMsPeakCluster[] _majorTaggedFeatures;
        private int _countMajorTaggedFeatures;


    }

}
