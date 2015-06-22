using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.MassFeature;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Ms1Peak : Peak
    {
        public Ms1Peak(double mz, double intensity, int indexInSpec)
            : base(mz, intensity)
        {
            IndexInSpectrum = indexInSpec;
            Active = true;
            Quantified = false;
        }

        public int IndexInSpectrum { get; private set; }
        public bool Active { get; private set; }
        public bool Quantified { get; internal set; }

        public void InActivate()
        {
            Active = false;
        }


        public void TagMajorPeakOf(LcMsPeakCluster feature)
        {
            if (_majorIsotopeFeatures == null) _majorTaggedFeatures = new List<LcMsPeakCluster>();
            _majorTaggedFeatures.Add(feature);
        }
        public void TagMinorPeakOf(LcMsPeakCluster feature)
        {
            if (_minorTaggedFeatures == null) _minorTaggedFeatures = new List<LcMsPeakCluster>();
            _minorTaggedFeatures.Add(feature);
        }

        public IEnumerable<LcMsPeakCluster> GetAllTaggedFeatures()
        {
            if (_minorTaggedFeatures != null) foreach (var f in _minorTaggedFeatures) yield return f;
            if (_majorTaggedFeatures != null) foreach (var f in _majorTaggedFeatures) yield return f;
        }
        public IEnumerable<LcMsPeakCluster> GetMajorTaggedFeatures()
        {
            if (_majorTaggedFeatures != null) foreach (var f in _majorTaggedFeatures) yield return f;
        }

        
        private List<LcMsPeakCluster> _minorTaggedFeatures;
        private List<LcMsPeakCluster> _majorTaggedFeatures;


        public void TagMajorPeakOf(Ms1FeatureCluster feature)
        {
            if (_majorIsotopeFeatures == null) _majorIsotopeFeatures = new List<Ms1FeatureCluster>();
            _majorIsotopeFeatures.Add(feature);
        }
        public void TagMinorPeakOf(Ms1FeatureCluster feature)
        {
            if (_minorIsotopeFeatures == null) _minorIsotopeFeatures = new List<Ms1FeatureCluster>();
            _minorIsotopeFeatures.Add(feature);
        }

        public IEnumerable<Ms1FeatureCluster> GetTaggedAllFeatures()
        {
            if (_minorIsotopeFeatures != null) foreach (var f in _minorIsotopeFeatures) yield return f;
            if (_majorIsotopeFeatures != null) foreach (var f in _majorIsotopeFeatures) yield return f;
        }
        public IEnumerable<Ms1FeatureCluster> GetTaggedMajorFeatures()
        {
            if (_majorIsotopeFeatures != null) foreach (var f in _majorIsotopeFeatures) yield return f;
        }

        internal ushort Ms1SpecIndex { get; set; }
        private List<Ms1FeatureCluster> _minorIsotopeFeatures;
        private List<Ms1FeatureCluster> _majorIsotopeFeatures;
 
    }
   
}
