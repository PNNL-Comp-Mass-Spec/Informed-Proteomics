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
        }

        public Ms1Peak() { } // For use in generics

        // For setting values in generics
        internal void SetMzIntensityIndices(double mz, double intensity, int indexInSpec, ushort ms1SpecIndex)
        {
            SetMzAndIntensity(mz, intensity);
            IndexInSpectrum = indexInSpec;
            Ms1SpecIndex = ms1SpecIndex;
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
            if (_majorTaggedFeatures == null) _majorTaggedFeatures = new List<LcMsPeakCluster>();
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

        public int IndexInSpectrum { get; private set; }
        public bool Active { get; private set; }
        public ushort Ms1SpecIndex { get; set; }
        
        private List<LcMsPeakCluster> _minorTaggedFeatures;
        private List<LcMsPeakCluster> _majorTaggedFeatures;
    }
   
}
