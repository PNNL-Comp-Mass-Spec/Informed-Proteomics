using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.Backend.Quantification
{
    public class Ms2Feature
    {
        public int ScanNum { get; set; }
        public int Charge { get; set; }
        public double Mass { get; set; }

        public double Net { get; set; }

        public string DataSetId { get; set; }
        public int Id { get; set; }
        public string Sequence { get; set; }

        public bool LinkedToSameMs1Feature(Ms2Feature other)
        {
            return Ms1Feature.Ms2FeatureSet.Contains(other);
        }

        public QuantifiedMs1Feature Ms1Feature 
        {
            get { return _ms1Feature; }
            set
            {
                _ms1Feature = value;
                _ms1Feature.LinkMs2Feature(this);
            }
        }

        private QuantifiedMs1Feature _ms1Feature;
    }

    public class QuantifiedMs1Feature : Ms1FeatureCluster
    {
        public QuantifiedMs1Feature(LcMsRun run, byte minCharge, IsotopeList isoList, double repMass, int repCharge, double repMz, int repScanNum) : 
            base(run, minCharge, isoList, repMass, repCharge, repMz, repScanNum)
        {
            Ms2FeatureSet = new HashSet<Ms2Feature>();
        }

        public void LinkMs2Feature(Ms2Feature m2)
        {
            Ms2FeatureSet.Add(m2);
        }

        public readonly HashSet<Ms2Feature> Ms2FeatureSet;
    }
}
