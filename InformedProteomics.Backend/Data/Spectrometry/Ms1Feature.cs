using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public abstract class Ms1Feature
    {
        protected Ms1Feature(LcMsRun run)
        {
            Run = run;
            MinScanNum = 0;
            MaxScanNum = 0;
            MinCharge = 0;
            MaxCharge = 0;
        }
        
        public string DataSetId { get; protected set; }
        public int FeatureId { get; protected set; }
        
        public int MinCharge { get; protected set; }
        public int MaxCharge { get; protected set; }
        public int ChargeLength { get { return (MaxCharge == 0) ? 0 : MaxCharge - MinCharge + 1; } }

        public int MinScanNum { get; protected set; }
        public int MaxScanNum { get; protected set; }
        public int ScanLength { get { return (MaxScanNum == 0) ? 0 : MaxScanNum - MinScanNum + 1; } }

        public double Abundance { get; protected set; }
        
        public double Mass { get { return RepresentativeMass;  } }
        public double RepresentativeMass { get; protected set; }
        public int RepresentativeCharge { get; protected set; }
        public int RepresentativeScanNum { get; protected set; }
        public double RepresentativeMz { get; protected set; }

        public double MinElutionTime { get { return Run.GetElutionTime(MinScanNum); }  }
        public double MaxElutionTime { get { return Run.GetElutionTime(MaxScanNum); }  }
        public double ElutionTime { get { return 0.5 * (MinElutionTime + MaxElutionTime); } }
        public double ElutionLength { get { return MaxElutionTime - MinElutionTime; } }
        
        public double MaxNet { get { return MaxElutionTime / Run.GetElutionTime(Run.MaxLcScan); } }
        public double MinNet { get { return MinElutionTime / Run.GetElutionTime(Run.MaxLcScan); } }
        public double NetLength { get { return MaxNet - MinNet; } }
        public double Net { get { return 0.5 * (MinNet + MaxNet); } }

        protected readonly LcMsRun Run;
    }

    /*
    public class Ms1Feature : IMs1Feature, IComparable<Ms1Feature>, IEquatable<Ms1Feature>
    {
        public Ms1Feature(LcMsRun run, int dataId, int featureId, double mass, int minCharge, int maxCharge, int minScan, int maxScan,  double abundance)
        {
            DataSetId = dataId;
            FeatureId = featureId;
            Mass = mass;
            MinCharge = minCharge;
            MaxCharge = maxCharge;
            MinScanNum = minScan;
            MaxScanNum = maxScan;
            Abundance = abundance;
            Active = true;
            MergedFeatures = new List<Ms1Feature>();
            _elutionTimeFit = false;
        }

        public int DataSetId { get; private set; }
        public int FeatureId { get; private set; }
        public int MinCharge { get; private set; }
        public int MaxCharge { get; private set; }
        public int MinScanNum { get; private set; }
        public int MaxScanNum { get; private set; }
        public double Abundance { get; private set; }
        public double Mass { get; private set; }

        public bool Active;
        public readonly List<Ms1Feature> MergedFeatures;

        public double[] GetMasses()
        {
            var ret = new double[1 + MergedFeatures.Count];
            ret[0] = Mass;
            for (var j = 0; j < MergedFeatures.Count; j++)
            {
                ret[j+1] = MergedFeatures[j].Mass;
            }
            Array.Sort(ret);
            return ret;
        }

        public bool Equals(Ms1Feature other)
        {
            if (MinCharge > other.MaxCharge || other.MinCharge > MaxCharge) return false;

            if (DataSetId == other.DataSetId)
            {
                var sameMass = MassMatch(other, true);
                if (!sameMass) return false;
                
                if (MinNet > other.MinNet - 0.001 && MinNet < other.MaxNet + 0.001) return true;
                if (other.MinNet > MinNet - 0.001 && other.MinNet < MaxNet + 0.001) return true;
                // for histone
                //if (MinNet > other.MinNet - 0.00 && MinNet < other.MaxNet + 0.00) return true;
                //if (other.MinNet > MinNet - 0.00 && other.MinNet < MaxNet + 0.00) return true;                
            }
            else
            {
                var sameMass = MassMatch(other);
                if (!sameMass) return false;
                
                Ms1Feature largeFeature = null;
                Ms1Feature smallFeature = null;
                if (NetLength > other.NetLength)
                {
                    largeFeature = this;
                    smallFeature = other;
                }
                else
                {
                    largeFeature = other;
                    smallFeature = this;
                }
                
                
                var lenDiff = Math.Abs(NetLength - other.NetLength);
                if (lenDiff > largeFeature.NetLength*0.3333) return false;

                // when aligning biological replicate pair, do not consider elution length?
                if (largeFeature.MinNet - 0.33333*largeFeature.NetLength < smallFeature.MinNet &&
                    smallFeature.MinNet < largeFeature.MinNet + 0.33333*largeFeature.NetLength) return true;
                
                if  (largeFeature.MaxNet - 0.33333 * largeFeature.NetLength < smallFeature.MaxNet &&
                    smallFeature.MaxNet < largeFeature.MaxNet - 0.33333 * largeFeature.NetLength) return true;
                
            }

            return false;
        }

        public bool MassMatch(Ms1Feature other, bool oneDaltonOk = false)
        {
            var sameMass = false;
            var tolerance = new Tolerance(5);
            var massTol = tolerance.GetToleranceAsTh(other.Mass);

            foreach (var m1 in GetMasses())
            {
                foreach (var m2 in other.GetMasses())
                {
                    var massDiff = Math.Abs(m1 - m2);
                    sameMass = oneDaltonOk ? (massDiff < massTol || Math.Abs(massDiff - 1) < massTol) : (massDiff < massTol);
                    if (sameMass) break;
                }
            }
            return sameMass;
        }

        public double NetLength
        {
            get { return MaxNet - MinNet; }
        }

        public int CompareTo(Ms1Feature other)
        {
            return Mass.CompareTo(other.Mass);
        }

        public double Net { get { return 0.5*(MinNet + MaxNet);  } }
        public double ElutionTime { get { return 0.5*(MinElutionTime + MaxElutionTime);  } }

        public double MinElutionTime
        {
            get
            {
                if (_elutionTimeFit)
                {
                    return Math.Max(_elutionTimeFitParam.Item1 + Run.GetElutionTime(MinScanNum)*_elutionTimeFitParam.Item2,
                        Run.GetElutionTime(Run.MinLcScan));
                }
                return Run.GetElutionTime(MinScanNum);
            }
        }

        public double MaxElutionTime
        {
            get
            {
                if (_elutionTimeFit)
                {
                    return Math.Min(_elutionTimeFitParam.Item1 + Run.GetElutionTime(MaxScanNum) * _elutionTimeFitParam.Item2,
                                            Run.GetElutionTime(Run.MaxLcScan));
                }
                return Run.GetElutionTime(MaxScanNum);
            }
        }

        public double MaxNet
        {
            get { return MaxElutionTime / Run.GetElutionTime(Run.MaxLcScan); }
        }

        public double MinNet
        {
            get { return MinElutionTime / Run.GetElutionTime(Run.MaxLcScan); }
        }
        
        public void Merge(Ms1Feature other)
        {
            if (DataSetId != other.DataSetId) return;

            MergedFeatures.Add(other);

            MaxScanNum = Math.Max(MaxScanNum, other.MaxScanNum);
            MinScanNum = Math.Min(MinScanNum, other.MinScanNum);

            MaxCharge = Math.Max(MaxCharge, other.MaxCharge);
            MinCharge = Math.Min(MinCharge, other.MinCharge);
            //_maxNet = Math.Max(MaxNet, other.MaxNet);
            //_minNet = Math.Min(MinNet, other.MinNet);

            Abundance += other.Abundance;
        }

        public void SetElutionTimeFitParam(Tuple<double, double> param)
        {
            //MinScanNum = minScan;
            //MaxScanNum = maxScan;
            _elutionTimeFitParam = param;
            _elutionTimeFit = true;
        }

        
        //private double _minNet;
        //private double _maxNet;
        public LcMsRun Run;

        private bool _elutionTimeFit;
        private Tuple<double, double> _elutionTimeFitParam;
    }
    */
}
