using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{

    public class Ms1FeatureResult
    {
        public static List<Ms1Feature> LoadProMexResult(string rawFilePath, string featureFilePath, int dataid = 0)
        {
            var featureList = new List<Ms1Feature>();
            var tsvReader = new TsvFileParser(featureFilePath);
            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);
            var featureIds = tsvReader.GetData("FeatureID");
            
            var minScans = tsvReader.GetData("QMinScanNum");
            var maxScans = tsvReader.GetData("QMaxScanNum");
            var abundance = tsvReader.GetData("QAbundance");
            //var minScans = tsvReader.GetData("MinScan");
            //var maxScans = tsvReader.GetData("MaxScan");
            //var abundance = tsvReader.GetData("Abundance");


            var minCharges = tsvReader.GetData("MinCharge");
            var maxCharges = tsvReader.GetData("MaxCharge");
            var monoMass = tsvReader.GetData("MonoMass");
            
            var repCharges = tsvReader.GetData("RepCharge");
            var repScans = tsvReader.GetData("RepScan");
            var repMzs = tsvReader.GetData("RepMz");

            for (var i = 0; i < tsvReader.NumData; i++)
            {
                var abu = double.Parse(abundance[i]);
                var mass = double.Parse(monoMass[i]);
                var minChg = int.Parse(minCharges[i]);
                var maxChg = int.Parse(maxCharges[i]);
                var minScan = int.Parse(minScans[i]);
                var maxScan = int.Parse(maxScans[i]);
                var fid = int.Parse(featureIds[i]);
                var repChg = int.Parse(repCharges[i]);
                var repMz = double.Parse(repMzs[i]);
                var repScanNum = int.Parse(repScans[i]);

                var feature = new Ms1Feature(run, minChg, maxChg, minScan, maxScan, abu, mass, repChg, repMz, repScanNum)
                {
                    FeatureId = fid, 
                    DataSetId = dataid,
                };
                
                featureList.Add(feature);
            }

            //featureList.Sort(new Ms1FeatureMassComparer());
            return featureList;            
        }
    }
    
    /*
    public class Ms1FeatureResult : Ms1Feature, IComparable<Ms1FeatureResult>
    {
        
        
        
        public Ms1FeatureResult(LcMsRun run, int dataid, int fid, double mass, int minChg, int maxChg, int minScan,
            int maxScan, double abu)
            : base(run, minChg, maxChg, minScan, maxScan, abu, mass, 
        {
            
//e(LcMsRun run, int minCharge, int maxCharge, int minScan, int maxScan, double abundance,
  //                       double repMass, int repCharge, double repMz, int repScanNum)            
            
            DataSetId = dataid;
            FeatureId = fid;
            RepresentativeMass = mass;
            MinCharge = minChg;
            MaxCharge = maxChg;
            MinScanNum = minScan;
            MaxScanNum = maxScan;
            Abundance = abu;
        }

        public int CompareTo(Ms1FeatureResult other)
        {
            return Mass.CompareTo(other.Mass);
        }
        
        
        public double GetMinNetDiff(Ms1Feature other)
        {
            if (CoEluted(other)) return 0d;

            var nd = new double[]
            {
                Math.Abs(MinNet - other.MinNet),
                Math.Abs(MinNet - other.MaxNet),
                Math.Abs(MaxNet - other.MinNet),
                Math.Abs(MaxNet - other.MaxNet)
            };
            return nd.Min();
        }


        public bool MassMatch(Ms1FeatureResult other, bool oneDaltonOk = false)
        {
            var tolerance = other.DataSetId == DataSetId ? new Tolerance(2.5) : new Tolerance(5);
            var massTol = tolerance.GetToleranceAsTh(other.Mass);
            var massDiff = Math.Abs(Mass - other.Mass);
            var sameMass = oneDaltonOk ? (massDiff < massTol || Math.Abs(massDiff - 1) < massTol) : (massDiff < massTol);
            return sameMass;
        }

        /*
        public readonly List<Ms1FeatureResult> MergedFeatures;
        public double[] GetMasses()
        {
            var ret = new double[1 + MergedFeatures.Count];
            ret[0] = Mass;
            for (var j = 0; j < MergedFeatures.Count; j++)
            {
                ret[j + 1] = MergedFeatures[j].Mass;
            }
            Array.Sort(ret);
            return ret;
        }

        public bool Equals(Ms1FeatureResult other)
        {
            if (MinCharge > other.MaxCharge || other.MinCharge > MaxCharge) return false;

            if (DataSetId == other.DataSetId)
            {
                var sameMass = MassMatch(other, false);
                if (!sameMass) return false;

                if (MinNet > other.MinNet - 0.005 && MinNet < other.MaxNet + 0.005) return true;
                if (other.MinNet > MinNet - 0.005 && other.MinNet < MaxNet + 0.005) return true;
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

                //var lenDiff = Math.Abs(NetLength - other.NetLength);
                //if (lenDiff > largeFeature.NetLength) return false;

                // when aligning biological replicate pair, do not consider elution length?
                if (largeFeature.MinNet - 0.33333 * largeFeature.NetLength < smallFeature.MinNet &&
                    smallFeature.MinNet < largeFeature.MaxNet + 0.33333 * largeFeature.NetLength) return true;

                if (largeFeature.MaxNet - 0.33333 * largeFeature.NetLength < smallFeature.MaxNet &&
                    smallFeature.MaxNet < largeFeature.MaxNet + 0.33333 * largeFeature.NetLength) return true;

            }

            return false;
        }

        public bool MassMatch(Ms1FeatureResult other, bool oneDaltonOk = false)
        {
            var sameMass = false;
            var tolerance = other.DataSetId == DataSetId ? new Tolerance(2.5) : new Tolerance(5);

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

        public void Merge(Ms1FeatureResult other)
        {
            if (DataSetId != other.DataSetId) return;

            MergedFeatures.Add(other);

            MaxScanNum = Math.Max(MaxScanNum, other.MaxScanNum);
            MinScanNum = Math.Min(MinScanNum, other.MinScanNum);

            MaxCharge = Math.Max(MaxCharge, other.MaxCharge);
            MinCharge = Math.Min(MinCharge, other.MinCharge);
            
            Abundance += other.Abundance;
        }         
        
        private bool _elutionTimeFit;
        private Tuple<double, double> _elutionTimeFitParam;
    }
    */
}
