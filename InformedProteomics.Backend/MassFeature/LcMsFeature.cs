using System;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsFeature
    {
        public LcMsFeature(double repMass, int repCharge, double repMz, int repScanNum, double abundance,
                          int minCharge = -1, int maxCharge = -1, int minScan = -1, int maxScan = -1,
                          LcMsRun run = null)
        {
            Abundance = abundance;
            RepresentativeMass = repMass;
            RepresentativeCharge = repCharge;
            RepresentativeMz = repMz;
            RepresentativeScanNum = repScanNum;
            
            MinCharge = (minCharge > 0) ? minCharge : repCharge;
            MaxCharge = (maxCharge > 0) ? maxCharge : repCharge;
            MinScanNum = (minScan > 0) ? minScan : repScanNum;
            MaxScanNum = (maxScan > 0) ? maxScan : repScanNum;
            Run = run;
        }
        
        public int DataSetId { get; set; }
        public int FeatureId { get; set; }
        
        public int MinCharge { get; protected set; }
        public int MaxCharge { get; protected set; }
        public int ChargeLength { get { return MaxCharge - MinCharge + 1; } }

        public int MinScanNum { get; protected set; }
        public int MaxScanNum { get; protected set; }
        public int ScanLength { get { return MaxScanNum - MinScanNum + 1; } }

        public double Abundance { get; protected set; }
        
        public double Mass { get { return RepresentativeMass;  } }
        public double RepresentativeMass { get; protected set; }
        public int RepresentativeCharge { get; protected set; }
        public int RepresentativeScanNum { get; protected set; }
        public double RepresentativeMz { get; protected set; }

        public double MinElutionTime 
        { 
            get { return (Run == null) ? 0 : Run.GetElutionTime(MinScanNum); }  
        }
        public double MaxElutionTime { get { return (Run == null) ? 0 : Run.GetElutionTime(MaxScanNum); } }
        public double ElutionTime { get { return (Run == null) ? 0 : 0.5 * (MinElutionTime + MaxElutionTime); } }
        public double ElutionLength { get { return (Run == null) ? 0 : MaxElutionTime - MinElutionTime; } }

        public double MaxNet { get { return (Run == null) ? 0 : MaxElutionTime / Run.GetElutionTime(Run.MaxLcScan); } }
        public double MinNet { get { return (Run == null) ? 0 : MinElutionTime / Run.GetElutionTime(Run.MaxLcScan); } }
        public double NetLength { get { return (Run == null) ? 0 : MaxNet - MinNet; } }
        public double Net { get { return (Run == null) ? 0 : 0.5 * (MinNet + MaxNet); } }
        
        public bool CoElutedByNet(LcMsFeature other, double tolNet = 0d)
        {
            if (Run == null || other.Run == null)
            {
                throw new LcMsRunNullException();
            }
            
            if (MinNet - tolNet < other.MinNet || other.MinNet < MaxNet + tolNet) return true;
            if (MinNet - tolNet < other.MaxNet || other.MaxNet < MaxNet + tolNet) return true;
            if (other.MinNet - tolNet < MinNet || MinNet < other.MaxNet + tolNet) return true;
            if (other.MinNet - tolNet < MaxNet || MaxNet < other.MaxNet + tolNet) return true;
            return false;
        }
        public bool CoElutedByScanNum(LcMsFeature other, int tolScan = 0)
        {
            tolScan++;

            if (MinScanNum - tolScan < other.MinScanNum || other.MinScanNum < MaxScanNum + tolScan) return true;
            if (MinScanNum - tolScan < other.MaxScanNum || other.MaxScanNum < MaxScanNum + tolScan) return true;
            if (other.MinScanNum - tolScan < MinScanNum || MinScanNum < other.MaxScanNum + tolScan) return true;
            if (other.MinScanNum - tolScan < MaxScanNum || MaxScanNum < other.MaxScanNum + tolScan) return true;
            return false;
        }

        public LcMsRun Run { get; protected set; }

        [Serializable]
        public class LcMsRunNullException : Exception { }

    }
}
