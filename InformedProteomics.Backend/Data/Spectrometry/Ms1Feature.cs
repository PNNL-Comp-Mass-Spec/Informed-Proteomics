using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Ms1Feature
    {
        public Ms1Feature(LcMsRun run)
        {
            Run = run;
            MinScanNum = 0;
            MaxScanNum = 0;
            MinCharge = 0;
            MaxCharge = 0;
        }

        public Ms1Feature(LcMsRun run, int minCharge, int maxCharge, int minScan, int maxScan, double abundance,
                         double repMass, int repCharge, double repMz, int repScanNum) : this(run)
        {
            MinCharge = minCharge;
            MaxCharge = maxCharge;
            MinScanNum = minScan;
            MaxScanNum = maxScan;
            Abundance = abundance;
            RepresentativeMass = repMass;
            RepresentativeCharge = repCharge;
            RepresentativeMz = repMz;
            RepresentativeScanNum = repScanNum;
        }

        public List<ObservedEnvelope> EnvelopeList;
        
        public int DataSetId { get; set; }
        public int FeatureId { get; set; }

        //public string Desc;
        
        public int MinCharge { get; protected set; }
        public int MaxCharge { get; protected set; }
        public int ChargeLength { get { return (MaxCharge == 0) ? 0 : MaxCharge - MinCharge + 1; } }

        public int MinScanNum { get; protected set; }
        public int MaxScanNum { get; protected set; }
        public int ScanLength { get { return (MaxScanNum == 0) ? 0 : MaxScanNum - MinScanNum + 1; } }

        public double Abundance { get; set; }
        
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
        
        public bool CoEluted(Ms1Feature other, double tolNet = 0d)
        {
            if (MinNet - tolNet < other.MinNet || other.MinNet < MaxNet + tolNet) return true;
            if (MinNet - tolNet < other.MaxNet || other.MaxNet < MaxNet + tolNet) return true;
            if (other.MinNet - tolNet < MinNet || MinNet < other.MaxNet + tolNet) return true;
            if (other.MinNet - tolNet < MaxNet || MaxNet < other.MaxNet + tolNet) return true;
            return false;
        }
        
        protected readonly LcMsRun Run;
    }

}
