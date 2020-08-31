using System;
using InformedProteomics.FeatureFinding.SpectrumMatching;

namespace InformedProteomics.FeatureFinding.Data
{
    public class LcMsFeature
    {
        // Ignore Spelling: Da

        public LcMsFeature(double repMass, int repCharge, double repMz, int repScanNum, double abundance)
            : this(repMass, repCharge, repMz, repScanNum, abundance,
            repCharge, repCharge, repScanNum, repScanNum, 0, 0)
        {
        }

        public LcMsFeature(double repMass, int repCharge, double repMz, int repScanNum, double abundance,
            int minCharge, int maxCharge, int minScan, int maxScan,
            double minElution, double maxElution, double minNet = 0, double maxNet = 0)
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

            MinElutionTime = minElution;
            MaxElutionTime = maxElution;

            MaxNet = maxNet;
            MinNet = minNet;

            ProteinSpectrumMatches = new ProteinSpectrumMatchSet(0);
        }

        public int DataSetId
        {
            get => ProteinSpectrumMatches.DataId;
            set => ProteinSpectrumMatches.SetDataId(value);
        }

        public int FeatureId { get; set; }

        public int MinCharge { get; protected set; }
        public int MaxCharge { get; protected set; }
        public int ChargeLength => MaxCharge - MinCharge + 1;

        public int MinScanNum { get; protected set; }
        public int MaxScanNum { get; protected set; }
        public int ScanLength => MaxScanNum - MinScanNum + 1;

        public double Abundance { get; protected set; }

        /// <summary>
        /// LikelihoodRatio
        /// </summary>
        public double Score { get; set; }

        public double Mass => RepresentativeMass;
        public double RepresentativeMass { get; protected set; }
        public int RepresentativeCharge { get; protected set; }
        public int RepresentativeScanNum { get; protected set; }
        public double RepresentativeMz { get; protected set; }

        public virtual double MinElutionTime { get; }
        public virtual double MaxElutionTime { get; }
        public double ElutionTime => 0.5 * (MinElutionTime + MaxElutionTime);
        public double ElutionLength => MaxElutionTime - MinElutionTime;

        /*
        public double MinElutionTime { get { return (Run == null) ? _minElution : Run.GetElutionTime(MinScanNum); } }
        public double MaxElutionTime { get { return (Run == null) ? _maxElution : Run.GetElutionTime(MaxScanNum); } }
        public double ElutionTime { get { return 0.5*(MinElutionTime + MaxElutionTime); } }
        public double ElutionLength { get { return MaxElutionTime - MinElutionTime; } }

        public double MaxNet { get { return (Run == null) ? 0 : MaxElutionTime / Run.GetElutionTime(Run.MaxLcScan); } }
        public double MinNet { get { return (Run == null) ? 0 : MinElutionTime / Run.GetElutionTime(Run.MaxLcScan); } }
        public double NetLength { get { return (Run == null) ? 0 : MaxNet - MinNet; } }
        public double Net { get { return (Run == null) ? 0 : 0.5 * (MinNet + MaxNet); } }
        */
        public virtual double MaxNet { get; private set; }
        public virtual double MinNet { get; private set; }
        public double NetLength => MaxNet - MinNet;
        public double Net => 0.5 * (MinNet + MaxNet);

        public void SetNet(int minNet, int maxNet)
        {
            MaxNet = maxNet;
            MinNet = minNet;
        }

        public double CoElutionLength(LcMsFeature other)
        {
            if (other.MaxScanNum >= MinScanNum && other.MaxScanNum <= MaxScanNum)
            {
                return (other.MaxElutionTime - Math.Max(MinElutionTime, other.MinElutionTime));
            }

            if (MaxScanNum >= other.MinScanNum && MaxScanNum <= other.MaxScanNum)
            {
                return (MaxElutionTime - Math.Max(other.MinElutionTime, MinElutionTime));
            }

            return 0d;
        }

        public double CoElutionNetLength(LcMsFeature other)
        {
            if (other.MaxNet >= MinNet && other.MaxNet<= MaxNet)
            {
                return (other.MaxNet - Math.Max(MinNet, other.MinNet));
            }

            if (MaxNet >= other.MinNet && MaxNet <= other.MaxNet)
            {
                return (MaxNet - Math.Max(other.MinNet, MinNet));
            }

            return 0d;
        }

        public bool CoElutedByTime(LcMsFeature other, double tolTime = 0d)
        {
            if (MinElutionTime - tolTime < other.MinElutionTime && other.MinElutionTime < MaxElutionTime + tolTime) return true;
            if (MinElutionTime - tolTime < other.MaxElutionTime && other.MaxElutionTime < MaxElutionTime + tolTime) return true;
            if (other.MinElutionTime - tolTime < MinElutionTime && MinElutionTime < other.MaxElutionTime + tolTime) return true;
            if (other.MinElutionTime - tolTime < MaxElutionTime && MaxElutionTime < other.MaxElutionTime + tolTime) return true;
            return false;
        }

        public bool CoElutedByNet(LcMsFeature other, double tolNet = 0d)
        {
            if (MinNet - tolNet < other.MinNet && other.MinNet < MaxNet + tolNet) return true;
            if (MinNet - tolNet < other.MaxNet && other.MaxNet < MaxNet + tolNet) return true;
            if (other.MinNet - tolNet < MinNet && MinNet < other.MaxNet + tolNet) return true;
            if (other.MinNet - tolNet < MaxNet && MaxNet < other.MaxNet + tolNet) return true;
            return false;
        }

        public bool CoElutedByScanNum(LcMsFeature other, int tolScan = 0)
        {
            if (tolScan == 0) tolScan++;

            if (MinScanNum - tolScan < other.MinScanNum && other.MinScanNum < MaxScanNum + tolScan) return true;
            if (MinScanNum - tolScan < other.MaxScanNum && other.MaxScanNum < MaxScanNum + tolScan) return true;
            if (other.MinScanNum - tolScan < MinScanNum && MinScanNum < other.MaxScanNum + tolScan) return true;
            if (other.MinScanNum - tolScan < MaxScanNum && MaxScanNum < other.MaxScanNum + tolScan) return true;
            return false;
        }

        public readonly ProteinSpectrumMatchSet ProteinSpectrumMatches;

        /// <inheritdoc />
        public override string ToString()
        {
            return RepresentativeMass.ToString("0.0") + " Da, " +
                RepresentativeCharge + "+, " +
                RepresentativeMz.ToString("0.0") + " m/z, scan " +
                RepresentativeScanNum;
        }
    }
}
