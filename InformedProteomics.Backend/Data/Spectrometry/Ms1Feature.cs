namespace InformedProteomics.Backend.Data.Spectrometry
{
    public abstract class Ms1Feature
    {
        public abstract int MinCharge { get; }
        public abstract int MaxCharge { get; }
        public abstract int MinScanNum { get; }
        public abstract int MaxScanNum { get; }
        
        //public double Score { get; protected set; }
        public int RepresentativeScanNum { get; protected set; }
        public double RepresentativeMass { get; protected set; }
        public double RepresentativeMz { get; protected set; }
        public int RepresentativeCharge { get; protected set; }
    }
}
