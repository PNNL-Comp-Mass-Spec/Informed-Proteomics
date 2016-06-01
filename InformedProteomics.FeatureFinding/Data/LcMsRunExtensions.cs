using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.FeatureFinding.Data
{
    public static class LcMsRunExtensions
    {
        public static Ms1Spectrum GetMs1Spectrum(this LcMsRun run, int scanNum)
        {
            int ms1ScanIndex;
            var spec = run.GetMs1Spectrum(scanNum, out ms1ScanIndex);
            var ms1Spec = new Ms1Spectrum(scanNum, ms1ScanIndex, spec.Peaks);
            return ms1Spec;
        }
    }
}
