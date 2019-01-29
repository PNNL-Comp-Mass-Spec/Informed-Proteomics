using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.FeatureFinding.Data
{
    public static class LcMsRunExtensions
    {
        public static Ms1Spectrum GetMs1Spectrum(this LcMsRun run, int scanNum)
        {
            var spec = run.GetMs1Spectrum(scanNum, out var ms1ScanIndex);
            var ms1Spec = new Ms1Spectrum(scanNum, ms1ScanIndex, spec.Peaks);
            return ms1Spec;
        }
    }
}
