using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public class LazyLcMsRun
    {
        public LazyLcMsRun(string rawFileName, double snr1, double snr2)
        {
            _rawFileName = rawFileName;
            _snr1 = snr1;
            _snr2 = snr2;
            _lcms = null;
        }
        public Spectrum GetSpectrum(int scanNum)
        {
            return Lcms.GetSpectrum(scanNum);
        }

        private LcMsRun Lcms
        {
            get { return _lcms ?? (_lcms = InMemoryLcMsRun.GetLcMsRun(_rawFileName, MassSpecDataType.XCaliburRun, _snr1, _snr2)); }
        }
        private readonly string _rawFileName;
        private readonly double _snr1;
        private readonly double _snr2;
        private LcMsRun _lcms;
    }
}
