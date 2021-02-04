using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Scoring
{
    public class Ms1ContainsIonFilter
    {
        public Ms1ContainsIonFilter(LcMsRun run, Tolerance mzTolerance)
        {
            Run = run;
            MzTolerance = mzTolerance;
        }

        public LcMsRun Run { get; }

        public Tolerance MzTolerance { get; }

        public const double RelativeIsotopeIntensityThreshold = 0.8;    // 0.5

        public bool IsValid(Ion precursorIon, int scanNum)
        {
            if (Run.GetMsLevel(scanNum) != 2)
            {
                return false;
            }

            var precursorScanNum = Run.GetPrecursorScanNum(scanNum);
            var nextMs1ScanNum = Run.GetNextScanNum(scanNum);

            var isValid =
                IsValidForMs1Scan(precursorIon, precursorScanNum)
                || IsValidForMs1Scan(precursorIon, nextMs1ScanNum);

            return isValid;
        }

        private bool IsValidForMs1Scan(Ion precursorIon, int scanNum)
        {
            if (scanNum < Run.MinLcScan || scanNum > Run.MaxLcScan)
            {
                return false;
            }

            if (Run.GetMsLevel(scanNum) != 1)
            {
                return true;
            }

            var spec = Run.GetSpectrum(scanNum);
            return spec?.ContainsIon(precursorIon, MzTolerance, RelativeIsotopeIntensityThreshold) == true;
        }
    }
}
