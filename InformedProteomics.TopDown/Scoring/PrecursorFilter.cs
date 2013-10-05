using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Scoring
{
    public class PrecursorFilter
    {
        public PrecursorFilter(LcMsRun run, Tolerance toleranceForBaseXic, Tolerance toleranceFromBasicXic)
        {
            Run = run;
            ToleranceForBaseXic = toleranceForBaseXic;
            ToleranceFromBasicXic = toleranceFromBasicXic;
        }

        public LcMsRun Run { get; private set; }
        public Tolerance ToleranceForBaseXic { get; private set; }
        public Tolerance ToleranceFromBasicXic { get; private set; }

        public const double RelativeIsotopeIntensityThreshold = 0.5;

        public bool IsValid(Ion precursorIon, int scanNum)
        {
            if(Run.GetMsLevel(scanNum) != 2) return false;

            var isotopes = precursorIon.GetIsotopes(RelativeIsotopeIntensityThreshold).ToArray();

            var precursorScanNum = Run.GetPrecursorScanNum(scanNum);
            var nextMs1ScanNum = Run.GetNextScanNum(scanNum, 1);

            return
                IsValidForMs1Scan(precursorIon, precursorScanNum, isotopes)
                || IsValidForMs1Scan(precursorIon, nextMs1ScanNum, isotopes);
        }

        private bool IsValidForMs1Scan(Ion precursorIon, int scanNum, IEnumerable<Tuple<int,float>> isotopes)
        {
            if (Run.GetMsLevel(scanNum) != 1) return false;
            var spec = Run.GetSpectrum(scanNum);
            if (spec == null) return true;

            // check whether the most abundant isotope peak exists
            foreach (var isotope in isotopes)
            {
                var isotopeIndex = isotope.Item1;
                var isotopeMz = precursorIon.GetIsotopeMz(isotopeIndex);
                if (spec.FindPeak(isotopeMz, ToleranceForBaseXic) == null) return false;
            }

            return true;
        }
    }
}
