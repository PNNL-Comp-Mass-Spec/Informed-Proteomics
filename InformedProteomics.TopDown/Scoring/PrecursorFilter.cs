using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Scoring
{
    public class PrecursorFilter
    {
        public PrecursorFilter(LcMsRun run, Tolerance mzTolerance)
        {
            Run = run;
            MzTolerance = mzTolerance;
        }

        public LcMsRun Run { get; private set; }
        public Tolerance MzTolerance { get; private set; }

        public const double RelativeIsotopeIntensityThreshold = 0.8;    // 0.5

        public bool IsValid(Ion precursorIon, int scanNum)
        {
            if (Run.GetMsLevel(scanNum) != 2) return false;
            var precursorScanNum = Run.GetPrecursorScanNum(scanNum);
            var nextMs1ScanNum = Run.GetNextScanNum(scanNum);

            var isValid =
                IsValidForMs1Scan(precursorIon, precursorScanNum)
                || IsValidForMs1Scan(precursorIon, nextMs1ScanNum);

            //Console.WriteLine("{0}\t{1}\t{2}\t{3}", precursorIon.GetMonoIsotopicMz(), precursorIon.precursorCharge, scanNum, isValid);
            return isValid;
        }

        private bool IsValidForMs1Scan(Ion precursorIon, int scanNum)
        {
            if (scanNum < Run.MinLcScan || scanNum > Run.MaxLcScan) return false;
            if (Run.GetMsLevel(scanNum) != 1) return true;
            var spec = Run.GetSpectrum(scanNum);
            return spec != null && spec.ContainsIon(precursorIon, MzTolerance, RelativeIsotopeIntensityThreshold);
        }

        //public bool IsValid2(Ion precursorIon, int scanNum)
        //{
        //    if (Run.GetMsLevel(scanNum) != 2) return false;

        //    var isotopes = precursorIon.GetIsotopes(RelativeIsotopeIntensityThreshold).ToArray();

        //    var precursorScanNum = Run.GetPrecursorScanNum(scanNum);
        //    var nextMs1ScanNum = Run.GetNextScanNum(scanNum, 1);

        //    var isValid =
        //        IsValidForMs1Scan2(precursorIon, precursorScanNum, isotopes)
        //        || IsValidForMs1Scan2(precursorIon, nextMs1ScanNum, isotopes);

        //    //Console.WriteLine("{0}\t{1}\t{2}\t{3}", precursorIon.GetMonoIsotopicMz(), precursorIon.precursorCharge, scanNum, isValid);
        //    return isValid;
        //}

        //private bool IsValidForMs1Scan2(Ion precursorIon, int scanNum, IEnumerable<Tuple<int, float>> isotopes)
        //{
        //    if (Run.GetMsLevel(scanNum) != 1) return false;
        //    var spec = Run.GetSpectrum(scanNum);
        //    if (spec == null) return true;

        //    // check whether the most abundant isotope peak exists
        //    foreach (var isotope in isotopes)
        //    {
        //        var isotopeIndex = isotope.Item1;
        //        var isotopeMz = precursorIon.GetIsotopeMz(isotopeIndex);
        //        if (spec.FindPeak(isotopeMz, MzTolerance) == null) return false;
        //    }

        //    return true;
        //}
    }
}
