using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Extension methods for ILcMsRun and ISpectrumAccessor
    /// </summary>
    public static class LcMsRunExtensions
    {
        /// <summary>
        /// Produce a summed spectrum using the data in the scans specified by <paramref name="scanNums"/>
        /// </summary>
        /// <param name="lcmsRun"></param>
        /// <param name="scanNums"></param>
        /// <param name="repScanNum">Representative scan number</param>
        /// <returns>Summed spectrum</returns>
        public static SummedSpectrum GetSummedSpectrum(this ISpectrumAccessor lcmsRun, IList<int> scanNums, int repScanNum = 0)
        {
            var mzComparer = new MzComparerWithBinning();

            var mzDic = new Dictionary<int, List<double>>();
            var intDic = new Dictionary<int, double>();

            foreach (var scanNum in scanNums)
            {
                var spec = lcmsRun.GetSpectrum(scanNum);
                if (spec == null)
                {
                    continue;
                }

                foreach (var peak in spec.Peaks)
                {
                    var mzBin = mzComparer.GetBinNumber(peak.Mz);

                    if (mzDic.TryGetValue(mzBin, out var mzList))
                    {
                        mzList.Add(peak.Mz);
                        intDic[mzBin] += peak.Intensity;
                    }
                    else
                    {
                        mzDic.Add(mzBin, new List<double> { peak.Mz });
                        intDic.Add(mzBin, peak.Intensity);
                    }
                }
            }
            var summedPeakList = new List<Peak>();
            foreach (var entry in mzDic)
            {
                var binNum = entry.Key;
                var mzList = entry.Value;
                var medianMz = mzList.Median();
                var intensity = intDic[binNum];
                summedPeakList.Add(new Peak(medianMz, intensity));
            }
            summedPeakList.Sort();

            var summedSpec = new SummedSpectrum(summedPeakList, repScanNum) { ScanNums = scanNums };
            return summedSpec;
        }

        /// <summary>
        /// Create a summed MS1 spectrum from the scans within <paramref name="elutionTimeTolerance"/> of <paramref name="scanNum"/>
        /// </summary>
        /// <param name="lcmsRun"></param>
        /// <param name="scanNum"></param>
        /// <param name="elutionTimeTolerance"></param>
        /// <returns>Summed spectrum</returns>
        public static Spectrum GetSummedMs1Spectrum(this ISpectrumAccessor lcmsRun, int scanNum, double elutionTimeTolerance)
        {
            var elutionTime = lcmsRun.GetElutionTime(scanNum);
            var minScanNum = scanNum;
            var curScanNum = scanNum;
            while (minScanNum >= lcmsRun.MinLcScan)
            {
                curScanNum = lcmsRun.GetPrevScanNum(curScanNum, 1);
                var curElutionTime = lcmsRun.GetElutionTime(curScanNum);
                if (elutionTime - curElutionTime < elutionTimeTolerance)
                {
                    minScanNum = curScanNum;
                }
                else
                {
                    break;
                }
            }
            var maxScanNum = scanNum;
            curScanNum = scanNum;
            while (maxScanNum <= lcmsRun.MaxLcScan)
            {
                curScanNum = lcmsRun.GetNextScanNum(curScanNum, 1);
                var curElutionTime = lcmsRun.GetElutionTime(curScanNum);
                if (curElutionTime - elutionTime < elutionTimeTolerance)
                {
                    maxScanNum = curScanNum;
                }
                else
                {
                    break;
                }
            }
            return lcmsRun.GetSummedMs1Spectrum(minScanNum, maxScanNum);
        }

        /// <summary>
        /// Create a summed MS1 spectrum from the scans in the supplied range
        /// </summary>
        /// <param name="lcmsRun"></param>
        /// <param name="minScanNum">min scan number, inclusive</param>
        /// <param name="maxScanNum">max scan number, inclusive</param>
        /// <returns>Summed spectrum</returns>
        public static SummedSpectrum GetSummedMs1Spectrum(this ISpectrumAccessor lcmsRun, int minScanNum, int maxScanNum)
        {
            if (minScanNum < lcmsRun.MinLcScan)
            {
                minScanNum = lcmsRun.MinLcScan;
            }

            if (maxScanNum > lcmsRun.MaxLcScan)
            {
                maxScanNum = lcmsRun.MaxLcScan;
            }

            var scanNums = new List<int>();
            for (var scanNum = minScanNum; scanNum <= maxScanNum; scanNum++)
            {
                if (lcmsRun.GetMsLevel(scanNum) != 1)
                {
                    continue;
                }

                scanNums.Add(scanNum);
            }

            return lcmsRun.GetSummedSpectrum(scanNums, minScanNum);
        }

        /// <summary>
        /// Get a summed MS2 spectrum from the dataset, with the provided limits
        /// </summary>
        /// <param name="lcmsRun"></param>
        /// <param name="monoIsotopicMass"></param>
        /// <param name="minScanNum">min scan number, inclusive</param>
        /// <param name="maxScanNum">max scan number, inclusive</param>
        /// <param name="minCharge">min charge, inclusive</param>
        /// <param name="maxCharge">max charge, inclusive</param>
        /// <param name="activationMethod"></param>
        /// <returns>Summed MS2 spectrum</returns>
        public static ProductSpectrum GetSummedMs2Spectrum(this ILcMsRun lcmsRun, double monoIsotopicMass,
            int minScanNum, int maxScanNum, int minCharge, int maxCharge, ActivationMethod activationMethod = ActivationMethod.Unknown)
        {
            var isoEnv = Averagine.GetIsotopomerEnvelope(monoIsotopicMass);
            var ms2ScanNumbers = new List<int>();
            for (var charge = minCharge; charge <= maxCharge; charge++)
            {
                var mostAbundantIsotopeMz = Ion.GetIsotopeMz(monoIsotopicMass, charge, isoEnv.MostAbundantIsotopeIndex);
                ms2ScanNumbers.AddRange(lcmsRun.GetFragmentationSpectraScanNums(mostAbundantIsotopeMz)
                    .Where(ms2ScanNum => ms2ScanNum >= minScanNum && ms2ScanNum <= maxScanNum &&
                        (activationMethod == ActivationMethod.Unknown ||
                        ((ProductSpectrum)lcmsRun.GetSpectrum(ms2ScanNum)).ActivationMethod == activationMethod)))
                        ;
            }
            var summedSpec = lcmsRun.GetSummedSpectrum(ms2ScanNumbers);
            return new ProductSpectrum(summedSpec.Peaks, 0) { ActivationMethod = activationMethod };
        }
    }
}
