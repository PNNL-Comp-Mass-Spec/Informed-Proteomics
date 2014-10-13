using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.Backend.MassSpecData
{
    public abstract class LcMsRun: ILcMsRun
    {
        public const int NumUniqueIsolationWindowThresholdForDia = 1000;
        public const double IsolationWindowBinningFactor = 10;
        //public const double DefaultSpectrumFilteringWindowSize = 100.0;

        public int MinLcScan { get; protected set; }
        public int MaxLcScan { get; protected set; }

        public int MinMsLevel { get; protected set; }
        public int MaxMsLevel { get; protected set; }

        public abstract double MinMs1Mz { get; }
        public abstract double MaxMs1Mz { get; }

        public bool IsDia
        {
            get
            {
                return IsDiaOrNull ?? (bool)(IsDiaOrNull = GetNumUniqueIsoWindows() < NumUniqueIsolationWindowThresholdForDia);
            }
        }

        public abstract Spectrum GetSpectrum(int scanNum);

        public Spectrum GetSummedMs1Spectrum(int scanNum, double elutionTimeTolerance)
        {
            var elutionTime = GetElutionTime(scanNum);
            var minScanNum = scanNum;
            var curScanNum = scanNum;
            while (minScanNum >= MinLcScan)
            {
                curScanNum = GetPrevScanNum(curScanNum, 1);
                var curElutionTime = GetElutionTime(curScanNum);
                if (elutionTime - curElutionTime < elutionTimeTolerance) minScanNum = curScanNum;
                else break;
            }
            var maxScanNum = scanNum;
            curScanNum = scanNum;
            while (maxScanNum <= MaxLcScan)
            {
                curScanNum = GetNextScanNum(curScanNum, 1);
                var curElutionTime = GetElutionTime(curScanNum);
                if (curElutionTime - elutionTime < elutionTimeTolerance) maxScanNum = curScanNum;
                else break;
            }
            return GetSummedMs1Spectrum(minScanNum, maxScanNum);
        }

        // minScanNum, maxScanNum: inclusive
        public SummedSpectrum GetSummedMs1Spectrum(int minScanNum, int maxScanNum)
        {
            if (minScanNum < MinLcScan) minScanNum = MinLcScan;
            if (maxScanNum > MaxLcScan) maxScanNum = MaxLcScan;

            var scanNums = new List<int>();
            for (var scanNum = minScanNum; scanNum <= maxScanNum; scanNum++)
            {
                if (GetMsLevel(scanNum) != 1) continue;
                scanNums.Add(scanNum);
            }

            return GetSummedSpectrum(scanNums, minScanNum);
        }

        public SummedSpectrum GetSummedSpectrum(IList<int> scanNums, int repScanNum = 0)
        {
            var mzComparer = new MzComparerWithBinning();

            var mzDic = new Dictionary<int, List<double>>();
            var intDic = new Dictionary<int, double>();

            foreach (var scanNum in scanNums)
            {
                var spec = GetSpectrum(scanNum);
                if (spec == null) continue;
                foreach (var peak in spec.Peaks)
                {
                    var mzBin = mzComparer.GetBinNumber(peak.Mz);

                    List<double> mzList;
                    if (mzDic.TryGetValue(mzBin, out mzList))
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

            var summedSpec = new SummedSpectrum(summedPeakList, repScanNum) {ScanNums = scanNums};
            return summedSpec;
        }

        // minScanNum, maxScanNum, minCharge, maxCharge: inclusive
        public ProductSpectrum GetSummedMs2Spectrum(double monoIsotopicMass,
            int minScanNum, int maxScanNum, int minCharge, int maxCharge, ActivationMethod activationMethod = ActivationMethod.Unknown)
        {
            var isoEnv = Averagine.GetIsotopomerEnvelope(monoIsotopicMass);
            var ms2ScanNums = new List<int>();
            for (var charge = minCharge; charge <= maxCharge; charge++)
            {
                var mostAbundantIsotopeMz = Ion.GetIsotopeMz(monoIsotopicMass, charge, isoEnv.MostAbundantIsotopeIndex);
                ms2ScanNums.AddRange(GetFragmentationSpectraScanNums(mostAbundantIsotopeMz)
                    .Where(ms2ScanNum => ms2ScanNum >= minScanNum && ms2ScanNum <= maxScanNum && 
                        (activationMethod == ActivationMethod.Unknown || 
                        ((ProductSpectrum) GetSpectrum(ms2ScanNum)).ActivationMethod == activationMethod)))
                        ;
            }
            var summedSpec = GetSummedSpectrum(ms2ScanNums);
            return new ProductSpectrum(summedSpec.Peaks, 0) {ActivationMethod = activationMethod};
        }

        public ProductSpectrum GetSummedMs2Spectrum(double monoIsotopicMass, ChargeScanRange range,
            ActivationMethod activationMethod = ActivationMethod.Unknown)
        {
            return GetSummedMs2Spectrum(monoIsotopicMass, range.MinScanNum, range.MaxScanNum, range.MinCharge,
                range.MaxCharge, activationMethod);
        }

        /// <summary>
        /// Gets the MS level of the specified scan
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>MS level</returns>
        public int GetMsLevel(int scanNum)
        {
            int msLevel;
            return ScanNumToMsLevel.TryGetValue(scanNum, out msLevel) ? msLevel : 0;
        }

        public double GetElutionTime(int scanNum)
        {
            double elutionTime;
            return ScanNumElutionTimeMap.TryGetValue(scanNum, out elutionTime) ? elutionTime : 0.0;
        }

        /// <summary>
        /// Gets the precursor scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns>precursor scan number or 0 for MS1</returns>
        public int GetPrecursorScanNum(int scanNum)
        {
            return _precursorScan[scanNum - MinLcScan];
        }

        /// <summary>
        /// Gets the next scan number whose ms level is smaller by 1
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>next scan number or MaxLc for MS1</returns>
        public int GetNextScanNum(int scanNum)
        {
            return _nextScan[scanNum - MinLcScan];
        }

        /// <summary>
        /// Gets the greatest scan number smaller than ms2ScanNum
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <param name="msLevel">MS level</param>
        /// <returns>previous scan number at the specified level</returns>
        public int GetPrevScanNum(int scanNum, int msLevel)
        {
            for (var curScanNum = scanNum - 1; curScanNum >= MinLcScan; curScanNum--)
            {
                if (GetMsLevel(curScanNum) == msLevel) return curScanNum;
            }
            return MinLcScan - 1;
        }

        /// <summary>
        /// Gets the smallest scan number larger than ms2ScanNum
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <param name="msLevel">MS level</param>
        /// <returns>next scan number at the specified level</returns>
        public int GetNextScanNum(int scanNum, int msLevel)
        {
            for (var curScanNum = scanNum + 1; curScanNum <= MaxLcScan; curScanNum++)
            {
                if (GetMsLevel(curScanNum) == msLevel)
                {
                    return curScanNum;
                }
            }
            return MaxLcScan + 1;
        }

        /// <summary>
        /// Gets the scan numbers of the specified msLevel
        /// </summary>
        /// <param name="msLevel">MS level</param>
        /// <returns>scan numbers of the specified msLevel</returns>
        public IList<int> GetScanNumbers(int msLevel)
        {
            var scanNumbers = new List<int>();
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                if (GetMsLevel(scanNum) == msLevel) scanNumbers.Add(scanNum);
            }
            return scanNumbers;
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z (using only MS1 spectra)
        /// XicPoint is created for every MS1 scan.
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <returns>XIC as an Xic object</returns>
        public Xic GetFullPrecursorIonExtractedIonChromatogram(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetFullPrecursorIonExtractedIonChromatogram(minMz, maxMz);
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS1 spectra)
        /// XicPoint is created for every MS1 scan.
        /// </summary>
        /// <param name="minMz">min m/z</param>
        /// <param name="maxMz">max m/z</param>
        /// <returns>XIC as an Xic object</returns>
        public Xic GetFullPrecursorIonExtractedIonChromatogram(double minMz, double maxMz)
        {
            var xic = GetPrecursorExtractedIonChromatogram(minMz, maxMz);

            var hasXicPoint = new bool[MaxLcScan - MinLcScan + 1];
            foreach (var xicPoint in xic) hasXicPoint[xicPoint.ScanNum - MinLcScan] = true;

            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                if (GetMsLevel(scanNum) == 1 && !hasXicPoint[scanNum - MinLcScan]) xic.Add(new XicPoint(scanNum, 0, 0));
            }
            xic.Sort();
            return xic;
        }

        private int[] _ms1ScanVector;

        public int[] GetMs1ScanVector()
        {
            return _ms1ScanVector ?? (_ms1ScanVector = GetScanNumbers(1).ToArray());
        }

        public double[] GetFullPrecursorIonExtractedIonChromatogramVector(double minMz, double maxMz)
        {
            var xic = GetPrecursorExtractedIonChromatogram(minMz, maxMz);

            var xicIntensityVector = new double[MaxLcScan - MinLcScan + 1];
            foreach (var xicPoint in xic) xicIntensityVector[xicPoint.ScanNum - MinLcScan] = xicPoint.Intensity;
            var intVector = GetMs1ScanVector().Select(s => xicIntensityVector[s - MinLcScan]).ToArray();
            return intVector;
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS2 spectra)
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <param name="precursorIonMz">precursor m/z of the precursor ion</param>
        /// <returns>XIC as an Xic object</returns>
        public Xic GetFullProductExtractedIonChromatogram(double mz, Tolerance tolerance, double precursorIonMz)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetFullProductExtractedIonChromatogram(minMz, maxMz, precursorIonMz);
        }

        public abstract Xic GetFullProductExtractedIonChromatogram(double minMz, double maxMz, double precursorMz);

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z (using only MS1 spectra)
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <returns>XIC as an Xic object</returns>
        public Xic GetPrecursorExtractedIonChromatogram(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetPrecursorExtractedIonChromatogram(minMz, maxMz);
        }

        public abstract Xic GetPrecursorExtractedIonChromatogram(double minMz, double maxMz);

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS1 spectra)
        /// Only XicPeaks around the targetScanNum are returned 
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <param name="targetScanNum">target scan number to generate xic</param>
        /// <param name="maxNumConsecutiveScansWithoutPeak">maximum number of consecutive scans with a peak</param>
        /// <returns>XIC around targetScanNum</returns>
        public Xic GetPrecursorExtractedIonChromatogram(double mz, Tolerance tolerance, int targetScanNum, int maxNumConsecutiveScansWithoutPeak = 3)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            if (targetScanNum < 0) return GetPrecursorExtractedIonChromatogram(minMz, maxMz);
            return GetPrecursorExtractedIonChromatogram(minMz, maxMz, targetScanNum);
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS1 spectra)
        /// Only XicPeaks around the targetScanNum are returned 
        /// </summary>
        /// <param name="minMz">min m/z</param>
        /// <param name="maxMz">max m/z</param>
        /// <param name="targetScanNum">target scan number to generate xic</param>
        /// <param name="tolerance">max number of consecutive scans without a peak</param>
        /// <returns>XIC around targetScanNum</returns>
        public Xic GetPrecursorExtractedIonChromatogram(double minMz, double maxMz, int targetScanNum, int tolerance = 3)
        {
            if (GetMsLevel(targetScanNum) > 1)
                targetScanNum = GetPrecursorScanNum(targetScanNum);
            var xic = GetPrecursorExtractedIonChromatogram(minMz, maxMz);
            return GetTrimmedXic(xic, targetScanNum);
        }

        /// <summary>
        /// Get a segment of Xic containing the targetScanNum
        /// </summary>
        /// <param name="xic">xic to be trimmed</param>
        /// <param name="targetScanNum">target scan number to generate xic</param>
        /// <param name="tolerance">number of scans that can be tolerated</param>
        /// <returns>Trimmed XIC around targetScanNum</returns>
        public Xic GetTrimmedXic(Xic xic, int targetScanNum, int tolerance = 3)
        {
            var index = xic.BinarySearch(new XicPoint(targetScanNum, 0, 0));
            if (index < 0) index = ~index;

            var xicSegment = new Xic();

            var curScanNum = targetScanNum;

            // go down
            var i = index - 1;
            while (i >= 0 && i < xic.Count)
            {
                var xicPeak = xic[i];
                // check whether there's no MS1 scan in between
                var isConsecutive = true;
                var numMissingScans = 0;
                for (var scanNum = xicPeak.ScanNum + 1; scanNum < curScanNum; scanNum++)
                {
                    if (GetMsLevel(scanNum) == 1) ++numMissingScans;
                    if (numMissingScans > tolerance)
                    {
                        isConsecutive = false;
                        break;
                    }
                }
                if (isConsecutive) xicSegment.Add(xicPeak);
                else break;
                curScanNum = xicPeak.ScanNum;
                --i;
            }

            // go up
            i = index;
            curScanNum = targetScanNum;
            while (i >= 0 && i < xic.Count)
            {
                var xicPeak = xic[i];
                // check whether there's no MS1 scan in between
                var numMissingScans = 0;
                var isConsecutive = true;
                for (var scanNum = curScanNum + 1; scanNum < xicPeak.ScanNum; scanNum++)
                {
                    if (GetMsLevel(scanNum) == 1) ++numMissingScans;
                    if (numMissingScans > tolerance)
                    {
                        isConsecutive = false;
                        break;
                    }
                }
                if (isConsecutive) xicSegment.Add(xicPeak);
                else break;
                curScanNum = xicPeak.ScanNum;
                ++i;
            }

            xicSegment.Sort();
            return xicSegment;
        }

        public static void WriteSpectrum(Spectrum spec, BinaryWriter writer)
        {
            // scan number: 4
            writer.Write(spec.ScanNum);

            // ms level: 1
            writer.Write(Convert.ToByte(spec.MsLevel));

            // elution time: 4
            writer.Write(spec.ElutionTime);

            var productSpec = spec as ProductSpectrum;
            if (productSpec != null)    // product spectrum
            {
                var isolationWindow = productSpec.IsolationWindow;
                // precursor mass: 8
                writer.Write(isolationWindow.MonoisotopicMz ?? 0.0);
                // precursor charge: 4
                writer.Write(isolationWindow.Charge ?? 0);
                // Activation method: 1
                writer.Write((byte)productSpec.ActivationMethod);
                // Isolation window target m/z: 8
                writer.Write(productSpec.IsolationWindow.IsolationWindowTargetMz);
                // Isolation window lower offset: 8
                writer.Write(productSpec.IsolationWindow.IsolationWindowLowerOffset);
                // Isolation window uppoer offset: 8
                writer.Write(productSpec.IsolationWindow.IsolationWindowUpperOffset);
            }
            // Number of peaks: 4
            writer.Write(spec.Peaks.Length);
            foreach (var peak in spec.Peaks)
            {
                // m/z: 8
                writer.Write(peak.Mz);
                // intensity: 4
                writer.Write(Convert.ToSingle(peak.Intensity));
            }
        }


        public int GetNumUniqueIsoWindows()
        {
            var isoWindowSet = new HashSet<IsolationWindow>();
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                var productSpec = GetSpectrum(scanNum) as ProductSpectrum;
                if (productSpec == null) continue;
                isoWindowSet.Add(productSpec.IsolationWindow);
            }
            return isoWindowSet.Count;
        }

        public double GetMinIsolationWindowWidth()
        {
            var minWidth = Double.MaxValue;
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                var productSpec = GetSpectrum(scanNum) as ProductSpectrum;
                if (productSpec == null) continue;
                if (productSpec.IsolationWindow.Width < minWidth) minWidth = productSpec.IsolationWindow.Width;
            }
            return minWidth;
        }

        /// <summary>
        /// Gets scan numbers of the fragmentation spectra whose isolation window contains the precursor ion specified
        /// </summary>
        /// <param name="precursorIon"></param>
        /// <returns>scan numbers of fragmentation spectra</returns>
        public int[] GetFragmentationSpectraScanNums(Ion precursorIon)
        {
            var mostAbundantIsotopeMz = precursorIon.GetMostAbundantIsotopeMz();
            return GetFragmentationSpectraScanNums(mostAbundantIsotopeMz);
        }

        /// <summary>
        /// Gets scan numbers of the fragmentation spectra whose isolation window contains the precursor ion specified
        /// </summary>
        /// <param name="mostAbundantIsotopeMz"></param>
        /// <returns>scan numbers of fragmentation spectra</returns>
        public int[] GetFragmentationSpectraScanNums(double mostAbundantIsotopeMz)
        {
            var targetIsoBin = (int)Math.Round(mostAbundantIsotopeMz * IsolationWindowBinningFactor);
            int[] scanNums;
            return IsolationMzBinToScanNums.TryGetValue(targetIsoBin, out scanNums) ? scanNums : new int[0];
        }

        //public IEnumerable<int> GetMs2ScansForPrecursorMz(double precursorMz)
        //{
        //    return
        //        from ms2ScanNum in GetScanNumbers(2)
        //        let productSpec = GetSpectrum(ms2ScanNum) as ProductSpectrum
        //        where productSpec != null
        //        where productSpec.IsolationWindow.Contains(precursorMz)
        //        select ms2ScanNum;
        //}

        // Fields to be defined in a child
        protected Dictionary<int, int[]> IsolationMzBinToScanNums;
        protected Dictionary<int, int> ScanNumToMsLevel;
        protected Dictionary<int, double> ScanNumElutionTimeMap;
        protected bool? IsDiaOrNull;

        private int[] _precursorScan;
        private int[] _nextScan;

        protected void CreatePrecursorNextScanMap()
        {
            _precursorScan = new int[MaxLcScan - MinLcScan + 1];
            _nextScan = new int[MaxLcScan - MinLcScan + 1];

            var precursorMap = new Dictionary<int, int>();
            var nextScanMap = new Dictionary<int, int>();
            for (var msLevel = MinMsLevel; msLevel <= MaxMsLevel; msLevel++)
            {
                precursorMap[msLevel] = 0;
                nextScanMap[msLevel] = MaxLcScan + 1;
            }

            var prevMsLevel = int.MinValue;
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                var index = scanNum - MinLcScan;

                var msLevel = GetMsLevel(scanNum);
                // determine precursor scan
                if (msLevel == prevMsLevel)
                {
                    _precursorScan[index] = _precursorScan[index - 1];
                }
                else if (msLevel > prevMsLevel)
                {
                    _precursorScan[index] = scanNum - 1;
                    precursorMap[msLevel] = scanNum - 1;
                }
                else // msLevel < prevMsLevel
                {
                    _precursorScan[index] = precursorMap[msLevel];
                }
                prevMsLevel = msLevel;
            }

            var nextMsLevel = int.MinValue;
            for (var scanNum = MaxLcScan; scanNum >= MinLcScan; scanNum--)
            {
                var index = scanNum - MinLcScan;
                var msLevel = GetMsLevel(scanNum);

                // determine next scan
                if (msLevel == nextMsLevel)
                {
                    _nextScan[index] = _nextScan[index + 1];
                }
                else if (msLevel > nextMsLevel)
                {
                    _nextScan[index] = scanNum + 1;
                    nextScanMap[msLevel] = scanNum + 1;
                }
                else // msLevel < prevMsLevel
                {
                    _nextScan[index] = nextScanMap[msLevel];
                }
                nextMsLevel = msLevel;
            }
        }
    }
}
