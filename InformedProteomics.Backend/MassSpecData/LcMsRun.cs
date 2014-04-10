using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public class LcMsRun
    {
        public static LcMsRun GetLcMsRun(string specFilePath, MassSpecDataType dataType)
        {
            return GetLcMsRun(specFilePath, dataType, 0.0, 0.0);
        }

        public static LcMsRun GetLcMsRun(string specFilePath, MassSpecDataType dataType, 
            double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold)
        {
            var pbfFilePath = Path.ChangeExtension(specFilePath, ".pbf");

            if (!File.Exists(pbfFilePath))
            {
                LcMsRun run;
                if (dataType == MassSpecDataType.XCaliburRun)
                {
                    run = new LcMsRun(new XCaliburReader(specFilePath),
                        0, 0);
                }
                else run = null;
                if (run != null) run.WriteTo(pbfFilePath);
                else throw new Exception("Unsupported raw file format!");
            }

            return new LcMsRun(new PbfReader(pbfFilePath),
                precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold);
        }

        public const double IsolationWindowBinningFactor = 10;

        private LcMsRun(IMassSpecDataReader massSpecDataReader, double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold)
        {
            _scanNumSpecMap = new Dictionary<int, Spectrum>();
            _ms1PeakList = new List<LcMsPeak>();
            var msLevelMap = new Dictionary<int, int>();

            var isolationMzBinToScanNums = new Dictionary<int, List<int>>();

            // Read all spectra
            var minScanNum = int.MaxValue;
            var maxScanNum = int.MinValue;
            var minMsLevel = int.MaxValue;
            var maxMsLevel = int.MinValue;

            foreach (var spec in massSpecDataReader.ReadAllSpectra())
            {
                //Console.WriteLine("Reading Scan {0}", spec.ScanNum);
                msLevelMap[spec.ScanNum] = spec.MsLevel;
                if (spec.MsLevel == 1)
                {
                    if (precursorSignalToNoiseRatioThreshold > 0.0) spec.FilterNoise(precursorSignalToNoiseRatioThreshold);
                    foreach (var peak in spec.Peaks)
                    {
                        _ms1PeakList.Add(new LcMsPeak(peak.Mz, peak.Intensity, spec.ScanNum));
                    }
                    _scanNumSpecMap.Add(spec.ScanNum, spec);
                }
                else if(spec.MsLevel == 2)
                {
                    var productSpec = spec as ProductSpectrum;

                    if (productSpec != null)
                    {
                        if (productSignalToNoiseRatioThreshold > 0.0) productSpec.FilterNoise(productSignalToNoiseRatioThreshold);
                        var isolationInfo = productSpec.IsolationWindow;
                        var minBinNum = (int)(isolationInfo.MinMz*IsolationWindowBinningFactor);
                        var maxBinNum = (int)((isolationInfo.MaxMz+0.49999999)*IsolationWindowBinningFactor);
                        for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
                        {
                            List<int> scanNumList;
                            if (!isolationMzBinToScanNums.TryGetValue(binNum, out scanNumList))
                            {
                                scanNumList = new List<int>();
                                isolationMzBinToScanNums[binNum] = scanNumList;
                            }
                            scanNumList.Add(productSpec.ScanNum);
                        }
                        _scanNumSpecMap.Add(spec.ScanNum, productSpec);
                    }
                }

                if (spec.ScanNum < minScanNum) minScanNum = spec.ScanNum;
                if (spec.ScanNum > maxScanNum) maxScanNum = spec.ScanNum;

                if (spec.MsLevel < minMsLevel) minMsLevel = spec.MsLevel;
                if (spec.MsLevel > maxMsLevel) maxMsLevel = spec.MsLevel;

            }

            _isolationMzBinToScanNums = new Dictionary<int, int[]>();
            foreach (var entry in isolationMzBinToScanNums)
            {
                var binNum = entry.Key;
                entry.Value.Sort();
                var scanNumList = entry.Value.ToArray();
                _isolationMzBinToScanNums[binNum] = scanNumList;
            }

           _ms1PeakList.Sort();

            // Read MS levels and precursor information

           MinLcScan = minScanNum;
           MaxLcScan = maxScanNum;

            MinMsLevel = minMsLevel;
            MaxMsLevel = maxMsLevel;

            var precursorMap = new Dictionary<int, int>();
            var nextScanMap = new Dictionary<int, int>();

            for (var msLevel = MinMsLevel; msLevel <= maxMsLevel; msLevel++)
            {
                precursorMap[msLevel] = 0;
                nextScanMap[msLevel] = MaxLcScan + 1;
            }

            _msLevel = new int[MaxLcScan - MinLcScan + 1];
            _precursorScan = new int[MaxLcScan - MinLcScan + 1];

            var prevMsLevel = int.MinValue;
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                var index = scanNum - MinLcScan;
                int msLevel;
                if (!msLevelMap.TryGetValue(scanNum, out msLevel)) msLevel = 0;
                _msLevel[index] = msLevel;

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

            _nextScan = new int[MaxLcScan - MinLcScan + 1];
            var nextMsLevel = int.MinValue;
            for (var scanNum = MaxLcScan; scanNum >= MinLcScan; scanNum--)
            {
                var index = scanNum - MinLcScan;
                var msLevel = _msLevel[index];

                // determine precursor scan
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

        public int MinLcScan { get; private set; }
        public int MaxLcScan { get; private set; }

        public int MinMsLevel { get; private set; }
        public int MaxMsLevel { get; private set; }

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
            return _isolationMzBinToScanNums.TryGetValue(targetIsoBin, out scanNums) ? scanNums : new int[0];
        }

        /// <summary>
        /// Gets the spectrum of the specified scan number
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>spectrum</returns>
        public Spectrum GetSpectrum(int scanNum)
        {
            Spectrum spec;
            return _scanNumSpecMap.TryGetValue(scanNum, out spec) ? spec : null;
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z (using only MS1 spectra)
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <returns>XIC as an Xic object</returns>
        public Xic GetExtractedIonChromatogram(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetExtractedIonChromatogram(minMz, maxMz);
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS1 spectra)
        /// Only XicPeaks around the targetScanNum are returned 
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <param name="targetScanNum">target scan number to generate xic</param>
        /// <returns>XIC around targetScanNum</returns>
        public Xic GetExtractedIonChromatogram(double mz, Tolerance tolerance, int targetScanNum)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            if (targetScanNum < 0) return GetExtractedIonChromatogram(minMz, maxMz);
            return GetExtractedIonChromatogram(minMz, maxMz, targetScanNum);
        }

        ///// <summary>
        ///// Gets the extracted ion chromatogram corresponding to an isotope of another XIC
        ///// </summary>
        ///// <param name="xic">base XIC</param>
        ///// <param name="mzDifference">m/z difference</param>
        ///// <param name="tolerance">tolerance</param>
        ///// <returns>XIC corresponding to an isotope of the input XIC</returns>
        //public Xic GetIsotopeExtractedIonChromatogram(Xic xic, double mzDifference, Tolerance tolerance)
        //{
        //    var isotopeXic = new Xic();
        //    foreach (var xicPoint in xic)
        //    {
        //        var spec = _scanNumSpecMap[xicPoint.ScanNum];
        //        var peak = spec.FindPeak(xicPoint.Mz + mzDifference, tolerance);
        //        if (peak != null) isotopeXic.Add(new XicPoint(xicPoint.ScanNum, peak.Mz, peak.Intensity));
        //    }
        //    return isotopeXic;
        //}

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS1 spectra)
        /// </summary>
        /// <param name="minMz">min m/z</param>
        /// <param name="maxMz">max m/z</param>
        /// <returns>XIC as an Xic object</returns>
        public Xic GetExtractedIonChromatogram(double minMz, double maxMz)
        {
            var xic = new Xic();

            var index = _ms1PeakList.BinarySearch(new LcMsPeak((minMz + maxMz) / 2, 0, 0));
            if (index < 0) index = ~index;

            // go down
            var i = index - 1;
            while (i >= 0 && i < _ms1PeakList.Count)
            {
                var peak = _ms1PeakList[i];
                if (peak.Mz <= minMz) break;
                xic.Add(new XicPoint(peak.ScanNum, peak.Intensity));
                --i;
            }

            // go up
            i = index;
            while (i >= 0 && i < _ms1PeakList.Count)
            {
                var peak = _ms1PeakList[i];
                if (peak.Mz >= maxMz) break;
                xic.Add(new XicPoint(peak.ScanNum, peak.Intensity));
                ++i;
            }

            if (xic.Count == 0) return xic;

            xic.Sort();
            // select one best peak for each scan
            var newXic = new Xic();

            var prevScanNum = xic[0].ScanNum;
            var bestPeak = xic[0];
            for (i = 1; i < xic.Count; i++)
            {
                var xicPeak = xic[i];
                if (xicPeak.ScanNum > prevScanNum)
                {
                    newXic.Add(bestPeak);
                    bestPeak = xicPeak;
                }
                else
                {
                    if (xicPeak.Intensity > bestPeak.Intensity)
                    {
                        bestPeak = xicPeak;
                    }
                }
                prevScanNum = xicPeak.ScanNum;
            }
            newXic.Add(bestPeak);
            return newXic;
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z (using only MS1 spectra)
        /// XicPoint is created for every MS1 scan.
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <returns>XIC as an Xic object</returns>
        public Xic GetFullExtractedIonChromatogram(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetFullExtractedIonChromatogram(minMz, maxMz);
        }


        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS1 spectra)
        /// XicPoint is created for every MS1 scan.
        /// </summary>
        /// <param name="minMz">min m/z</param>
        /// <param name="maxMz">max m/z</param>
        /// <returns>XIC as an Xic object</returns>
        public Xic GetFullExtractedIonChromatogram(double minMz, double maxMz)
        {
            var xic = GetExtractedIonChromatogram(minMz, maxMz);

            var hasXicPoint = new bool[MaxLcScan - MinLcScan + 1];
            foreach (var xicPoint in xic) hasXicPoint[xicPoint.ScanNum - MinLcScan] = true;

            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                if(GetMsLevel(scanNum) == 1 && !hasXicPoint[scanNum-MinLcScan]) xic.Add(new XicPoint(scanNum, 0));
            }
            xic.Sort();
            return xic;
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS1 spectra)
        /// Only XicPeaks around the targetScanNum are returned 
        /// </summary>
        /// <param name="minMz">min m/z</param>
        /// <param name="maxMz">max m/z</param>
        /// <param name="targetScanNum">target scan number to generate xic</param>
        /// <returns>XIC around targetScanNum</returns>
        public Xic GetExtractedIonChromatogram(double minMz, double maxMz, int targetScanNum)
        {
            if (GetMsLevel(targetScanNum) > 1)
                targetScanNum = GetPrecursorScanNum(targetScanNum);
            var xic = GetExtractedIonChromatogram(minMz, maxMz);
            return GetTrimmedXic(xic, targetScanNum);
        }

        /// <summary>
        /// Get a segment of Xic containing the targetScanNum
        /// </summary>
        /// <param name="xic">xic to be trimmed</param>
        /// <param name="targetScanNum">target scan number to generate xic</param>
        /// <returns>Trimmed XIC around targetScanNum</returns>
        public Xic GetTrimmedXic(Xic xic, int targetScanNum)
        {
            //TODO: this tolerance value is not optimal for all data (revisit required)
            const int tolerance = 3;

            var index = xic.BinarySearch(new XicPoint(targetScanNum, 0));
            if(index < 0) index = ~index;

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
                    if(numMissingScans > tolerance)
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
                for (var scanNum = curScanNum+1; scanNum < xicPeak.ScanNum; scanNum++)
                {
                    if (GetMsLevel(scanNum) == 1) ++numMissingScans;
                    if(numMissingScans > tolerance)
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

        public Xic GetProductIonChromatogram(double productIonMz, double precursorIonMz, Tolerance tolerance, int minScanNum, int maxScanNum)
        {
            var productXic = new Xic();

            for (var scanNum = minScanNum; scanNum <= maxScanNum; scanNum++)
            {
                if (_msLevel[scanNum] == 1) continue;

                var productSpec = _scanNumSpecMap[scanNum] as ProductSpectrum;
                if (productSpec == null) continue;

                if (!productSpec.IsolationWindow.Contains(precursorIonMz)) continue;

                var peak = productSpec.FindPeak(productIonMz, tolerance);
                if(peak != null) productXic.Add(new XicPoint(scanNum, peak.Intensity));
            }

            return productXic;
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
        /// Gets the MS level of the specified scan
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>MS level</returns>
        public int GetMsLevel(int scanNum)
        {
            return _msLevel[scanNum - MinLcScan];
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

        private class ComparerWithTolerance : IComparer<double>
        {
            private readonly double _tolerance;
            public ComparerWithTolerance(double tolerance)
            {
                _tolerance = tolerance;
            }
            public int Compare(double x, double y)
            {
                if (Math.Abs((x - y)/x*1E6) <= _tolerance) return 0;
                return x.CompareTo(y);
            }
        }
        private readonly ComparerWithTolerance _comparer = new ComparerWithTolerance(10);

        public bool CheckMs1Signature(Ion precursorIon, int ms2ScanNum, Tolerance tolerance)
        {
            if (GetMsLevel(ms2ScanNum) != 2) return false;

            var prevMs1ScanNum = GetPrecursorScanNum(ms2ScanNum);
            var nextMs1ScanNum = GetNextScanNum(ms2ScanNum);

            if (_ms1Features != null)
            {
                var mostAbundantIsotopeMz = precursorIon.GetMostAbundantIsotopeMz();
                if (prevMs1ScanNum > 0 &&
                    _ms1Features[prevMs1ScanNum][precursorIon.Charge - _ms1FeatureMinCharge]
                    .BinarySearch(mostAbundantIsotopeMz, _comparer) >= 0) return true;
                if (nextMs1ScanNum <= MaxLcScan
                    && _ms1Features[nextMs1ScanNum][precursorIon.Charge - _ms1FeatureMinCharge]
                    .BinarySearch(mostAbundantIsotopeMz, _comparer) >= 0) return true;
                return false;
            }

            var isotopes = precursorIon.GetTop3Isotopes();
            var prevMs1Spec = GetSpectrum(prevMs1ScanNum);
            if (prevMs1Spec != null)
            {
                //if (prevMs1Spec.ContainsIon(precursorIon, ToleranceForBaseXic, 0.9)) return true;
                var isPrevIsotopeValid = false;
                foreach (var isotope in isotopes)
                {
                    var mz = precursorIon.GetIsotopeMz(isotope.Index);
                    if (prevMs1Spec.FindPeak(mz, tolerance) != null)  // match
                    {
                        if (isPrevIsotopeValid) return true;
                        isPrevIsotopeValid = true;
                    }
                }
            }

            var nextMs1Spec = GetSpectrum(nextMs1ScanNum);
            if (nextMs1Spec != null)
            {
                var isPrevIsotopeValid = false;
                foreach (var isotope in isotopes)
                {
                    var mz = precursorIon.GetIsotopeMz(isotope.Index);
                    if (nextMs1Spec.FindPeak(mz, tolerance) != null)  // match
                    {
                        if (isPrevIsotopeValid) return true;
                        isPrevIsotopeValid = true;
                    }
                }
            }                
            return false;                
        }


        // These writing/reading methods are added to speed up the testing
        public void WriteTo(string outputFilePath)
        {
            using (var writer = new BinaryWriter(File.Open(outputFilePath, FileMode.Create)))
            {
                for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++ )
                {
                    var spec = GetSpectrum(scanNum);
                    PbfReader.WriteSpectrumAsPbf(spec, writer);
                }
            }
        }

        public void ComputeMs1Features(int minCharge, int maxCharge, Tolerance tolerance)
        {
            _ms1Features = new Dictionary<int, List<double>[]>();

            foreach (var scanNum in GetScanNumbers(msLevel: 1))
            {
                var spec = GetSpectrum(scanNum);
                var ms1FeaturesForScan = new List<double>[maxCharge - minCharge + 1];
                for (var charge = minCharge; charge <= maxCharge; charge++)
                {
                    ms1FeaturesForScan[charge - minCharge] = new List<double>();
                    foreach (var p in spec.Peaks)
                    {
                        var mz = p.Mz;
                        if (spec.FindPeak(mz + Constants.C13MinusC12 / charge, tolerance) != null)
                        {
                            ms1FeaturesForScan[charge - minCharge].Add(mz);
                        }
                    }
                }
                _ms1Features[scanNum] = ms1FeaturesForScan;
            }
            _ms1FeatureMinCharge = minCharge;
        }

        //private readonly Dictionary<int, int> _msLevel;
        private readonly int[] _msLevel;
        private readonly int[] _precursorScan;
        private readonly int[] _nextScan;

        private readonly List<LcMsPeak> _ms1PeakList;
        private readonly Dictionary<int, Spectrum> _scanNumSpecMap;  // scan number -> spectrum
        private readonly Dictionary<int, int[]> _isolationMzBinToScanNums;
        private Dictionary<int, List<double>[]> _ms1Features;
        private int _ms1FeatureMinCharge;
    }

    class LcMsPeak: Peak
    {
        public LcMsPeak(double mz, double intensity, int scanNum) : base(mz, intensity)
        {
            ScanNum = scanNum;
        }

        public int ScanNum { get; private set; }
    }
}
