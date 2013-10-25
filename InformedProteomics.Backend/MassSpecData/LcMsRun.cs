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
            if (File.Exists(pbfFilePath))
            {
                return new LcMsRun(new PbfReader(pbfFilePath),
                    precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold);
            }

            LcMsRun run;
            if (dataType == MassSpecDataType.XCaliburRun)
            {
                run = new LcMsRun(new XCaliburReader(specFilePath), 
                    precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold);
            }
            else run = null;

            if(run != null) run.WriteTo(pbfFilePath);

            return run;
        }

        public const double IsolationWindowBinningFactor = 10;

        public LcMsRun(IMassSpecDataReader massSpecDataReader) : this(massSpecDataReader, 0.0, 0.0)
        {
            
        }

        public LcMsRun(IMassSpecDataReader massSpecDataReader, double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold)
        {
            _scanNumSpecMap = new Dictionary<int, Spectrum>();
            _ms1PeakList = new List<LcMsPeak>();
            _msLevel = new Dictionary<int, int>();

            var isolationMzBinToScanNums = new Dictionary<int, List<int>>();

            // Read all spectra
            var minScanNum = int.MaxValue;
            var maxScanNum = int.MinValue;

            foreach (var spec in massSpecDataReader.ReadAllSpectra())
            {
                //Console.WriteLine("Reading Scan {0}", spec.ScanNum);
                _msLevel[spec.ScanNum] = spec.MsLevel;
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

                _isolationMzBinToScanNums = new Dictionary<int, IEnumerable<int>>();
                foreach (var entry in isolationMzBinToScanNums)
                {
                    var binNum = entry.Key;
                    var scanNumList = entry.Value;
                    _isolationMzBinToScanNums[binNum] = scanNumList;
                }

                if (spec.ScanNum < minScanNum) minScanNum = spec.ScanNum;
                if (spec.ScanNum > maxScanNum) maxScanNum = spec.ScanNum;
            }

           _ms1PeakList.Sort();

            // Read MS levels and precursor information

           MinLcScan = minScanNum;
           MaxLcScan = maxScanNum;

            _precursorScan = new int[MaxLcScan - MinLcScan + 1];
            var precursorMap = new Dictionary<int, int>();
            precursorMap[1] = 0;

            var prevMsLevel = int.MinValue;
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                var index = scanNum - MinLcScan;
                var msLevel = _msLevel[scanNum];
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

            massSpecDataReader.Close();
        }

        public int MinLcScan { get; private set; }
        public int MaxLcScan { get; private set; }

        /// <summary>
        /// Gets scan numbers of the fragmentation spectra whose isolation window contains the precursor ion specified
        /// </summary>
        /// <param name="precursorIon"></param>
        /// <returns>scan numbers of fragmentation spectra</returns>
        public IEnumerable<int> GetFragmentationSpectraScanNums(Ion precursorIon)
        {
            var precursorBaseMz = precursorIon.GetBaseIsotopeMz();
            var targetIsoBin = (int)Math.Round(precursorBaseMz*IsolationWindowBinningFactor);
            IEnumerable<int> scanNums;
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

        /// <summary>
        /// Gets the extracted ion chromatogram corresponding to an isotope of another XIC
        /// </summary>
        /// <param name="xic">base XIC</param>
        /// <param name="mzDifference">m/z difference</param>
        /// <param name="tolerance">tolerance</param>
        /// <returns>XIC corresponding to an isotope of the input XIC</returns>
        public Xic GetIsotopeExtractedIonChromatogram(Xic xic, double mzDifference, Tolerance tolerance)
        {
            var isotopeXic = new Xic();
            foreach (var xicPeak in xic)
            {
                var spec = _scanNumSpecMap[xicPeak.ScanNum];
                var peak = spec.FindPeak(xicPeak.Mz + mzDifference, tolerance);
                if (peak != null) isotopeXic.Add(new XicPeak(xicPeak.ScanNum, peak.Mz, peak.Intensity));
            }
            return isotopeXic;
        }

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
                xic.Add(new XicPeak(peak.ScanNum, peak.Mz, peak.Intensity));
                --i;
            }

            // go up
            i = index;
            while (i >= 0 && i < _ms1PeakList.Count)
            {
                var peak = _ms1PeakList[i];
                if (peak.Mz >= maxMz) break;
                xic.Add(new XicPeak(peak.ScanNum, peak.Mz, peak.Intensity));
                ++i;
            }

            if (xic.Count == 0) return xic;

            xic.Sort();
            // select one best peak for each scan
            var newXic = new Xic();

            var prevScanNum = -1;
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

            var index = xic.BinarySearch(new XicPeak(targetScanNum, 0, 0));
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
                if(peak != null) productXic.Add(new XicPeak(scanNum, peak.Mz, peak.Intensity));
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
        /// Gets the MS level of the specified scan
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>MS level</returns>
        public int GetMsLevel(int scanNum)
        {
            int msLevel;
            if (_msLevel.TryGetValue(scanNum, out msLevel))
            {
                return msLevel;
            }
            Console.WriteLine("Strange: {0}" + scanNum);
            return -1;
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
        /// Gets the greatest scan number smaller than scanNum
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
        /// Gets the smallest scan number larger than scanNum
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
        /// Gets MS/MS spectra whose isolation windows contain the most abundant peak of the precursorIon
        /// </summary>
        /// <param name="precursorIon"></param>
        /// <returns></returns>
        public IEnumerable<ProductSpectrum> GetMatchingMs2Spectra(Ion precursorIon)
        {
            //var baseIsotopeIndex = precursorIon.Composition.GetMostAbundantIsotopeZeroBasedIndex();
            //var baseXic = GetExtractedIonChromatogram()
            //var relativeIsotopeIntensityThreshold = 0.5;
            //precursorIon.GetIsotopes(relativeIsotopeIntensityThreshold);
            return null;
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

        private readonly Dictionary<int, int> _msLevel;
        private readonly int[] _precursorScan;

        private readonly List<LcMsPeak> _ms1PeakList;
        private readonly Dictionary<int, Spectrum> _scanNumSpecMap;  // scan number -> spectrum
        private readonly Dictionary<int, IEnumerable<int>> _isolationMzBinToScanNums;

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
