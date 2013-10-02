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
            var pgfFilePath = Path.ChangeExtension(specFilePath, ".pgf");
            if(File.Exists(pgfFilePath)) return new LcMsRun(new PgfReader(pgfFilePath));

            LcMsRun run;
            if (dataType == MassSpecDataType.XCaliburRun) run = new LcMsRun(new XCaliburReader(specFilePath));
            else run = null;

            if(run != null) run.WriteTo(pgfFilePath);

            return run;
        }

        public LcMsRun(IMassSpecDataReader massSpecDataReader)
        {
            _scanNumSpecMap = new Dictionary<int, Spectrum>();
            _ms1PeakList = new List<LcMsPeak>();
            _msLevel = new Dictionary<int, int>();

            // Read all spectra
            var minScanNum = int.MaxValue;
            var maxScanNum = int.MinValue;

            foreach (var spec in massSpecDataReader.ReadAllSpectra())
            {
                _scanNumSpecMap.Add(spec.ScanNum, spec);
                _msLevel[spec.ScanNum] = spec.MsLevel;
                if (spec.MsLevel == 1)
                {
                    foreach (var peak in spec.Peaks)
                    {
                        _ms1PeakList.Add(new LcMsPeak(peak.Mz, peak.Intensity, spec.ScanNum));
                    }
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
        /// Gets the spectrum of the specified scan number
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>spectrum</returns>
        public Spectrum GetSpectrum(int scanNum)
        {
            return _scanNumSpecMap[scanNum];
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
                xic.Add(new XicPeak(peak.ScanNum, peak.Intensity));
                --i;
            }

            // go up
            i = index;
            while (i >= 0 && i < _ms1PeakList.Count)
            {
                var peak = _ms1PeakList[i];
                if (peak.Mz >= maxMz) break;
                xic.Add(new XicPeak(peak.ScanNum, peak.Intensity));
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
            //TODO: this tolerance value is not optimal for all data (revisit required)
            const int tolerance = 3;

            if (GetMsLevel(targetScanNum) > 1)
                targetScanNum = GetPrecursorScanNum(targetScanNum);
            var xic = GetExtractedIonChromatogram(minMz, maxMz);
            var index = xic.BinarySearch(new XicPeak(targetScanNum, 0));
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
            return _msLevel[scanNum];
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
            for (var curScanNum = scanNum + 1; curScanNum >= MaxLcScan; curScanNum++)
            {
                if (GetMsLevel(curScanNum) == msLevel) return curScanNum;
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
            using (var writer = new StreamWriter(outputFilePath))
            {
                for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++ )
                {
                    var spec = GetSpectrum(scanNum);
                    spec.WriteTo(writer);
                }
            }
        }

        private readonly Dictionary<int, int> _msLevel;
        private readonly int[] _precursorScan;

        private readonly Dictionary<int, Spectrum> _scanNumSpecMap;  // scan number -> spectrum
        private readonly List<LcMsPeak> _ms1PeakList;
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
