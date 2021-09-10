using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using PSI_Interface.CV;

// ReSharper disable UnusedMember.Global

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Base class for objects that provide access to LCMS run data
    /// </summary>
    public abstract class LcMsRun : ILcMsRun, IMassSpecDataReader
    {
        // Ignore Spelling: xic

        /// <summary>
        /// Number of unique isolation windows kept for DIA data
        /// </summary>
        public const int NumUniqueIsolationWindowThresholdForDia = 1000;

        /// <summary>
        /// Factor used to bin isolation window data
        /// </summary>
        public const double IsolationWindowBinningFactor = 10;

        /// <summary>
        /// Index of first LC scan in the dataset
        /// </summary>
        public int MinLcScan { get; protected set; }

        /// <summary>
        /// Index of last LC scan in the dataset
        /// </summary>
        public int MaxLcScan { get; protected set; }

        /// <summary>
        /// The number of spectra in the file.
        /// </summary>
        public int NumSpectra { get; protected set; }

        /// <summary>
        /// Lowest MS Level in the dataset. Usually 1.
        /// </summary>
        public int MinMsLevel { get; protected set; }

        /// <summary>
        /// Highest MS Level in the dataset.
        /// </summary>
        public int MaxMsLevel { get; protected set; }

        /// <summary>
        /// List of all scan numbers in the dataset
        /// </summary>
        public IEnumerable<int> AllScanNumbers { get { return ScanNumToMsLevel.Select(x => x.Key); } }

        /// <summary>
        /// The smallest MS1 m/z
        /// </summary>
        public abstract double MinMs1Mz { get; }

        /// <summary>
        /// The largest MS1 m/z
        /// </summary>
        public abstract double MaxMs1Mz { get; }

        /// <summary>
        /// True if the dataset is DIA data
        /// </summary>
        public bool IsDia => IsDiaOrNull ?? (bool)(IsDiaOrNull = GetNumUniqueIsoWindows() < NumUniqueIsolationWindowThresholdForDia);

        /// <summary>
        /// Default values for configuration properties
        /// </summary>
        protected LcMsRun()
        {
            HigherPrecursorChromatogramCacheSize = 0;
            LowerPrecursorChromatogramCacheSize = 0;
        }

        #region IMassSpecDataReader

        /// <summary>
        /// Close the reader
        /// </summary>
        public abstract void Close();

        /// <summary>
        /// Properly dispose of all unmanaged resources (specifically, file handles)
        /// </summary>
        public abstract void Dispose();

        /// <summary>
        /// Gets all spectra
        /// </summary>
        /// <param name="includePeaks"></param>
        /// <returns>all spectra</returns>
        public IEnumerable<Spectrum> ReadAllSpectra(bool includePeaks = true)
        {
            return AllScanNumbers.OrderBy(x => x).Select(x => ReadMassSpectrum(x, includePeaks));
        }

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000767, native spectrum identifier format
        /// </summary>
        public CV.CVID NativeIdFormat { get; protected set; }

        /// <summary>
        /// The Native Format of the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </summary>
        public CV.CVID NativeFormat { get; protected set; }

        /// <summary>
        /// Path to the file; is <see cref="string.Empty"/> if the reader is in-memory
        /// </summary>
        public abstract string FilePath { get; protected set; }

        /// <summary>
        /// SHA-1 Checksum of the original input file (raw, mzML, .d folder, etc.)
        /// </summary>
        public abstract string SrcFileChecksum { get; protected set; }

        /// <summary>
        /// Version of the immediate prior input file (raw, mzML, .d folder, etc.)
        /// </summary>
        public abstract string FileFormatVersion { get; protected set; }

        /// <summary>
        /// Returns the spectrum specified by the scan number.
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="includePeaks"></param>
        /// <returns>Spectrum object</returns>
        public Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true)
        {
            return GetSpectrum(scanNum, includePeaks);
        }

        /// <summary>
        /// Try to make the reader random access capable
        /// </summary>
        /// <returns>true if is random access capable, false if not</returns>
        public abstract bool TryMakeRandomAccessCapable();

        #endregion

        #region Spectra and scan operations

        /// <summary>
        /// Return the specified spectrum from the file, optionally reading only the metadata
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="includePeaks"></param>
        /// <returns>Spectrum object</returns>
        public abstract Spectrum GetSpectrum(int scanNum, bool includePeaks = true);

        /// <summary>
        /// If <paramref name="scanNum"/> is a MS1 scan, return it; otherwise, return null.
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="ms1ScanIndex"></param>
        /// <returns>Spectrum object</returns>
        public abstract Spectrum GetMs1Spectrum(int scanNum, out int ms1ScanIndex);

        /// <summary>
        /// Return the isolation window for the specified scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns>Isolation window object</returns>
        public abstract IsolationWindow GetIsolationWindow(int scanNum);

        /*
        public ProductSpectrum GetSummedMs2Spectrum(double monoIsotopicMass, Ms1Feature range,
            ActivationMethod activationMethod = ActivationMethod.Unknown)
        {
            return GetSummedMs2Spectrum(monoIsotopicMass, range.MinScanNum, range.MaxScanNum, range.MinCharge,
                range.MaxCharge, activationMethod);
        }*/

        #endregion

        #region Detail data by scan number

        /// <summary>
        /// Gets the MS Level of the specified scan
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>MS Level</returns>
        public int GetMsLevel(int scanNum)
        {
            return ScanNumToMsLevel.TryGetValue(scanNum, out var msLevel) ? msLevel : 0;
        }

        /// <summary>
        /// Get the elution time of the specified scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns>Elution time</returns>
        public double GetElutionTime(int scanNum)
        {
            return ScanNumElutionTimeMap.TryGetValue(scanNum, out var elutionTime) ? elutionTime : 0.0;
        }

        /// <summary>
        /// Get the native ID of the specified scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns>Native ID string</returns>
        public virtual string GetNativeId(int scanNum)
        {
            return ScanNumNativeIdMap.TryGetValue(scanNum, out var nativeId) ? nativeId : string.Empty;
        }

        /// <summary>
        /// Get the ion mobility drift time of the specified scan number (if data is ion mobility)
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns>Drift time</returns>
        public double GetDriftTime(int scanNum)
        {
            // TODO: Not implemented for LcMsRun, not stored in PbfLcMsRun, so return 0.
            return 0;
        }

        #endregion

        #region Scan numbers

        /// <summary>
        /// Gets the precursor scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns>precursor scan number or 0 for MS1</returns>
        public int GetPrecursorScanNum(int scanNum)
        {
            if (_precursorScan.ContainsKey(scanNum))
            {
                return _precursorScan[scanNum];
            }
            return 0;
        }

        /// <summary>
        /// Gets the next scan number whose MS Level is smaller by 1
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>next scan number or MaxLc for MS1</returns>
        public int GetNextScanNum(int scanNum)
        {
            if (_nextScan.ContainsKey(scanNum))
            {
                return _nextScan[scanNum];
            }
            return MaxLcScan;
        }

        /// <summary>
        /// Gets the greatest scan number smaller than ms2ScanNum
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <param name="msLevel">MS Level</param>
        /// <returns>previous scan number at the specified level</returns>
        public int GetPrevScanNum(int scanNum, int msLevel)
        {
            return ScanNumToMsLevel.Where(x => x.Value == msLevel && x.Key < scanNum).DefaultIfEmpty(new KeyValuePair<int, int>(MinLcScan - 1, 0)).Max(x => x.Key);
        }

        /// <summary>
        /// Gets the smallest scan number larger than ms2ScanNum, for the given MS Level
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <param name="msLevel">MS Level</param>
        /// <returns>next scan number at the specified level</returns>
        public int GetNextScanNum(int scanNum, int msLevel)
        {
            return ScanNumToMsLevel.Where(x => x.Value == msLevel && x.Key > scanNum).DefaultIfEmpty(new KeyValuePair<int, int>(MaxLcScan + 1, 0)).Min(x => x.Key);
        }

        /// <summary>
        /// Gets the scan numbers of the specified msLevel
        /// </summary>
        /// <param name="msLevel">MS Level</param>
        /// <returns>scan numbers of the specified msLevel</returns>
        public IList<int> GetScanNumbers(int msLevel)
        {
            return ScanNumToMsLevel.Where(x => x.Value == msLevel).Select(x => x.Key).ToList();
        }

        private int[] _ms1ScanVector;

        /// <summary>
        /// An array of all of the MS1 scan numbers
        /// </summary>
        /// <returns>MS1 scan numbers</returns>
        public int[] GetMs1ScanVector()
        {
            return _ms1ScanVector ??= GetScanNumbers(1).ToArray();
        }

        private int[] _ms1ScanNumToIndex;

        /// <summary>
        /// Array of length MaxLcScan where entries that are non-zero are the scan index of the given scan number
        /// </summary>
        /// <returns>Scan number</returns>
        /// <remarks>
        /// For example, if scan 7 is the 5th MS1 scan in the dataset, then _ms1ScanNumToIndex[7] is 4
        /// Entries in the array that are 0 mean that MS1 scan does not map to an index
        /// (exception: scan 1 is listed as index 0)
        /// </remarks>
        public int[] GetMs1ScanNumToIndex()
        {
            if (_ms1ScanNumToIndex != null)
            {
                return _ms1ScanNumToIndex;
            }

            var ms1ScanNumbers = GetMs1ScanVector();
            _ms1ScanNumToIndex = new int[MaxLcScan + 1];
            for (var i = 0; i < ms1ScanNumbers.Length; i++)
            {
                _ms1ScanNumToIndex[ms1ScanNumbers[i]] = i;
            }

            return _ms1ScanNumToIndex;
        }

        #endregion

        #region Precursor XICs

        /// <summary>
        /// Number of extra precursor chromatogram points to cache on the higher-m/z side of a requested XIC
        /// </summary>
        /// <remarks>This will be ignored for any size less than 20 (i.e., no caching on the higher-m/z side will occur)
        /// This property is specifically designed for use in getting large numbers of XICs that are very close in m/z; a high value will cause degraded performance when reading XICs at random</remarks>
        public int HigherPrecursorChromatogramCacheSize { get; set; }

        /// <summary>
        /// Number of extra precursor chromatogram points to cache on the lower-m/z side of a requested XIC
        /// </summary>
        /// <remarks>This will be ignored for any size less than 20 (i.e., no caching on the lower-m/z side will occur)
        /// This property is specifically designed for use in getting large numbers of XICs that are very close in m/z; a high value will cause degraded performance when reading XICs at random</remarks>
        public int LowerPrecursorChromatogramCacheSize { get; set; }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z (using only MS1 spectra)
        /// XicPoint is created for every MS1 scan.
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <returns>XIC object</returns>
        public Xic GetFullPrecursorIonExtractedIonChromatogram(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsMz(mz);
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
        /// <returns>XIC object</returns>
        public Xic GetFullPrecursorIonExtractedIonChromatogram(double minMz, double maxMz)
        {
            var xic = GetPrecursorExtractedIonChromatogram(minMz, maxMz);

            var hasXicPoint = new bool[MaxLcScan - MinLcScan + 1];
            foreach (var xicPoint in xic)
            {
                hasXicPoint[xicPoint.ScanNum - MinLcScan] = true;
            }

            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                if (GetMsLevel(scanNum) == 1 && !hasXicPoint[scanNum - MinLcScan])
                {
                    xic.Add(new XicPoint(scanNum, 0, 0));
                }
            }
            xic.Sort();
            return xic;
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS1 spectra)
        /// </summary>
        /// <param name="minMz">min m/z</param>
        /// <param name="maxMz">max m/z</param>
        /// <returns>An array of doubles, with every intensity value in the provided m/z range</returns>
        public double[] GetFullPrecursorIonExtractedIonChromatogramVector(double minMz, double maxMz)
        {
            var xic = GetPrecursorExtractedIonChromatogram(minMz, maxMz);

            var xicIntensityVector = new double[MaxLcScan - MinLcScan + 1];
            foreach (var xicPoint in xic)
            {
                xicIntensityVector[xicPoint.ScanNum - MinLcScan] = xicPoint.Intensity;
            }

            var intVector = GetMs1ScanVector().Select(s => xicIntensityVector[s - MinLcScan]).ToArray();
            return intVector;
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z (using only MS1 spectra)
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <returns>XIC object</returns>
        public Xic GetPrecursorExtractedIonChromatogram(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsMz(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetPrecursorExtractedIonChromatogram(minMz, maxMz);
        }

        /// <summary>
        /// Returns selected peaks between minMz and maxMz. The biggest peak per scan is selected.
        /// </summary>
        /// <param name="minMz"></param>
        /// <param name="maxMz"></param>
        /// <returns>Extracted ion chromatogram data</returns>
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
            var tolTh = tolerance.GetToleranceAsMz(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            if (targetScanNum < 0)
            {
                return GetPrecursorExtractedIonChromatogram(minMz, maxMz);
            }

            return GetPrecursorExtractedIonChromatogram(minMz, maxMz, targetScanNum, maxNumConsecutiveScansWithoutPeak);
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
            {
                targetScanNum = GetPrecursorScanNum(targetScanNum);
            }

            var xic = GetPrecursorExtractedIonChromatogram(minMz, maxMz);
            return GetTrimmedXic(xic, targetScanNum, tolerance);
        }

        /// <summary>
        /// Returns all precursor peaks between minMz and maxMz, including multiple peaks per scan
        /// </summary>
        /// <param name="minMz"></param>
        /// <param name="maxMz"></param>
        /// <returns>Extracted ion chromatogram data</returns>
        public abstract Xic GetPrecursorChromatogramRange(double minMz, double maxMz);

        #endregion

        #region Product XICs

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS2 spectra)
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <param name="precursorIonMz">precursor m/z of the precursor ion</param>
        /// <returns>XIC object</returns>
        public Xic GetFullProductExtractedIonChromatogram(double mz, Tolerance tolerance, double precursorIonMz)
        {
            var tolTh = tolerance.GetToleranceAsMz(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetFullProductExtractedIonChromatogram(minMz, maxMz, precursorIonMz);
        }

        /// <summary>
        /// Returns a xic for the chosen range that covers the entire run.
        /// </summary>
        /// <param name="minMz"></param>
        /// <param name="maxMz"></param>
        /// <param name="precursorMz"></param>
        /// <returns>Extracted ion chromatogram data</returns>
        public abstract Xic GetFullProductExtractedIonChromatogram(double minMz, double maxMz, double precursorMz);

        #endregion

        #region XIC general

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
            if (index < 0)
            {
                index = ~index;
            }

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
                    if (GetMsLevel(scanNum) == 1)
                    {
                        ++numMissingScans;
                    }

                    if (numMissingScans > tolerance)
                    {
                        isConsecutive = false;
                        break;
                    }
                }
                if (isConsecutive)
                {
                    xicSegment.Add(xicPeak);
                }
                else
                {
                    break;
                }

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
                    if (GetMsLevel(scanNum) == 1)
                    {
                        ++numMissingScans;
                    }

                    if (numMissingScans > tolerance)
                    {
                        isConsecutive = false;
                        break;
                    }
                }
                if (isConsecutive)
                {
                    xicSegment.Add(xicPeak);
                }
                else
                {
                    break;
                }

                curScanNum = xicPeak.ScanNum;
                ++i;
            }

            xicSegment.Sort();
            return xicSegment;
        }

        #endregion

        #region Isolation Windows

        /// <summary>
        /// Return the number of unique isolation windows in the dataset
        /// </summary>
        /// <returns>Isolation window count</returns>
        public int GetNumUniqueIsoWindows()
        {
            var isoWindowSet = new HashSet<IsolationWindow>();
            foreach (var scanNum in ScanNumToMsLevel.Where(x => x.Value > 1).Select(x => x.Key))
            {
                if (GetSpectrum(scanNum) is not ProductSpectrum productSpec)
                {
                    continue;
                }

                isoWindowSet.Add(productSpec.IsolationWindow);
            }
            return isoWindowSet.Count;
        }

        /// <summary>
        /// Get the narrowest isolation window width
        /// </summary>
        /// <returns>Narrowest window width</returns>
        public double GetMinIsolationWindowWidth()
        {
            var minWidth = double.MaxValue;
            foreach (var scanNum in ScanNumToMsLevel.Where(x => x.Value > 1).Select(x => x.Key))
            {
                if (GetSpectrum(scanNum) is not ProductSpectrum productSpec)
                {
                    continue;
                }

                if (productSpec.IsolationWindow.Width < minWidth)
                {
                    minWidth = productSpec.IsolationWindow.Width;
                }
            }
            return minWidth;
        }

        #endregion

        #region Fragmentation spectra scan numbers

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
            return IsolationMzBinToScanNums.TryGetValue(targetIsoBin, out var scanNums) ? scanNums : Array.Empty<int>();
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

        #endregion

        #region Protected members and functions

        // Fields to be defined in a child
        /// <summary>
        /// Dictionary to map IsolationMzBins to scan numbers
        /// </summary>
        protected Dictionary<int, int[]> IsolationMzBinToScanNums;

        /// <summary>
        /// Dictionary to map scan numbers to MS Levels
        /// </summary>
        protected Dictionary<int, int> ScanNumToMsLevel;

        /// <summary>
        /// Dictionary to map scan numbers to elution times
        /// </summary>
        protected Dictionary<int, double> ScanNumElutionTimeMap;

        /// <summary>
        /// Dictionary to map scan numbers to native IDs
        /// </summary>
        protected Dictionary<int, string> ScanNumNativeIdMap = new();

        /// <summary>
        /// True if DIA data, false if not, null if unknown
        /// </summary>
        protected bool? IsDiaOrNull;

        /// <summary>
        /// This dictionary tracks, for each scan, its corresponding precursor scan number (if an MS2 scan); for MS1 scans, it stores 0
        /// </summary>
        private Dictionary<int, int> _precursorScan;

        /// <summary>
        /// This dictionary tracks, for each scan, the next scan number whose MS Level is smaller by 1
        /// </summary>
        /// <remarks>
        /// For all MS1 scans, it stores one more than the final scan number in the file
        /// </remarks>
        private Dictionary<int, int> _nextScan;

        /// <summary>
        /// Create the maps for linking MSn scans to their precursors, and for getting the next MS1 scan number given a scan number
        /// </summary>
        protected void CreatePrecursorNextScanMap()
        {
            _precursorScan = new Dictionary<int, int>(NumSpectra + 1);
            _nextScan = new Dictionary<int, int>(NumSpectra + 1);

            // This dictionary tracks the most recent precursor scan number for the given MS Level
            var precursorMap = new Dictionary<int, int>();

            // This dictionary tracks the next scan number for the given MS Level
            var nextScanMap = new Dictionary<int, int>();

            // Initialize the dictionaries to defaults
            for (var msLevel = MinMsLevel; msLevel <= MaxMsLevel; msLevel++)
            {
                precursorMap[msLevel] = 0;
                nextScanMap[msLevel] = MaxLcScan + 1;
            }

            var prevMsLevel = int.MinValue;
            var prevScanNum = -1;
            foreach (var scan in ScanNumToMsLevel.OrderBy(x => x.Key))
            {
                var scanNum = scan.Key;
                var msLevel = scan.Value;

                if (msLevel == 0)
                {
                    continue; // corrupted scan
                }

                // determine precursor scan
                if (msLevel == prevMsLevel)
                {
                    if (prevScanNum < 0)
                    {
                        // We shouldn't encounter this, but better safe than sorry.
                        _precursorScan[scanNum] = MinLcScan;
                    }
                    else
                    {
                        _precursorScan[scanNum] = _precursorScan[prevScanNum];
                    }
                }
                else if (msLevel > prevMsLevel)
                {
                    _precursorScan[scanNum] = scanNum - 1;
                    precursorMap[msLevel] = scanNum - 1;
                }
                else // msLevel < prevMsLevel
                {
                    _precursorScan[scanNum] = precursorMap[msLevel];
                }
                prevMsLevel = msLevel;
                prevScanNum = scanNum;
            }

            var nextMsLevel = int.MinValue;
            var nextScanNum = int.MaxValue;
            foreach (var scan in ScanNumToMsLevel.OrderByDescending(x => x.Key))
            {
                var scanNum = scan.Key;
                var msLevel = scan.Value;

                if (msLevel == 0)
                {
                    continue; // corrupted scan
                }

                // determine next scan
                if (msLevel == nextMsLevel)
                {
                    if (nextScanNum == int.MaxValue)
                    {
                        // We shouldn't encounter this, but better safe than sorry.
                        _nextScan[scanNum] = MaxLcScan;
                    }
                    else
                    {
                        _nextScan[scanNum] = _nextScan[nextScanNum];
                    }
                }
                else if (msLevel > nextMsLevel)
                {
                    _nextScan[scanNum] = scanNum + 1;
                    nextScanMap[msLevel] = scanNum + 1;
                }
                else // msLevel < prevMsLevel
                {
                    _nextScan[scanNum] = nextScanMap[msLevel];
                }
                nextMsLevel = msLevel;
                nextScanNum = scanNum;
            }
        }

        #endregion
    }
}
