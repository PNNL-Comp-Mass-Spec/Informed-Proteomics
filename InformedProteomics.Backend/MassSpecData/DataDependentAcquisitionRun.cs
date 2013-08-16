using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public class DataDependentAcquisitionRun
    {
        public DataDependentAcquisitionRun(IRawDataReader rawDataReader)
        {
            _rawDataReader = rawDataReader;

            //_ms1List = new LinkedList<Spectrum>();
            _scanNumSpecMap = new Dictionary<int, Spectrum>();
            _ms1PeakList = new List<LcMsPeak>();

            // Parse all spectra
            var numMs1Scans = 0;
            foreach (var spec in _rawDataReader.ReadAllSpectra())
            {
                _scanNumSpecMap.Add(spec.ScanNum, spec);
                if (spec.MsLevel == 1)
                {
                    //_ms1List.AddLast(spec); // MS1 spectrum                    
                    ++numMs1Scans;
                    foreach (var peak in spec.Peaks)
                    {
                        _ms1PeakList.Add(new LcMsPeak(peak.Mz, peak.Intensity, spec.ScanNum));
                    }
                }
            }

            _numMs1Scans = numMs1Scans;
            _ms1PeakList.Sort();
        }

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
        /// <returns>XIC as a list of XICPeaks</returns>
        public IList<XicPeak> GetExtractedIonChromatogram(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetExtractedIonChromatogram(minMz, maxMz);
        }

        //public IList<XicPeak> GetExtractedIonChromatogramOld(double mz, Tolerance tolerance)
        //{
        //    var xic = new List<XicPeak>();
        //    foreach (var spec in _ms1List)
        //    {
        //        var peak = spec.FindPeak(mz, tolerance);
        //        if (peak != null)
        //        {
        //            xic.Add(new XicPeak(spec.ScanNum, peak.Intensity));
        //        }
        //    }
        //    return xic;
        //}


        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z (using only MS1 spectra)
        /// </summary>
        /// <param name="minMz">min m/z</param>
        /// <param name="maxMz">max m/z</param>
        /// <returns>XIC as a list of XICPeaks</returns>
        public IList<XicPeak> GetExtractedIonChromatogram(double minMz, double maxMz)
        {
            var xicPeaks = new List<XicPeak>();

            var index = _ms1PeakList.BinarySearch(new LcMsPeak((minMz + maxMz) / 2, 0, 0));
            if (index < 0) index = ~index;

            // go down
            var i = index - 1;
            while (i >= 0 && i < _ms1PeakList.Count)
            {
                var peak = _ms1PeakList[i];
                if (peak.Mz <= minMz) break;
                xicPeaks.Add(new XicPeak(peak.ScanNum, peak.Intensity));
                --i;
            }

            // go up
            i = index;
            while (i >= 0 && i < _ms1PeakList.Count)
            {
                var peak = _ms1PeakList[i];
                if (peak.Mz >= maxMz) break;
                xicPeaks.Add(new XicPeak(peak.ScanNum, peak.Intensity));
                ++i;
            }

            if (xicPeaks.Count == 0) return xicPeaks;

            xicPeaks.Sort();
            var xic = new List<XicPeak>();
            var prevScanNum = -1;
            var bestPeak = xicPeaks[0];
            for (i=1; i<xicPeaks.Count; i++)
            {
                var xicPeak = xicPeaks[i];
                if (xicPeak.ScanNum > prevScanNum)
                {
                    xic.Add(bestPeak);
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
            xic.Add(bestPeak);
            return xic;
        }

        /// <summary>
        /// Gets the scan numbers of the specified msLevel
        /// </summary>
        /// <param name="msLevel">MS level</param>
        /// <returns>scan numbers of the specified msLevel</returns>
        public IEnumerable<int> GetScanNumbers(int msLevel)
        {
            return _rawDataReader.GetScanNumbers(msLevel);
        }

        /// <summary>
        /// Gets MS/MS spectra whose isolation windows contain the most abundant peak of the precursorIon
        /// </summary>
        /// <param name="precursorIon"></param>
        /// <returns></returns>
        public IEnumerable<ProductSpectrum> GetMatchingMs2Spectra(Ion precursorIon)
        {
            throw new System.NotImplementedException();
        }

        private readonly IRawDataReader _rawDataReader;
        private readonly Dictionary<int, Spectrum> _scanNumSpecMap;  // scan number -> spectrum
        //private readonly LinkedList<Spectrum> _ms1List; // linked list of MS1 spectra
        private readonly List<LcMsPeak> _ms1PeakList;
        private readonly int _numMs1Scans;
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
