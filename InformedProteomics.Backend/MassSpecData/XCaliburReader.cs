using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Spectrometry;
using ThermoRawFileReader;

namespace InformedProteomics.Backend.MassSpecData
{
    public class XCaliburReader : IMassSpecDataReader
    {
        // Parameters for centroiding spectra
        public const int PeakToBackgroundRatio = 0;

        public XCaliburReader(string filePath)
        {
            _cachedScanInfo = new clsScanInfo(-1);

            _msfileReader = new XRawFileIO();

            var dataFile = new FileInfo(filePath);
            if (!dataFile.Exists)
                throw new FileNotFoundException("Thermo .raw file not found: " + filePath, dataFile.FullName);

            _msfileReader.OpenRawFile(filePath);

            _minLcScan = 1;
            _numSpectra = ReadNumSpectra();
            _maxLcScan = _numSpectra;

            _msLevel = new Dictionary<int, int>();
        }

        /// <summary>
        /// Reads all spectra
        /// </summary>
        /// <returns>all spectra</returns>
        [System.Runtime.ExceptionServices.HandleProcessCorruptedStateExceptions]
        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            for (var scanNum = _minLcScan; scanNum <= _maxLcScan; scanNum++)
            {
                Spectrum spec = null;
                try
                {
                    spec = ReadMassSpectrum(scanNum);
                }
                catch (System.Runtime.InteropServices.COMException/* ex*/)
                {
                    Console.WriteLine("[Warning] Ignore corrupted spectrum Scan={0}", scanNum);
                }

                if (spec != null) yield return spec;
            }
        }

        /// <summary>
        /// Always random-access capable.
        /// </summary>
        /// <returns></returns>
        public bool TryMakeRandomAccessCapable()
        {
            return true;
        }

        /// <summary>
        /// Reads the mass spectrum with the specified scanNum from the raw file
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <param name="includePeaks">whether to include peak data</param>
        /// <returns>mass spectrum</returns>
        public Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true)
        {
            var scanInfo = GetScanInfo(scanNum);

            // default empty arrays, if peak data not requested.
            double[] mzArr = new double[]{};
            double[] intensityArr = new double[]{};

            if (includePeaks)
            {
                _msfileReader.GetScanData(scanNum, out mzArr, out intensityArr, 0, true);
            }

            var elutionTime = RtFromScanNum(scanNum);
            var nativeId = "controllerType=0 controllerNumber=1 scan=" + scanNum;

            // Call scanInfo.MSLevel in order to update dictionary _msLevel
            var msLevel = ReadMsLevel(scanNum);

            if (msLevel == 1)
                return new Spectrum(mzArr, intensityArr, scanNum)
                {
                    ElutionTime = elutionTime,
                    TotalIonCurrent = scanInfo.TotalIonCurrent,
                    NativeId = nativeId,
                };

            var isolationWindow = ReadPrecursorInfo(scanNum);

            var productSpec = new ProductSpectrum(mzArr, intensityArr, scanNum)
            {
                MsLevel = scanInfo.MSLevel,
                ElutionTime = elutionTime,
                TotalIonCurrent = scanInfo.TotalIonCurrent,
                NativeId = nativeId,
                ActivationMethod = GetActivationMethod(scanNum),
                IsolationWindow = isolationWindow
            };
            return productSpec;
        }

        /// <summary>
        /// Reads the precursor information of the specified scan
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>precursor information</returns>
        public IsolationWindow ReadPrecursorInfo(int scanNum)
        {
            if (ReadMsLevel(scanNum) <= 1) return null;

            // Get Isolation Target m/z
            var isolationTargetMz = ReadIsolationWindowTargetMz(scanNum);

            // Get isolation window width
            var precursorInfo = ReadPrecursorInfoFromTrailerExtra(scanNum);

            var isolationWindowWidth = precursorInfo.IsolationWidth;
            var monoisotopicMz = precursorInfo.MonoisotopicMz;
            var charge = precursorInfo.Charge;

            return new IsolationWindow(
                isolationTargetMz,
                isolationWindowWidth / 2,
                isolationWindowWidth / 2,
                monoisotopicMz,
                charge
                );
        }

        public int GetMaxScanNum()
        {
            return _maxLcScan;
        }

        public int NumSpectra
        {
            get { return _numSpectra; }
        }

        public int GetMinScanNum()
        {
            return _minLcScan;
        }

        public int ReadMsLevel(int scanNum)
        {
            int msLevel;
            if (_msLevel.TryGetValue(scanNum, out msLevel))
                return msLevel;

            var scanInfo = GetScanInfo(scanNum);
            _msLevel[scanNum] = scanInfo.MSLevel;

            return scanInfo.MSLevel;
        }

        public double RtFromScanNum(int scanNum)
        {
            var scanInfo = GetScanInfo(scanNum);
            return scanInfo.RetentionTime;
        }

        // Obsolete: Now using ThermoRawFileReader")>
        //private readonly MSFileReaderLib.IXRawfile5 _msfileReader;
        private readonly XRawFileIO _msfileReader;

        private clsScanInfo _cachedScanInfo;

        private readonly int _minLcScan;
        private readonly int _maxLcScan;
        private readonly int _numSpectra;
        private readonly Dictionary<int, int> _msLevel;

        /// <summary>
        /// Reads the isolation window target m/z
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>isolation window target m/z</returns>
        private double ReadIsolationWindowTargetMz(int scanNum)
        {
            // string scanFilterString = null;
            // _msfileReader.GetFilterForScanNum(scanNum, ref scanFilterString);

            var scanInfo = GetScanInfo(scanNum);

            var isolationTargetMz = -1.0;
            if (!string.IsNullOrWhiteSpace(scanInfo.FilterText))
            {
                isolationTargetMz = ParseMzValueFromThermoScanInfo(scanInfo.FilterText);
            }
            return isolationTargetMz;
        }

        private PrecursorInfo ReadPrecursorInfoFromTrailerExtra(int scanNum)
        {
            if (ReadMsLevel(scanNum) == 1) return null;

            var isolationWidth = 0.0;
            double? monoIsotopicMz = 0.0;
            int? charge = null;

            var scanInfo = GetScanInfo(scanNum);
            string valueText;

            if (scanInfo.TryGetScanEvent("Monoisotopic M/Z:", out valueText))
            {
                monoIsotopicMz = Convert.ToDouble(valueText);
                if (monoIsotopicMz == 0.0)
                    monoIsotopicMz = null;
            }

            if (scanInfo.TryGetScanEvent("Charge State:", out valueText))
            {
                charge = Convert.ToInt32(valueText);
                if (charge == 0)
                    charge = null;
            }

            if (scanInfo.TryGetScanEvent("MS2 Isolation Width:", out valueText))
            {
                isolationWidth = Convert.ToDouble(valueText);
            }

            return new PrecursorInfo(isolationWidth, monoIsotopicMz, charge);
        }

        internal class PrecursorInfo
        {
            public PrecursorInfo(double isolationWidth, double? monoisotopicMz, int? charge)
            {
                MonoisotopicMz = monoisotopicMz;
                IsolationWidth = isolationWidth;
                Charge = charge;
            }

            internal double IsolationWidth { get; private set; }
            internal double? MonoisotopicMz { get; private set; }
            internal int? Charge { get; private set; }
        }

        private static double ParseMzValueFromThermoScanInfo(string scanFilterString)
        {
            const string pattern = @"(?<mz>[0-9.]+)@";

            var matchCid = Regex.Match(scanFilterString, pattern);

            if (matchCid.Success)
            {
                return Convert.ToDouble(matchCid.Groups["mz"].Value);
            }
            return -1;
        }

        private int ReadNumSpectra()
        {
            var numSpectra = _msfileReader.GetNumScans();

            return numSpectra;
        }

        /// <summary>
        /// Gets the activation method
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>activation method</returns>
        private ActivationMethod GetActivationMethod(int scanNum)
        {
            var scanInfo = GetScanInfo(scanNum);

            switch (scanInfo.ActivationType)
            {
                case ActivationTypeConstants.CID:
                    return ActivationMethod.CID;

                case ActivationTypeConstants.ECD:
                    return ActivationMethod.ECD;

                case ActivationTypeConstants.PQD:
                    return ActivationMethod.PQD;

                case ActivationTypeConstants.ETD:
                    return ActivationMethod.ETD;

                case ActivationTypeConstants.HCD:
                    return ActivationMethod.HCD;

                default:
                    return ActivationMethod.Unknown;
            }
        }

        private clsScanInfo GetScanInfo(int scanNum)
        {
            if (_cachedScanInfo == null || _cachedScanInfo.ScanNumber != scanNum)
                _msfileReader.GetScanInfo(scanNum, out _cachedScanInfo);

            return _cachedScanInfo;
        }

        public void Close()
        {
            if (_msfileReader != null)
                _msfileReader.CloseRawFile();
        }
    }
}
