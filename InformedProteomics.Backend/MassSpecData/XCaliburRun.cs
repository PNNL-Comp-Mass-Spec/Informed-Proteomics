using System;
using System.Collections.Generic;
using System.Text.RegularExpressions;
using DeconTools.Backend.ProcessingTasks.PeakDetectors;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public class XCaliburRun
    {
        // Parameters for centroiding spectra
        public const int PeakToBackgroundRatio = 0;
        public const int SignalToNoiseThreshold = 0;

        public XCaliburRun(string filePath)
        {
            _msfileReader = (MSFileReaderLib.IXRawfile5) new MSFileReaderLib.MSFileReader_XRawfile();

            _msfileReader.Open(filePath);
            _msfileReader.SetCurrentController(0, 1);

            MinLcScan = 1;
            MaxLcScan = GetNumSpectra();

            // determine all MS levels and precursor scans
            _msLevel = new int[MaxLcScan-MinLcScan+1];
            _precursorScan = new int[MaxLcScan - MinLcScan + 1];
            var precursorMap = new Dictionary<int, int>();

            var prevMsLevel = int.MinValue;
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                int index = scanNum - MinLcScan;
                var msLevel = GetMsLevelFromRawData(scanNum);
                _msLevel[index] = msLevel;
                // determine precursor scan
                if (msLevel == prevMsLevel)
                {
                    _precursorScan[index] = _precursorScan[index - 1];
                }
                else if(msLevel > prevMsLevel)
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
        }

        public void Close()
        {
            if (_msfileReader != null)
            {
                _msfileReader.Close();
            }
        }

        public int MinLcScan { get; private set; }
        public int MaxLcScan { get; private set; }

        public Spectrum GetMassSpectrum(int scanNum)
        {
            object peakArr = null;
            object peakFlags = null;
            var arraySize = 0;

            _msfileReader.GetMassListFromScanNum(
                ref scanNum,
                null, // optional scan scanFilterString
                0, // intensity cutoff type
                0, // intensity cutoff value
                0, // max number of peaks
                0, // return centroided spectrum if nonzero
                0, // centroid peak width
                ref peakArr,
                ref peakFlags,
                ref arraySize);

            var vals = peakArr as double[,];

            if (vals == null) return null;

            var mzArr = new double[vals.GetLength(1)];
            var intensityArr = new double[vals.GetLength(1)];

            var sortRequired = false;

            for (var i = 0; i < vals.GetLength(1); i++)
            {
                mzArr[i] = vals[0, i];
                intensityArr[i] = vals[1, i];

                if (i > 0 && mzArr[i] < mzArr[i - 1])
                    sortRequired = true;
            }

            if (sortRequired)
            {
                Array.Sort(mzArr, intensityArr);
            }

            // Centroid spectrum
            if (!IsCentroidScan(scanNum))
            {
                if (_peakDetector == null)
                {
                    _peakDetector = new DeconToolsPeakDetectorV2(PeakToBackgroundRatio, SignalToNoiseThreshold);
                }
                var peaks = _peakDetector.FindPeaks(mzArr, intensityArr);
                mzArr = new double[peaks.Count];
                intensityArr = new double[peaks.Count];

                var index = -1;
                foreach (var peak in peaks)
                {
                    mzArr[++index] = peak.XValue;
                    intensityArr[index] = peak.Height;
                }
            }

            var msLevel = GetMsLevel(scanNum);
            if(msLevel == 1) return new Spectrum(mzArr, intensityArr);

            return new ProductSpectrum(mzArr, intensityArr)
                {
                    MsLevel = msLevel,
                    ActivationMethod = GetActivationMethod(scanNum),
                    PrecursorInfo = GetPrecursorInfo(scanNum)
                };
        }

        public PrecursorInfo GetPrecursorInfo(int scanNum)
        {
            // Get precursor scan number
            var precursorScanNum = GetPrecursorScanNum(scanNum);

            string scanFilterString = null;
            _msfileReader.GetFilterForScanNum(scanNum, ref scanFilterString);

            // Get Isolation Target m/z
            var isolationTargetMz = -1.0;
            if (scanFilterString != null)
            {
                isolationTargetMz = ParseMzValueFromThermoScanInfo(scanFilterString);
            }

            // Get isolation window width
            var isolationWindowWidth = GetIsolationWidth(scanNum);

            return new PrecursorInfo(precursorScanNum, isolationTargetMz, isolationWindowWidth / 2, isolationWindowWidth / 2);
        }

        public int GetPrecursorScanNum(int scanNum)
        {
            return _precursorScan[scanNum - MinLcScan];
        }

        public int GetMsLevel(int scanNum)
        {
            return _msLevel[scanNum - MinLcScan];
        }

        public int GetNumSpectra()
        {
            var numSpectra = 0;

            _msfileReader.GetNumSpectra(ref numSpectra);
            return numSpectra;
        }

        public double GetIsolationWidth(int scanNum)
        {
            object value = null;
            _msfileReader.GetTrailerExtraValueForScanNum(scanNum, "MS2 Isolation Width:", ref value);

            return Convert.ToDouble(value);
        }

        public ActivationMethod GetActivationMethod(int scanNum)
        {
            var activationType = 0;
            _msfileReader.GetActivationTypeForScanNum(scanNum, GetMsLevel(scanNum), ref activationType);
            switch (activationType)
            {
                case 0: 
                    return ActivationMethod.CID;
                case 2:
                    return ActivationMethod.ECD;
                case 3:
                    return ActivationMethod.PQD;
                case 4:
                    return ActivationMethod.ETD;
                case 5:
                    return ActivationMethod.HCD;
                default:
                    return ActivationMethod.Unknown;
            }
        }

        public bool IsCentroidScan(int scanNum)
        {
            var isCentroid = 0;
            _msfileReader.IsCentroidScanForScanNum(scanNum, ref isCentroid);
            return isCentroid != 0;
        }

        private readonly MSFileReaderLib.IXRawfile5 _msfileReader;
        private readonly int[] _msLevel;
        private readonly int[] _precursorScan;
        private static DeconToolsPeakDetectorV2 _peakDetector;  // basic peak detector

        private int GetMsLevelFromRawData(int scanNum)
        {
            var msLevel = 0;
            _msfileReader.GetMSOrderForScanNum(scanNum, ref msLevel);

            return msLevel;
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
    }
}
