using System;
using System.Collections.Generic;
using System.IO;
using System.Security.Cryptography;
using InformedProteomics.Backend.Data.Spectrometry;
using PRISM;
using PSI_Interface.CV;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.FilterEnums;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Class for reading spectra from Thermo .RAW files, using an installed MSFileReader DLL
    /// </summary>
    public class XcaliburReader : IMassSpecDataReader
    {
        /// <summary>
        /// Parameters for centroiding spectra
        /// </summary>
        public const int PeakToBackgroundRatio = 0;

        /// <summary>
        /// Constructor - open the file, and prepare to read
        /// </summary>
        /// <param name="filePath"></param>
        public XcaliburReader(string filePath)
        {
            var dataFile = new FileInfo(filePath);
            if (!dataFile.Exists)
            {
                throw new FileNotFoundException("Thermo .raw file not found: " + filePath, dataFile.FullName);
            }

            FilePath = filePath;

            _checkSum = string.Empty;

            _rawFileReader = RawFileReaderAdapter.FileFactory(filePath);
            _rawFileReader.SelectInstrument(Device.MS, 1);

            _minLcScan = 1;
            NumSpectra = ReadNumSpectra();
            _maxLcScan = NumSpectra;

            _msLevel = new Dictionary<int, int>();

            FileFormatVersion = _rawFileReader.FileHeader.Revision.ToString();
        }

        /// <summary>
        /// Reads all spectra
        /// </summary>
        /// <param name="includePeaks"></param>
        /// <returns>all spectra</returns>
        [System.Runtime.ExceptionServices.HandleProcessCorruptedStateExceptions]
        public IEnumerable<Spectrum> ReadAllSpectra(bool includePeaks = true)
        {
            for (var scanNum = _minLcScan; scanNum <= _maxLcScan; scanNum++)
            {
                Spectrum spec = null;
                try
                {
                    spec = ReadMassSpectrum(scanNum, includePeaks);
                }
                catch (System.Runtime.InteropServices.COMException/* ex*/)
                {
                    ConsoleMsgUtils.ShowWarning(string.Format("[Warning] Ignore corrupted spectrum Scan={0}", scanNum));
                }

                if (spec != null)
                {
                    yield return spec;
                }
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
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000767, native spectrum identifier format
        /// </summary>
        public CV.CVID NativeIdFormat => CV.CVID.MS_Thermo_nativeID_format;

        /// <summary>
        /// The Native Format of the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </summary>
        public CV.CVID NativeFormat => CV.CVID.MS_Thermo_RAW_format;

        /// <summary>
        /// Path to the file
        /// </summary>
        public string FilePath { get; }

        /// <summary>
        /// SHA-1 Checksum of the raw file
        /// </summary>
        /// <remarks>
        /// It can take some time to compute this value for .raw files over 500 MB,
        /// particularly if they're located on a remote share
        /// </remarks>
        public string SrcFileChecksum
        {
            get
            {
                if (string.IsNullOrEmpty(_checkSum))
                {
                    ComputeChecksum();
                }

                return _checkSum;
            }
        }

        /// <summary>
        /// Version of the file format
        /// </summary>
        public string FileFormatVersion { get; }

        /// <summary>
        /// Reads the mass spectrum with the specified scanNum from the raw file
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <param name="includePeaks">whether to include peak data</param>
        /// <returns>mass spectrum</returns>
        public Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true)
        {
            var scanFilter = GetScanFilter(scanNum);
            var scanStats = _rawFileReader.GetScanStatsForScanNumber(scanNum);

            // default empty arrays, if peak data not requested.
            var peaks = new List<Data.Spectrometry.Peak>();

            if (includePeaks)
            {
                var data = _rawFileReader.GetCentroidStream(scanNum, false);
                if (data.Masses?.Length > 0)
                {
                    // Always get the high-res centroid information if we can
                    peaks.Capacity = data.Masses.Length + 2;
                    for (var i = 0; i < data.Masses.Length; i++)
                    {
                        peaks.Add(new Data.Spectrometry.PeakDetailed(data.Masses[i], data.Intensities[i], data.Noises[i]));
                    }
                }
                else if (scanStats.IsCentroidScan)
                {
                    // Centroided, but doesn't have the centroid stream
                    var data2 = _rawFileReader.GetSegmentedScanFromScanNumber(scanNum, null);
                    peaks.Capacity = data2.Positions.Length + 2;
                    for (var i = 0; i < data2.Positions.Length; i++)
                    {
                        peaks.Add(new Data.Spectrometry.Peak(data2.Positions[i], data2.Intensities[i]));
                    }
                }
                else
                {
                    // Not centroided, no centroid data available - perform the centroiding operation from the vendor
                    var scanProf = Scan.FromFile(_rawFileReader, scanNum);
                    var centroided = Scan.ToCentroid(scanProf);

                    peaks.Capacity = centroided.PreferredMasses.Length + 2;
                    for (var i = 0; i < centroided.PreferredMasses.Length; i++)
                    {
                        peaks.Add(new Data.Spectrometry.Peak(centroided.PreferredMasses[i], centroided.PreferredIntensities[i]));
                    }
                }
            }

            var elutionTime = RtFromScanNum(scanNum);
            var nativeId = "controllerType=0 controllerNumber=1 scan=" + scanNum;

            // Call scanInfo.MSLevel in order to update dictionary _msLevel
            var msLevel = ReadMsLevel(scanNum);

            if (msLevel == 1)
            {
                return new Spectrum(peaks, scanNum)
                {
                    ElutionTime = elutionTime,
                    TotalIonCurrent = scanStats.TIC,
                    NativeId = nativeId,
                };
            }

            var isolationWindow = ReadPrecursorInfo(scanNum);

            var productSpec = new ProductSpectrum(peaks, scanNum)
            {
                MsLevel = (int)scanFilter.MSOrder, // MSOrder generally matches MSLevel
                ElutionTime = elutionTime,
                TotalIonCurrent = scanStats.TIC,
                NativeId = nativeId,
                ActivationMethod = GetActivationMethod(scanNum),
                IsolationWindow = isolationWindow
            };
            return productSpec;
        }

        /// <summary>
        /// Read the specified spectrum from the file, optionally reading only the metadata
        /// </summary>
        /// <param name="scanNum"></param>
        /// <param name="includePeaks"></param>
        /// <returns></returns>
        public Spectrum GetSpectrum(int scanNum, bool includePeaks = true)
        {
            return ReadMassSpectrum(scanNum, includePeaks);
        }

        /// <summary>
        /// Reads the precursor information of the specified scan
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>precursor information</returns>
        public IsolationWindow ReadPrecursorInfo(int scanNum)
        {
            if (ReadMsLevel(scanNum) <= 1)
            {
                return null;
            }

            // Get precursor info
            var preInfo = ReadPrecursorData(scanNum);

            return new IsolationWindow(
                preInfo.IsolationWindowTargetMz,
                preInfo.IsolationWidth / 2,
                preInfo.IsolationWidth / 2,
                preInfo.MonoisotopicMz,
                preInfo.Charge
                );
        }

        /// <summary>
        /// Get the maximum scan number in the file
        /// </summary>
        /// <returns></returns>
        public int GetMaxScanNum()
        {
            return _maxLcScan;
        }

        /// <summary>
        /// The number of spectra in the file.
        /// </summary>
        public int NumSpectra { get; }

        /// <summary>
        /// Get the minimum scan number in the file
        /// </summary>
        /// <returns></returns>
        public int GetMinScanNum()
        {
            return _minLcScan;
        }

        /// <summary>
        /// Read the MS Level of the specified scan number from the file
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns></returns>
        public int ReadMsLevel(int scanNum)
        {
            if (_msLevel.TryGetValue(scanNum, out var msLevel))
            {
                return msLevel;
            }

            var scanFilter = GetScanFilter(scanNum);
            _msLevel[scanNum] = (int)scanFilter.MSOrder; // MSOrder generally matches MSLevel

            return (int)scanFilter.MSOrder;
        }

        /// <summary>
        /// Get the retention time for the scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns></returns>
        public double RtFromScanNum(int scanNum)
        {
            return _rawFileReader.RetentionTimeFromScanNumber(scanNum);
        }

        private readonly IRawDataPlus _rawFileReader;

        private IScanFilter _cachedScanFilter = null;

        private int _cachedScanFilterScanNumber = -1;

        private readonly int _minLcScan;
        private readonly int _maxLcScan;
        private readonly Dictionary<int, int> _msLevel;
        private string _checkSum;

        /// <summary>
        /// Compute the SHA-1 checksum for this data file
        /// </summary>
        private void ComputeChecksum()
        {
            using (var fs = new FileStream(FilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
            using (var sha1 = new SHA1Managed())
            {
                var hash = sha1.ComputeHash(fs);
                _checkSum = BitConverter.ToString(hash).ToLower().Replace("-", "");
            }
        }

        /// <summary>
        /// Read the isolation window and precursor information
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns></returns>
        private PrecursorInfo ReadPrecursorData(int scanNum)
        {
            if (ReadMsLevel(scanNum) == 1)
            {
                return null;
            }

            var scanFilter = GetScanFilter(scanNum);
            var reactions = scanFilter.MassCount;

            var index = reactions - 1;
            if (index > 0 && scanFilter.GetIsMultipleActivation(index))
            {
                // The last activation is part of a ETciD/EThcD pair
                index--;
            }

            var preInfo = new PrecursorInfo();

            var reaction = scanFilter.GetReaction(index);
            preInfo.IsolationWindowTargetMz = reaction.PrecursorMass;
            preInfo.IsolationWindowWidth = reaction.IsolationWidth;
            preInfo.IsolationWindowOffset = reaction.IsolationWidthOffset;

            var extras = _rawFileReader.GetTrailerExtraInformation(scanNum);
            var monoMzIndex = -1;
            var chargeIndex = -1;
            var ms2IsolationWidthIndex = -1;
            for (var i = 0; i < extras.Length; i++)
            {
                var label = extras.Labels[i];
                if (label.Equals("Monoisotopic M/Z:", StringComparison.OrdinalIgnoreCase))
                {
                    monoMzIndex = i;
                }
                else if (label.Equals("Charge State:", StringComparison.OrdinalIgnoreCase))
                {
                    chargeIndex = i;
                }
                else if (label.Equals("MS2 Isolation Width:", StringComparison.OrdinalIgnoreCase))
                {
                    ms2IsolationWidthIndex = i;
                }
            }

            if (monoMzIndex >= 0)
            {
                preInfo.MonoisotopicMz = Convert.ToDouble(extras.Values[monoMzIndex]);
                if (Math.Abs(preInfo.MonoisotopicMz.Value) < float.Epsilon)
                {
                    preInfo.MonoisotopicMz = null;
                }
            }

            if (chargeIndex >= 0)
            {
                preInfo.Charge = Convert.ToInt32(extras.Values[chargeIndex]);
                if (preInfo.Charge == 0)
                {
                    preInfo.Charge = null;
                }
            }

            if (preInfo.MonoisotopicMz < float.Epsilon && preInfo.Charge == 0)
            {
                // Thermo software probably couldn't determine a charge state, and set the monoisotopic m/z to 0
                preInfo.MonoisotopicMz = reaction.PrecursorMass;
            }

            if (ms2IsolationWidthIndex >= 0)
            {
                preInfo.IsolationWidth = Convert.ToDouble(extras.Values[ms2IsolationWidthIndex]);
            }

            return preInfo;
        }

        private class PrecursorInfo
        {
            public PrecursorInfo()
            {
                IsolationWidth = -1;
                MonoisotopicMz = null;
                Charge = null;
                IsolationWindowTargetMz = -1;
                IsolationWindowWidth = -1;
            }

            public double IsolationWidth { get; set; }
            public double? MonoisotopicMz { get; set; }
            public int? Charge { get; set; }
            public double IsolationWindowTargetMz { get; set; }
            public double IsolationWindowWidth { get; set; }
            public double IsolationWindowOffset { get; set; }
        }

        private int ReadNumSpectra()
        {
            var numSpectra = _rawFileReader.RunHeaderEx.SpectraCount;

            return numSpectra;
        }

        /// <summary>
        /// Gets the activation method
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>activation method</returns>
        private ActivationMethod GetActivationMethod(int scanNum)
        {
            try
            {
                var scanFilter = GetScanFilter(scanNum);
                var reactions = scanFilter.MassCount;

                var index = reactions - 1;
                if (index > 0 && scanFilter.GetIsMultipleActivation(index))
                {
                    // The last activation is part of a ETciD/EThcD pair
                    index--;
                }

                var activationTypeCode = scanFilter.GetActivation(index);

                switch (activationTypeCode)
                {
                    case ActivationType.CollisionInducedDissociation:
                        return ActivationMethod.CID;

                    case ActivationType.ElectronCaptureDissociation:
                        return ActivationMethod.ECD;

                    case ActivationType.PQD:
                        return ActivationMethod.PQD;

                    case ActivationType.ElectronTransferDissociation:
                        return ActivationMethod.ETD;

                    case ActivationType.HigherEnergyCollisionalDissociation:
                        return ActivationMethod.HCD;
                    case ActivationType.UltraVioletPhotoDissociation:
                        return ActivationMethod.UVPD;
                    default:
                        return ActivationMethod.Unknown;
                }
            }
            catch (Exception)
            {
                // TODO: Should report something somehow...
                return ActivationMethod.Unknown;
            }
        }

        private IScanFilter GetScanFilter(int scanNum)
        {
            if (_cachedScanFilter == null || _cachedScanFilterScanNumber != scanNum)
            {
                _cachedScanFilter = _rawFileReader.GetFilterForScanNumber(scanNum);
                _cachedScanFilterScanNumber = scanNum;
            }

            return _cachedScanFilter;
        }

        /// <summary>
        /// Close the reader
        /// </summary>
        public void Close()
        {
            _rawFileReader?.Dispose();
        }

        /// <summary>
        /// Clean up - close the file handle
        /// </summary>
        public void Dispose()
        {
            _rawFileReader?.Dispose();
        }
    }
}
