using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace ProMexAlign
{
    /// <summary>
    /// This class is used when a ProMex .ms1ft file is available, but the .raw or .pbf file for the dataset is not available
    /// </summary>
    internal class VirtualLcMsRun : LcMsRun
    {
        public sealed override double MinMs1Mz { get; }
        public sealed override double MaxMs1Mz { get; }

        public sealed override string FilePath { get; protected set; }
        public sealed override string SrcFileChecksum { get; protected set; }
        public sealed override string FileFormatVersion { get; protected set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="specFilePath"></param>
        /// <param name="minLcScan"></param>
        /// <param name="maxLcScan"></param>
        /// <param name="startTimeMinutes"></param>
        /// <param name="endTimeMinutes"></param>
        public VirtualLcMsRun(string specFilePath, int minLcScan = 1, int maxLcScan = 5000, double startTimeMinutes = 0, double endTimeMinutes = 120)
        {
            FilePath = specFilePath;
            SrcFileChecksum = string.Empty;
            FileFormatVersion = "1";

            MinMs1Mz = 100;
            MaxMs1Mz = 2000;

            MinLcScan = Math.Max(1, minLcScan);
            MaxLcScan = Math.Max(minLcScan, maxLcScan);

            PopulateDictionaries(startTimeMinutes, endTimeMinutes);
        }

        private void PopulateDictionaries(double startTimeMinutes, double endTimeMinutes)
        {
            ScanNumToMsLevel = new Dictionary<int, int>();

            ScanNumElutionTimeMap = new Dictionary<int, double>();

            var acqTimeMinutes = endTimeMinutes - startTimeMinutes;
            var scanCount = MaxLcScan - MinLcScan + 1;

            for (var scan = MinLcScan; scan <= MaxLcScan; scan++)
            {
                ScanNumToMsLevel.Add(scan, 1);

                var elutionTime = acqTimeMinutes * (scan - MinLcScan) / scanCount;
                ScanNumElutionTimeMap.Add(scan, elutionTime);
            }

            CreatePrecursorNextScanMap();
        }

        public override void Close()
        {
        }

        public override bool TryMakeRandomAccessCapable()
        {
            return true;
        }

        public override Spectrum GetSpectrum(int scanNum, bool includePeaks = true)
        {
            return null;
        }

        public override Spectrum GetMs1Spectrum(int scanNum, out int ms1ScanIndex)
        {
            var ms1ScanNums = GetMs1ScanVector();
            ms1ScanIndex = Array.BinarySearch(ms1ScanNums, scanNum);

            if (ms1ScanIndex < 1)
                ms1ScanIndex = 0;

            // Uncomment to create a placeholder spectrum with three data points
            //var mzList = new List<double> {100, 500, 1000};

            //var intensityList = new List<double> {750, 1000, 750};

            //return new Spectrum(mzList, intensityList, scanNum);

            return null;
        }

        public override IsolationWindow GetIsolationWindow(int scanNum)
        {
            return null;
        }

        public override Xic GetPrecursorExtractedIonChromatogram(double minMz, double maxMz)
        {
            return new();
        }

        public override Xic GetPrecursorChromatogramRange(double minMz, double maxMz)
        {
            return new();
        }

        public override Xic GetFullProductExtractedIonChromatogram(double minMz, double maxMz, double precursorMz)
        {
            return new();
        }

        public override void Dispose()
        {
        }
    }
}
