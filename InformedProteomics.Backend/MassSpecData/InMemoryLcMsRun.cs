using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.MassSpecData
{
    public class InMemoryLcMsRun : LcMsRun //LcMsRun 
    {
        [ObsoleteAttribute("Remove MassSpecDataType -> now uses MassSpecDataReaderFactory", true)]
        public static LcMsRun GetLcMsRun(string specFilePath, MassSpecDataType dataType, IProgress<ProgressData> progress = null)
        {
            return GetLcMsRun(specFilePath, dataType, 0.0, 0.0, progress);
        }

        public static LcMsRun GetLcMsRun(string specFilePath, IProgress<ProgressData> progress = null)
        {
            return GetLcMsRun(specFilePath, MassSpecDataReaderFactory.GetMassSpecDataReader(specFilePath), 0.0, 0.0, progress);
        }

        [ObsoleteAttribute("Remove MassSpecDataType -> now uses MassSpecDataReaderFactory", true)]
        public static LcMsRun GetLcMsRun(string specFilePath, MassSpecDataType dataType, double precursorSignalToNoiseRatioThreshold,
            double productSignalToNoiseRatioThreshold, IProgress<ProgressData> progress = null)
        {
            //var pbfFilePath = ConvertToPbf(specFilePath, dataType, precursorSignalToNoiseRatioThreshold,
            //    productSignalToNoiseRatioThreshold, null, progress);
            //
            //return new InMemoryLcMsRun(new PbfLcMsRun(pbfFilePath), precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold, progress);
            return GetLcMsRun(specFilePath, MassSpecDataReaderFactory.GetMassSpecDataReader(specFilePath), precursorSignalToNoiseRatioThreshold,
                productSignalToNoiseRatioThreshold, progress);
        }

        public static LcMsRun GetLcMsRun(string specFilePath, double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold,
            IProgress<ProgressData> progress = null)
        {
            return GetLcMsRun(specFilePath, MassSpecDataReaderFactory.GetMassSpecDataReader(specFilePath), precursorSignalToNoiseRatioThreshold,
                productSignalToNoiseRatioThreshold, progress);
        }

        public static LcMsRun GetLcMsRun(string specFilePath, IMassSpecDataReader specReader,
            double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold, IProgress<ProgressData> progress = null)
        {
            var pbfFilePath = ConvertToPbf(specFilePath, specReader, precursorSignalToNoiseRatioThreshold,
                productSignalToNoiseRatioThreshold, null, progress);

            return new InMemoryLcMsRun(new PbfLcMsRun(pbfFilePath), precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold, progress);
        }

        public static string ConvertToPbf(string specFilePath, double precursorSignalToNoiseRatioThreshold,
            double productSignalToNoiseRatioThreshold, string pbfFilePath = null, IProgress<ProgressData> progress = null)
        {
            return ConvertToPbf(specFilePath, MassSpecDataReaderFactory.GetMassSpecDataReader(specFilePath),
                precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold, pbfFilePath, progress);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="specFilePath"></param>
        /// <param name="specReader"></param>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        /// <param name="pbfFilePath">If supplied, file will be written to this path; otherwise the file will be written to the same directory as specFilePath, or to the temp directory if the user does not have write permissions</param>
        /// <param name="progress">Progress data, as a percentage</param>
        /// <returns></returns>
        public static string ConvertToPbf(string specFilePath, IMassSpecDataReader specReader,
            double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold, string pbfFilePath = null,
            IProgress<ProgressData> progress = null)
        {
            string pbfPath = pbfFilePath;
            string fileName = String.Empty;
            string tempPath = String.Empty;

            Progress<ProgressData> prog = new Progress<ProgressData>();
            var progData = new ProgressData();
            progData.IsPartialRange = true;
            progData.MaxPercentage = 75.0;
            if (progress != null)
            {
                prog = new Progress<ProgressData>(p =>
                {
                    progData.Status = p.Status;
                    progress.Report(progData.UpdatePercent(p.Percent));
                });
            }

            if (string.IsNullOrWhiteSpace(pbfFilePath))
            {
                pbfPath = PbfLcMsRun.GetPbfFileName(specFilePath);
                fileName = Path.GetFileName(pbfPath);
                if (String.IsNullOrEmpty(fileName))
                {
                    throw new ArgumentException("Cannot create .pbf cache file", "specFilePath");
                }

                tempPath = Path.Combine(Path.GetTempPath(), fileName);
                // Return the temp path if the pbf file of proper format already exists in the temp directory
                if (File.Exists(tempPath) && PbfLcMsRun.CheckFileFormatVersion(tempPath))
                {
                    return tempPath;
                }
            }

            if (!File.Exists(pbfPath) || !PbfLcMsRun.CheckFileFormatVersion(pbfPath))
            {
                if (specReader == null)
                {
                    throw new Exception("Unsupported file format!");
                }
                InMemoryLcMsRun run = new InMemoryLcMsRun(specReader, 0, 0, prog);
                try
                {
                    progData.StepRange(100.0);
                    run.WriteAsPbf(pbfPath, prog);
                }
                catch (UnauthorizedAccessException) // Cannot write to same directory, attempt to write to temp directory
                {
                    if (string.IsNullOrWhiteSpace(pbfFilePath))
                    {
                        throw;
                    }
                    //var fileName = Path.GetFileName(pbfFilePath);
                    if (String.IsNullOrEmpty(fileName)) throw; // invalid path?
                    //var tempPath = Path.Combine(Path.GetTempPath(), fileName);
                    if (!File.Exists(tempPath) || !PbfLcMsRun.CheckFileFormatVersion(tempPath))
                        run.WriteAsPbf(tempPath, prog);
                    pbfPath = tempPath;
                }
            }
            return pbfPath;
        }

        public InMemoryLcMsRun(IMassSpecDataReader massSpecDataReader, double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold, IProgress<ProgressData> progress = null)
        {
            ScanNumElutionTimeMap = new Dictionary<int, double>();
            ScanNumToMsLevel = new Dictionary<int, int>();
            IsolationMzBinToScanNums = new Dictionary<int, int[]>();

            _ms1PeakList = new List<LcMsPeak>();
            _scanNumSpecMap = new Dictionary<int, Spectrum>();

            var isolationMzBinToScanNums = new Dictionary<int, List<int>>();

            // Read all spectra
            var minScanNum = int.MaxValue;
            var maxScanNum = int.MinValue;
            var minMsLevel = int.MaxValue;
            var maxMsLevel = int.MinValue;

            if (progress == null)
            {
                progress = new Progress<ProgressData>();
            }
            var progressData = new ProgressData();
            progressData.Status = "Reading spectra from file";

            NumSpectra = massSpecDataReader.NumSpectra;
            int specRead = 0;
            progressData.IsPartialRange = true;
            progressData.MaxPercentage = 95.0;

            foreach (var spec in massSpecDataReader.ReadAllSpectra())
            {
                progress.Report(progressData.UpdatePercent((double)specRead / NumSpectra * 100.0));
                specRead++;
                //Console.WriteLine("Reading Scan {0}; {1} peaks", spec.ScanNum, spec.Peaks.Length);
                ScanNumToMsLevel[spec.ScanNum] = spec.MsLevel;
                ScanNumElutionTimeMap[spec.ScanNum] = spec.ElutionTime;
                if (spec.MsLevel == 1)
                {
                    if (precursorSignalToNoiseRatioThreshold > 0.0)
                        spec.FilterNoise(precursorSignalToNoiseRatioThreshold);

                    //foreach (var peak in spec.Peaks)
                    //{
                    //    _ms1PeakList.Add(new LcMsPeak(peak.Mz, peak.Intensity, spec.ScanNum));
                    //}
                    _ms1PeakList.AddRange(spec.Peaks.Select(peak => new LcMsPeak(peak.Mz, peak.Intensity, spec.ScanNum)));
                    _scanNumSpecMap.Add(spec.ScanNum, spec);
                }
                else if (spec.MsLevel == 2)
                {
                    var productSpec = spec as ProductSpectrum;

                    if (productSpec != null)
                    {
                        if (productSignalToNoiseRatioThreshold > 0.0)
                            productSpec.FilterNoise(productSignalToNoiseRatioThreshold);

                        var isolationWindow = productSpec.IsolationWindow;
                        var minBinNum = (int)Math.Round(isolationWindow.MinMz * IsolationWindowBinningFactor);
                        var maxBinNum = (int)Math.Round(isolationWindow.MaxMz * IsolationWindowBinningFactor);
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

            progressData.Status = "Processing Isolation Bins";
            progressData.IsPartialRange = false;
            progress.Report(progressData.UpdatePercent(95.1));

            foreach (var entry in isolationMzBinToScanNums)
            {
                var binNum = entry.Key;
                entry.Value.Sort();
                var scanNumList = entry.Value.ToArray();
                IsolationMzBinToScanNums[binNum] = scanNumList;
            }

            _ms1PeakList.Sort();
            //_ms2PeakList.Sort();

            progress.Report(progressData.UpdatePercent(99.5));
            // Read MS levels and precursor information

            MinLcScan = minScanNum;
            MaxLcScan = maxScanNum;

            MinMsLevel = minMsLevel;
            MaxMsLevel = maxMsLevel;

            //var precursorMap = new Dictionary<int, int>();
            //var nextScanMap = new Dictionary<int, int>();
            //
            //for (var msLevel = MinMsLevel; msLevel <= maxMsLevel; msLevel++)
            //{
            //    precursorMap[msLevel] = 0;
            //    nextScanMap[msLevel] = MaxLcScan + 1;
            //}
            //progress.Report(progressData.UpdatePercent(99.8));

            progress.Report(progressData.UpdatePercent(100.0));
            CreatePrecursorNextScanMap();
        }

        public List<LcMsPeak> Ms1PeakList { get { return _ms1PeakList; } }

        public override double MinMs1Mz { get { return _ms1PeakList[0].Mz; } }
        public override double MaxMs1Mz { get { return _ms1PeakList[_ms1PeakList.Count - 1].Mz; } }

        /// <summary>
        /// Gets the spectrum of the specified scan number
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>spectrum</returns>
        public override Spectrum GetSpectrum(int scanNum)
        {
            Spectrum spec;
            return _scanNumSpecMap.TryGetValue(scanNum, out spec) ? spec : null;
        }

        public override Ms1Spectrum GetMs1Spectrum(int scanNum)
        {
            //throw new NotImplementedException();
            var spec = GetSpectrum(scanNum);
            var ms1ScanNums = GetMs1ScanVector();
            var ms1ScanIndex = Array.BinarySearch(ms1ScanNums, scanNum);
            if (ms1ScanIndex < 0) return null;
            return new Ms1Spectrum(scanNum, (ushort)ms1ScanIndex, spec.Peaks);
        }

        public override IsolationWindow GetIsolationWindow(int scanNum)
        {
            var spec = GetSpectrum(scanNum);
            if (spec == null) return null;

            var productSpec = GetSpectrum(scanNum) as ProductSpectrum;
            return productSpec == null ? null : productSpec.IsolationWindow;
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS1 spectra)
        /// </summary>
        /// <param name="minMz">min m/z</param>
        /// <param name="maxMz">max m/z</param>
        /// <returns>XIC as an Xic object</returns>
        public override Xic GetPrecursorExtractedIonChromatogram(double minMz, double maxMz)
        {
            return GetExtractedIonChromatogram(minMz, maxMz, Ms1PeakList);
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS2 spectra)
        /// </summary>
        /// <param name="minMz">min m/z</param>
        /// <param name="maxMz">max m/z</param>
        /// <param name="precursorIonMz">precursor m/z of the precursor ion</param>
        /// <returns>XIC as an Xic object</returns>
        public override Xic GetFullProductExtractedIonChromatogram(double minMz, double maxMz, double precursorIonMz)
        {
            return GetProductExtractedIonChromatogram(minMz, maxMz, precursorIonMz, MinLcScan, MaxLcScan);
        }


        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS2 spectra)
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <param name="precursorIonMz">precursor m/z of the precursor ion</param>
        /// <param name="minScanNum">minimum scan number (inclusive)</param>
        /// <param name="maxScanNum">maximum scan number (inclusive)</param>
        /// <returns>XIC as an Xic object</returns>
        public Xic GetProductExtractedIonChromatogram(double mz, Tolerance tolerance, double precursorIonMz, int minScanNum, int maxScanNum)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetProductExtractedIonChromatogram(minMz, maxMz, precursorIonMz, minScanNum, maxScanNum);
        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS2 spectra)
        /// </summary>
        /// <param name="minMz">min m/z</param>
        /// <param name="maxMz">max m/z</param>
        /// <param name="precursorIonMz">precursor m/z of the precursor ion</param>
        /// <param name="minScanNum">minimum scan number (inclusive)</param>
        /// <param name="maxScanNum">maximum scan number (inclusive)</param>
        /// <returns>XIC as an Xic object</returns>
        public Xic GetProductExtractedIonChromatogram(double minMz, double maxMz, double precursorIonMz, int minScanNum, int maxScanNum)
        {
            var xic = new Xic();
            //for (var scanNum = minScanNum; scanNum <= maxScanNum; scanNum++)
            foreach (var scanNum in GetFragmentationSpectraScanNums(precursorIonMz))
            {
                if (scanNum < minScanNum || scanNum > maxScanNum) continue;
                var spec = GetSpectrum(scanNum) as ProductSpectrum;
                if (spec == null) continue;

                var peak = spec.FindPeak(minMz, maxMz);
                xic.Add(peak != null ? new XicPoint(scanNum, peak.Mz, peak.Intensity) : new XicPoint(scanNum, 0, 0));
            }
            return xic;
        }

        public Xic GetProductExtractedIonChromatogram(double productIonMz, double precursorIonMz, Tolerance tolerance, int minScanNum, int maxScanNum)
        {
            var productXic = new Xic();

            for (var scanNum = minScanNum; scanNum <= maxScanNum; scanNum++)
            {
                if (GetMsLevel(scanNum) == 1) continue;

                var productSpec = _scanNumSpecMap[scanNum] as ProductSpectrum;
                if (productSpec == null) continue;

                if (!productSpec.IsolationWindow.Contains(precursorIonMz)) continue;

                var peak = productSpec.FindPeak(productIonMz, tolerance);
                if (peak != null) productXic.Add(new XicPoint(scanNum, peak.Mz, peak.Intensity));
            }

            return productXic;
        }

        public void WriteAsPbf(string outputFilePath, IProgress<ProgressData> progress = null)
        {
            using (var writer = new BinaryWriter(File.Open(outputFilePath, FileMode.Create)))
            {
                WriteAsPbf(writer, progress);
            }
        }

        public void WriteAsPbf(BinaryWriter writer, IProgress<ProgressData> progress = null)
        {
            long countTotal = 1;
            long counter = 0;
            if (progress == null)
            {
                progress = new Progress<ProgressData>();
            }
            var progressData = new ProgressData();
            progressData.IsPartialRange = true;
            progressData.Status = "Writing spectra data";

            var scanNumToSpecOffset = new long[MaxLcScan - MinLcScan + 1];
            var scanNumToIsolationWindow = new IsolationWindow[MaxLcScan - MinLcScan + 1];

            // Spectra
            countTotal = MaxLcScan - MinLcScan;
            counter = 0;
            progressData.MaxPercentage = 42.9; // SpecData: Approximately 43% of total file size
            long countMS2Spec = 0;
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                progress.Report(progressData.UpdatePercent((double)counter / countTotal * 100.0));
                counter++;
                scanNumToSpecOffset[scanNum - MinLcScan] = writer.BaseStream.Position;
                var spec = GetSpectrum(scanNum);
                if (spec == null) continue;
                var productSpec = spec as ProductSpectrum;
                scanNumToIsolationWindow[scanNum - MinLcScan] = null;
                if (productSpec != null)
                {
                    scanNumToIsolationWindow[scanNum - MinLcScan] = productSpec.IsolationWindow;
                    countMS2Spec++;
                }
                WriteSpectrum(spec, writer);
            }

            // Precursor ion chromatograms
            var offsetBeginPrecursorChromatogram = writer.BaseStream.Position;

            var minMzIndex = Ms1PeakList.Any() ? PbfLcMsRun.GetMzBinIndex(Ms1PeakList[0].Mz) : 0;
            var maxMzIndex = Ms1PeakList.Any() ? PbfLcMsRun.GetMzBinIndex(Ms1PeakList[Ms1PeakList.Count - 1].Mz) : -1;

            var chromMzIndexToOffset = new long[maxMzIndex - minMzIndex + 1];
            var prevMzIndex = -1;
            counter = 0;
            countTotal = Ms1PeakList.Count;
            progressData.Status = "Writing precursor chromatograms";
            progressData.StepRange(42.9 + 15.7); // Approximately 16% of total file size
            foreach (var peak in Ms1PeakList)
            {
                progress.Report(progressData.UpdatePercent((double)counter / countTotal * 100.0));
                counter++;
                var mz = peak.Mz;
                var mzIndex = PbfLcMsRun.GetMzBinIndex(mz);
                if (mzIndex > prevMzIndex)
                {
                    chromMzIndexToOffset[mzIndex - minMzIndex] = writer.BaseStream.Position;
                    prevMzIndex = mzIndex;
                }
                writer.Write(peak.Mz);
                writer.Write((float)peak.Intensity);
                writer.Write(peak.ScanNum);
            }

            // Product ion chromatograms
            var ms2PeakList = new List<LcMsPeak>();
            counter = 0;
            countTotal = countMS2Spec;
            progressData.Status = "Processing product chromatograms";
            progressData.StepRange(42.9 + 15.7 + (41.2 / 2)); // Approximately 41% of total file size
            foreach (var ms2ScanNum in GetScanNumbers(2))
            {
                progress.Report(progressData.UpdatePercent((double)counter / countTotal * 100.0));
                counter++;
                var productSpec = GetSpectrum(ms2ScanNum) as ProductSpectrum;
                if (productSpec == null) continue;
                foreach (var peak in productSpec.Peaks)
                {
                    ms2PeakList.Add(new LcMsPeak(peak.Mz, peak.Intensity, ms2ScanNum));
                }
            }
            ms2PeakList.Sort();

            var offsetBeginProductChromatogram = writer.BaseStream.Position;
            counter = 0;
            countTotal = ms2PeakList.Count;
            progressData.Status = "Writing product chromatograms";
            progressData.StepRange(42.9 + 15.7 + 41.2); // Approximately 41% of total file size
            foreach (var peak in ms2PeakList)
            {
                progress.Report(progressData.UpdatePercent((double)counter / countTotal * 100.0));
                counter++;
                writer.Write(peak.Mz);
                writer.Write((float)peak.Intensity);
                writer.Write(peak.ScanNum);
            }

            // Meta information
            var offsetBeginMetaInformation = writer.BaseStream.Position;
            progressData.Status = "Writing metadata";
            progressData.IsPartialRange = false;
            progress.Report(progressData.UpdatePercent(99.8)); // Metadata: Approximately 0.2% of total file size

            var warnedInvalidScanNum = false;
            var warnedNullScanToIsolationWindow = false;

            writer.Write(MinLcScan);
            writer.Write(MaxLcScan);
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                var msLevel = GetMsLevel(scanNum);
                writer.Write(GetMsLevel(scanNum));
                writer.Write(GetElutionTime(scanNum));
             
                if (msLevel == 2)
                {
                    float minMz = 0;
                    float maxMz = 0;

                    if (scanNum - MinLcScan < 0 || scanNum - MinLcScan >= scanNumToIsolationWindow.Length)
                    {
                        if (!warnedInvalidScanNum)
                        {
                            Console.WriteLine("\nWriteAsPbf encountered an invalid scan number: " + scanNum + "; " +
                                              "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown");
                            warnedInvalidScanNum = true;
                        }
                    }
                    else
                    {
                        if (scanNumToIsolationWindow[scanNum - MinLcScan] == null)
                        {
                            if (!warnedNullScanToIsolationWindow)
                            {
                                Console.WriteLine("\nWriteAsPbf encountered a Null entry in scanNumToIsolationWindow for scan " + scanNum + "; " +
                                                  "MinMz and MaxMz will be 0 for this scan; subsequent warnings of this type will not be shown");
                                warnedNullScanToIsolationWindow = true;
                            }
                        }
                        else
                        {
                            minMz = (float)scanNumToIsolationWindow[scanNum - MinLcScan].MinMz;
                            maxMz = (float)scanNumToIsolationWindow[scanNum - MinLcScan].MaxMz;
                        }
                    }

                    writer.Write(minMz);
                    writer.Write(maxMz);
                }
                writer.Write(scanNumToSpecOffset[scanNum - MinLcScan]);
            }

            // Precursor chromatogram index
            writer.Write(minMzIndex);   // min index
            writer.Write(maxMzIndex);
            progress.Report(progressData.UpdatePercent(99.9)); // Metadata: Approximately 0.2% of total file size

            var prevOffset = offsetBeginMetaInformation;
            for (var i = chromMzIndexToOffset.Length - 1; i >= 0; i--)
            {
                if (chromMzIndexToOffset[i] < offsetBeginPrecursorChromatogram) chromMzIndexToOffset[i] = prevOffset;
                else prevOffset = chromMzIndexToOffset[i];
            }

            foreach (var offset in chromMzIndexToOffset)
            {
                writer.Write(offset);
            }

            writer.Write(offsetBeginPrecursorChromatogram); // 8
            writer.Write(offsetBeginProductChromatogram); // 8
            writer.Write(offsetBeginMetaInformation); // 8
            progress.Report(progressData.UpdatePercent(100.0));
            writer.Write(PbfLcMsRun.FileFormatId); // 4
        }

        private static Xic GetXicPointsWithin(double minMz, double maxMz, List<LcMsPeak> peakList)
        {
            var xic = new Xic();

            var index = peakList.BinarySearch(new LcMsPeak((minMz + maxMz) / 2, 0, 0));
            if (index < 0) index = ~index;

            // go down
            var i = index - 1;
            while (i >= 0 && i < peakList.Count)
            {
                var peak = peakList[i];
                if (peak.Mz <= minMz) break;
                xic.Add(new XicPoint(peak.ScanNum, peak.Mz, peak.Intensity));
                --i;
            }

            // go up
            i = index;
            while (i >= 0 && i < peakList.Count)
            {
                var peak = peakList[i];
                if (peak.Mz >= maxMz) break;
                xic.Add(new XicPoint(peak.ScanNum, peak.Mz, peak.Intensity));
                ++i;
            }

            return xic;
        }

        private static Xic GetExtractedIonChromatogram(double minMz, double maxMz, List<LcMsPeak> peakList)
        {
            var xic = GetXicPointsWithin(minMz, maxMz, peakList);
            xic.Sort();
            return Xic.GetSelectedXic(xic);
        }

        private readonly List<LcMsPeak> _ms1PeakList;
        private readonly Dictionary<int, Spectrum> _scanNumSpecMap;  // scan number -> spectrum
    }

    public class LcMsPeak : Peak
    {
        public LcMsPeak(double mz, double intensity, int scanNum)
            : base(mz, intensity)
        {
            ScanNum = scanNum;
        }
        public int ScanNum { get; private set; }
    }
}
