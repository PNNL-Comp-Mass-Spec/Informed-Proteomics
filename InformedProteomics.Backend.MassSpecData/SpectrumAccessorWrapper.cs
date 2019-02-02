using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using PRISM;

namespace InformedProteomics.Backend.MassSpecData
{
    public class SpectrumAccessorWrapper : ISpectrumAccessor, IDisposable
    {
        private readonly IMassSpecDataReader reader;
        private readonly Dictionary<int, ScanMetadata> scanMetadata = new Dictionary<int, ScanMetadata>();

        public SpectrumAccessorWrapper(IMassSpecDataReader msReader, IProgress<ProgressData> progress = null)
        {
            reader = msReader;
            MinLcScan = int.MaxValue;
            MaxLcScan = int.MinValue;
            NumSpectra = 0;

            var progressData = new ProgressData(progress);
            var expectedNumSpectra = reader.NumSpectra;
            progressData.Report(0, expectedNumSpectra);
            foreach (var spec in reader.ReadAllSpectra(false))
            {
                var metadata = new ScanMetadata(spec.MsLevel, spec.ElutionTime);
                scanMetadata.Add(spec.ScanNum, metadata);

                NumSpectra++;
                MinLcScan = Math.Min(MinLcScan, spec.ScanNum);
                MaxLcScan = Math.Max(MaxLcScan, spec.ScanNum);
                progressData.Report(NumSpectra, expectedNumSpectra);
            }
        }

        public Spectrum GetSpectrum(int scanNum, bool includePeaks = true)
        {
            if (!scanMetadata.ContainsKey(scanNum))
            {
                return null;
            }

            return reader.GetSpectrum(scanNum, includePeaks);
        }

        public double GetElutionTime(int scanNum)
        {
            return scanMetadata.TryGetValue(scanNum, out var metadata) ? metadata.ElutionTime : 0.0;
        }

        public int GetMsLevel(int scanNum)
        {
            return scanMetadata.TryGetValue(scanNum, out var metadata) ? metadata.MsLevel : 0;
        }

        public int GetPrevScanNum(int scanNum, int msLevel)
        {
            return scanMetadata.Where(x => x.Value.MsLevel == msLevel && x.Key < scanNum).Select(x => x.Key).DefaultIfEmpty(MinLcScan - 1).Max();
        }

        public int GetNextScanNum(int scanNum, int msLevel)
        {
            return scanMetadata.Where(x => x.Value.MsLevel == msLevel && x.Key > scanNum).Select(x => x.Key).DefaultIfEmpty(MaxLcScan + 1).Min();
        }

        public IList<int> GetScanNumbers(int msLevel)
        {
            return scanMetadata.Where(x => x.Value.MsLevel == msLevel).Select(x => x.Key).ToList();
        }

        public IsolationWindow GetIsolationWindow(int scanNum)
        {
            return (GetSpectrum(scanNum, false) as ProductSpectrum)?.IsolationWindow;
        }

        public int MinLcScan { get; }
        public int MaxLcScan { get; }
        public int NumSpectra { get; }

        public void Dispose()
        {
            reader?.Dispose();
        }

        private class ScanMetadata
        {
            public int MsLevel { get; }
            public double ElutionTime { get; }

            public ScanMetadata(int msLevel, double elutionTime)
            {
                MsLevel = msLevel;
                ElutionTime = elutionTime;
            }
        }
    }
}
