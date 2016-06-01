using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows.Forms;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public class DiaLcMsRun: LcMsRun
    {
        public DiaLcMsRun(IMassSpecDataReader massSpecDataReader, double precursorSignalToNoiseRatioThreshold, double productSignalToNoiseRatioThreshold)
            : base(massSpecDataReader, precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold)
        {
            _isolationWindows = new SortedSet<IsolationWindow>();
            _isoWindowToChromPeaks = new Dictionary<IsolationWindow, List<LcMsPeak>>();

            var variableWindows = false;
            var upper = -1.0;
            var lower = -1.0;
            for (var scanNum = MinLcScan; scanNum <= MaxLcScan; scanNum++)
            {
                var productSpec = GetSpectrum(scanNum) as ProductSpectrum;
                if (productSpec == null) continue;
                var isoWindow = productSpec.IsolationWindow;
                if (isoWindow == null) continue;

                if (_isolationWindows.Add(isoWindow))
                {
                    // new isolation window
                    if (!variableWindows)
                    {
                        if (upper > 0 && Math.Abs(upper - isoWindow.IsolationWindowUpperOffset) > 0.001)
                        {
                            variableWindows = true;
                        }
                        if (lower > 0 && Math.Abs(lower - isoWindow.IsolationWindowLowerOffset) > 0.001)
                        {
                            variableWindows = true;
                        }
                        if (upper < 0) upper = isoWindow.IsolationWindowUpperOffset;
                        if (lower < 0) lower = isoWindow.IsolationWindowLowerOffset;
                    }
                    _isoWindowToChromPeaks.Add(isoWindow, new List<LcMsPeak>());
                }
                var peakList = _isoWindowToChromPeaks[isoWindow];
                peakList.AddRange(productSpec.Peaks.Select(peak => new LcMsPeak(peak.Mz, peak.Intensity, scanNum)));
            }

            foreach (var peakList in _isoWindowToChromPeaks.Values)
            {
                peakList.Sort();
            }

            _variableWindows = variableWindows;
            _isolationWindowLowerOffset = lower;
            _isolationWindowUpperOffset = upper;
        }

        public Xic GetProductExtractedIonChromatogram(double minMz, double maxMz, double precursorIonMz)
        {
            var xic = new Xic();
            foreach (var isolationWindow in GetIsoWindows(precursorIonMz))
            {
                var peakList = _isoWindowToChromPeaks[isolationWindow];
                xic.AddRange(GetXicPointsWithin(minMz, maxMz, peakList));
            }
            xic.Sort();
            return Xic.GetSelectedXic(xic);
        }

        public new Xic GetFullProductExtractedIonChromatogram(double mz, Tolerance tolerance, double precursorIonMz)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;
            return GetFullProductExtractedIonChromatogram(minMz, maxMz, precursorIonMz);
        }

        public new Xic GetFullProductExtractedIonChromatogram(double minMz, double maxMz, double precursorIonMz)
        {
            var xic = GetProductExtractedIonChromatogram(minMz, maxMz, precursorIonMz);
            var scanToXicPoint = new XicPoint[MaxLcScan - MinLcScan + 1];
            foreach (var xicPoint in xic) scanToXicPoint[xicPoint.ScanNum - MinLcScan] = xicPoint;

            var newXic = new Xic();
            newXic.AddRange(GetFragmentationSpectraScanNums(precursorIonMz).Select(scanNum => scanToXicPoint[scanNum - MinLcScan] ?? new XicPoint(scanNum, 0, 0)));
            return newXic;
        }

        public new void WriteAsRaf(string outputFilePath)
        {
            using (var writer = new BinaryWriter(File.Open(outputFilePath, FileMode.Create)))
            {
                WriteAsPbf(writer);
            }
        }

        private IEnumerable<IsolationWindow> GetIsoWindows(double precursorMz)
        {
            var lower = new IsolationWindow(precursorMz - _isolationWindowUpperOffset, 0, 0);
            var higher = new IsolationWindow(precursorMz + _isolationWindowLowerOffset, 0, 0);
            return _isolationWindows.GetViewBetween(lower, higher);
        }

        private readonly bool _variableWindows;
        private readonly double _isolationWindowUpperOffset;
        private readonly double _isolationWindowLowerOffset;
        private readonly SortedSet<IsolationWindow> _isolationWindows;
        private readonly Dictionary<IsolationWindow, List<LcMsPeak>> _isoWindowToChromPeaks;
    }
}
