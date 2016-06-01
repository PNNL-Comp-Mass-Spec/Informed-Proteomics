using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.FeatureFinding.Data;
using InformedProteomics.FeatureFinding.Util;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.FeatureFinding.FeatureDetection
{
    public class LcMsPeakMatrixLowResolution
    {
        public LcMsPeakMatrixLowResolution(LcMsRun run, int maxFeaturesPerSpec = 100, double tolerancePpm = 64, int minScanCharge = 20, int maxScanCharge = 100, double minScanMass = 20000, double maxScanMass = 60000)
        {
            Run = run;
            _tolerance = new Tolerance(Math.Max(tolerancePpm, 32));
            _minCharge = minScanCharge;
            _maxCharge = maxScanCharge;
            _minMass = minScanMass;
            _maxMass = maxScanMass;
            _ms1Features= new NodeSet<Ms1Feature>();

            var nBits = (int)Math.Round(27 - Math.Log(_tolerance.GetValue() / 16, 2));
            _comparer = new MzComparerWithBinning(nBits);

            _maxFeaturesPerSpec = maxFeaturesPerSpec;
        }

        public List<LcMsFeature> GetLcMsFeatures()
        {
            return ClusterMs1Features();
        }

        public List<Ms1Feature> DetectMs1Features(int scanNum)
        {
            var features = new List<Ms1Feature>();
            var spec = Run.GetSpectrum(scanNum);
            var deconvoluter = new MaxEntDeconvoluter(spec, _tolerance, _comparer, _minCharge, _maxCharge, _minMass, _maxMass);
            var maxInt = -1.0d;


            foreach (var feature in deconvoluter.GetMassFeatures())
            {
                if (features.Count >= _maxFeaturesPerSpec) break;
                if (maxInt < 0) maxInt = feature.Abundance;
                if (feature.Abundance < AbundanceRatioCutoff*maxInt) break;
                features.Add(feature);
            }

            lock (_ms1Features)
            {
                _ms1Features.AddRange(features);
            }

            return features;
        }

        private List<LcMsFeature> ClusterMs1Features()
        {
            var clusters = _ms1Features.ConnnectedComponents(new Ms1FeatureComparer(Run, _tolerance));
            var ret = new List<LcMsFeature>();

            foreach (var cluster in clusters)
            {
                if (cluster.Count < 3) continue;

                var maxCharge = cluster.Select(m => m.MaxCharge).Max();
                var minCharge = cluster.Select(m => m.MinCharge).Min();
                var mass = cluster.Select(m => m.Mass).Median();
                var maxScan = cluster.Select(m => m.ScanNum).Max();
                var minScan = cluster.Select(m => m.ScanNum).Min();
                minScan = Run.GetPrevScanNum(minScan, 1);
                maxScan = Run.GetNextScanNum(maxScan, 1);

                var repScanNum = (int)cluster.Select(m => (double)m.ScanNum).Median();
                var repCharge = (int)(0.5*(minCharge + maxCharge));
                var abundance = cluster.Sum(m => m.Abundance);
                var minElution = Run.GetElutionTime(minScan);
                var maxElution = Run.GetElutionTime(maxScan);

                var repMz = (mass/repCharge) + Constants.Proton;

                ret.Add(new LcMsFeature(mass, repCharge, repMz, repScanNum, abundance, minCharge, maxCharge, minScan, maxScan, minElution, maxElution));
            }
            return ret;
        }

        public LcMsRun Run;
        private const double AbundanceRatioCutoff = 0.3;
        private readonly MzComparerWithBinning _comparer;
        private readonly Tolerance _tolerance;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly double _minMass;
        private readonly double _maxMass;
        private readonly int _maxFeaturesPerSpec;

        private readonly NodeSet<Ms1Feature> _ms1Features;
        internal class Ms1FeatureComparer : INodeComparer<Ms1Feature>
        {
            internal Ms1FeatureComparer(LcMsRun run, Tolerance tol)
            {
                _run = run;
                _tolerance = tol;
            }

            public bool SameCluster(Ms1Feature node1, Ms1Feature node2)
            {
                var massDiff = Math.Abs(node1.Mass - node2.Mass);
                if (massDiff > _tolerance.GetToleranceAsTh(node1.Mass)) return false;

                var etDiff = Math.Abs(_run.GetElutionTime(node1.ScanNum) - _run.GetElutionTime(node2.ScanNum));
                if (etDiff > ElutionTimeWindowSize) return false;

                return true;
            }

            private const double ElutionTimeWindowSize = 0.4;
            private readonly LcMsRun _run;
            private readonly Tolerance _tolerance;
        }

    }
}
