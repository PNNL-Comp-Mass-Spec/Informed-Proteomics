using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Ms1FeatureComparer : IComparer<Ms1FeatureCluster>
    {
        public int Compare(Ms1FeatureCluster x, Ms1FeatureCluster y)
        {
            return (y.Probability.CompareTo(x.Probability));
            /*
            if (y.GetScore(Ms1FeatureScore.BhattacharyyaDistance) < 0.05 ||
                x.GetScore(Ms1FeatureScore.BhattacharyyaDistance) < 0.05)
            {
                return
                    x.GetScore(Ms1FeatureScore.BhattacharyyaDistance)
                        .CompareTo(y.GetScore(Ms1FeatureScore.BhattacharyyaDistance));
            }
            return
                x.GetScore(Ms1FeatureScore.BhattacharyyaDistanceSummed)
                    .CompareTo(y.GetScore(Ms1FeatureScore.BhattacharyyaDistanceSummed));
            */
        }
    }
    
    public class Ms1FeatureContainer
    {
        public Ms1FeatureContainer(List<Ms1Spectrum> ms1Spectra)
        {
            //_peakToFeatureMap = new Dictionary<Ms1Peak, List<Ms1FeatureCluster>>();
            _featureList = new List<Ms1FeatureCluster>();
            _spectrums = ms1Spectra;
        }

        public bool Add(Ms1FeatureCluster newFeature)
        {
            
            for (var i = _featureList.Count - 1; i >= 0; i--)
            {
                if (!_featureList[i].CoEluted(newFeature)) continue;
                var massDiff = Math.Abs(_featureList[i].RepresentativeMass - newFeature.RepresentativeMass);

                //if (massDiff > 1.5) break;
                //var massDiffPpm = (1e6 * massDiff) / newFeature.RepresentativeMass;
                //if (massDiffPpm < 2.5)
                if (massDiff < 1e-6)
                {
                    // already exists, then skip!
                    _featureList[i].Merge(newFeature);
                    return false;
                }
            }
            /*
            
            for (var i = _featureList.Count - 1; i >= 0; i--)
            {
                var massDiff = Math.Abs(_featureList[i].RepresentativeMass - newFeature.RepresentativeMass);
                var massDiffPpm = (1e6*massDiff) / newFeature.RepresentativeMass;
                if (massDiffPpm > 20) break;
                if (massDiff < 1e-8)
                {
                    if ((_featureList[i].MinCol <= newFeature.MinCol && newFeature.MinCol <= _featureList[i].MaxCol) ||
                        (_featureList[i].MinCol <= newFeature.MaxCol && newFeature.MaxCol <= _featureList[i].MaxCol))
                    {
                        var set1 = new HashSet<Ms1Peak>(newFeature.GetMajorPeaks());
                        set1.ExceptWith(_featureList[i].GetMajorPeaks());

                        var set2 = new HashSet<Ms1Peak>(newFeature.GetMinorPeaks());
                        set2.ExceptWith(_featureList[i].GetMinorPeaks());

                        _featureList[i].Merge(newFeature);
                        _featureList[i].UpdateScores(_spectrums);

                        foreach (var peak in set1) peak.TagMajorPeakOf(_featureList[i]);
                        foreach (var peak in set2) peak.TagMinorPeakOf(_featureList[i]);
                        
                        // already exists, then skip!
                        return false;
                    }
                }
            }*/
            
            foreach (var peak in newFeature.GetMajorPeaks()) peak.TagMajorPeakOf(newFeature);

            foreach (var peak in newFeature.GetMinorPeaks()) peak.TagMinorPeakOf(newFeature);
            
            _featureList.Add(newFeature);
            return true;
        }
        
        public void Add(IEnumerable<Ms1FeatureCluster> features)
        {
            foreach(var f in features) Add(f);
        }

        public IEnumerable<Ms1FeatureCluster> GetFilteredFeatures(IList<SortedSet<Ms1FeatureCluster>> connectedFeatureList)
        {
            var filteredFeatures = new List<Ms1FeatureCluster>();
            foreach (var featureSet in connectedFeatureList)
            {
                var newList = RemoveOverlappedFeatures(featureSet);
                foreach (var f in newList) f.UpdateAbundance();
                filteredFeatures.AddRange(newList);
            }
            return filteredFeatures.OrderBy(f => f.RepresentativeMass);
        }


        private bool SimilarScore(Ms1FeatureCluster f1, Ms1FeatureCluster f2)
        {
            if (Math.Abs(f1.GetScore(Ms1FeatureScore.EnvelopeCorrelation) -
                         f2.GetScore(Ms1FeatureScore.EnvelopeCorrelation)) > 0.05) return false;

            if (Math.Abs(f1.GetScore(Ms1FeatureScore.BhattacharyyaDistance) -
                                     f2.GetScore(Ms1FeatureScore.BhattacharyyaDistance)) > 0.01) return false;
            return true;
        }
        

        private IList<Ms1FeatureCluster> RemoveOverlappedFeatures(SortedSet<Ms1FeatureCluster> featureSet)
        {
            var outFeatures = new List<Ms1FeatureCluster>();
            var tol = new Tolerance(5);
            while (true)
            {
                if (featureSet.Count < 1) break;

                var bestFeature = featureSet.First();
                featureSet.Remove(bestFeature);
                outFeatures.Add(bestFeature);
                var massTol = tol.GetToleranceAsTh(bestFeature.RepresentativeMass);
                
                var tempList = new List<Ms1FeatureCluster>();
                foreach (var f in bestFeature.OverlappedFeatures)
                {
                    if (featureSet.Remove(f))
                    {
                        var massDiff = Math.Abs(bestFeature.RepresentativeMass - f.RepresentativeMass);
                        if ((Math.Abs(massDiff - 1.0) < massTol || Math.Abs(massDiff - 2.0) < massTol) && SimilarScore(bestFeature, f))
                        {
                            outFeatures.Add(f);
                            continue;
                        }
                        
                        tempList.Add(f);
                    }
                }

                bestFeature.InActivateSignificantPeaks();
                foreach (var f in tempList)
                {
                    f.UpdateScores(_spectrums);
                    if (f.GoodEnough) featureSet.Add(f);
                }
            }

            return outFeatures;
        }
        
        public IList<SortedSet<Ms1FeatureCluster>> GetAllConnectedFeatures()
        {
            //Console.WriteLine("Generating feature graphs from {0} features", NumberOfFeatures);
            //var stopwatch = Stopwatch.StartNew();
            var ret = new List<SortedSet<Ms1FeatureCluster>>();
            var startIndex = 0;
            while (true)
            {
                var featureSet = GetConnectedFeatures(ref startIndex);
                if (featureSet.Count < 1) break;
                ret.Add(featureSet);
            }
            //stopwatch.Stop(); stopwatch.Reset();
            //var elapsed = (stopwatch.ElapsedMilliseconds) / 1000.0d;
            //Console.WriteLine("# of generated connected components = {0};  Elapsed Time = {1:0.000} sec", ret.Count, elapsed);
            return ret;
        }

        private SortedSet<Ms1FeatureCluster> GetConnectedFeatures(ref int startIndex)
        {
            var neighbors = new Queue<Ms1FeatureCluster>();
            for (var i = startIndex; i < _featureList.Count; i++)
            {
                if (_featureList[i].Flag == 1) continue;
                neighbors.Enqueue(_featureList[i]);
                startIndex = i + 1;
                break;
            }

            var featureSet = new SortedSet<Ms1FeatureCluster>(new Ms1FeatureComparer());
            while (true)
            {
                if (neighbors.Count < 1) break;
                var feature = neighbors.Dequeue();
                if (feature.Flag == 1) continue;

                feature.Flag = 1;
                featureSet.Add(feature);

                foreach (var f in feature.OverlappedFeatures)
                {
                    if (f.Flag == 1) continue;
                    neighbors.Enqueue(f);
                }
            }

            return featureSet;
        }

        public IEnumerable<Ms1FeatureCluster> GetAllFeatureClusters()
        {
            return _featureList; 
        }
        
        private readonly List<Ms1FeatureCluster> _featureList;
        private readonly List<Ms1Spectrum> _spectrums;
        public int NumberOfFeatures { get { return _featureList.Count;  } }
    }
}
