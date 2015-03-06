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

        public void Add(Ms1FeatureCluster newFeature)
        {
            foreach (var peak in newFeature.GetMajorPeaks())
                peak.TagMajorPeakOf(newFeature);

            foreach (var peak in newFeature.GetMajorPeaks())
                peak.TagMinorPeakOf(newFeature);

            _featureList.Add(newFeature);
        }

        public void Add(IEnumerable<Ms1FeatureCluster> features)
        {
            foreach(var f in features) Add(f);
        }

        public IList<Ms1FeatureCluster> GetFilteredFeatures(IList<SortedSet<Ms1FeatureCluster>> connectedFeatureList)
        {
            var filteredFeatures = new List<Ms1FeatureCluster>();
            //var stopwatch = Stopwatch.StartNew();
            var i = 0;
            foreach (var featureSet in connectedFeatureList)
            {
                i++;
                var n1 = featureSet.Count;
                //Console.Write("Processing {0} connected component; # of features = {1}", i, n1);
                var newList = RemoveOverlappedFeatures(featureSet);
                var n2 = newList.Count;
                filteredFeatures.AddRange(newList);
                //Console.WriteLine("...reduced to {0} ", n2);
            }

            //stopwatch.Stop();
            //var elapsed = (stopwatch.ElapsedMilliseconds) / 1000.0d;
            //Console.WriteLine("# of filtered features = {0};  Elapsed Time = {1:0.000} sec", filteredFeatures.Count, elapsed);
            return filteredFeatures;
        }

        private IList<Ms1FeatureCluster> RemoveOverlappedFeatures(SortedSet<Ms1FeatureCluster> featureSet)
        {
            var outFeatures = new List<Ms1FeatureCluster>();

            while (true)
            {
                if (featureSet.Count < 1) break;

                var bestFeature = featureSet.First();
                featureSet.Remove(bestFeature);
                outFeatures.Add(bestFeature);
                
                var tempList = new List<Ms1FeatureCluster>();
                foreach (var f in bestFeature.OverlappedFeatures)
                {
                    if (featureSet.Remove(f)) tempList.Add(f);
                }

                bestFeature.InActivateSignificantPeaks();

                foreach (var f in tempList)
                {
                    if (f.RepresentativeScanNum == 5726 && f.MinCharge == 9)
                    {
                        var debug = 0;
                    }
                    
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
