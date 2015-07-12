using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsFeatureContainer
    {
        public LcMsFeatureContainer(List<Ms1Spectrum> ms1Spectra, LcMsFeatureLikelihood scorer)
        {
            _featureList = new List<LcMsPeakCluster>();
            _spectra = ms1Spectra;
            _scorer = scorer;
        }


        public const double ScoreThreshold = 0;

        public bool Add(LcMsPeakCluster newFeature)
        {
            if (newFeature.Score < ScoreThreshold) return false;
            if (!newFeature.GoodEnougth) return false;
            
            for (var i = _featureList.Count - 1; i >= 0; i--)
            {
                var massDiff = Math.Abs(_featureList[i].RepresentativeMass - newFeature.RepresentativeMass);
                if (massDiff < 1e-6)
                {
                    var coeLen = _featureList[i].CoElutionLength(newFeature);
                    if (coeLen / _featureList[i].ElutionLength > 0.6 && coeLen / newFeature.ElutionLength > 0.6) return false;
                }
            }

            foreach (var peak in newFeature.GetMajorPeaks())
            {
                peak.TagMajorPeakOf(newFeature);
            }

            foreach (var peak in newFeature.GetMinorPeaks())
            {
                peak.TagMinorPeakOf(newFeature);
            }
                
            
            _featureList.Add(newFeature);
            return true;
        }

        public void Add(IEnumerable<LcMsPeakCluster> features)
        {
            foreach(var f in features) Add(f);
        }

        public IEnumerable<LcMsPeakCluster> GetFilteredFeatures(IList<SortedSet<LcMsPeakCluster>> connectedFeatureList)
        {
            var filteredFeatures = new List<LcMsPeakCluster>();
            foreach (var featureSet in connectedFeatureList)
            {
                var newList = RemoveOverlappedFeatures(featureSet);
                //foreach (var f in newList) f.UpdateAbundance();
                filteredFeatures.AddRange(newList);
            }
            return filteredFeatures.OrderBy(f => f.RepresentativeMass);
        }

        
        private bool SimilarScore(LcMsPeakCluster f1, LcMsPeakCluster f2)
        {
            /*var maxScore = Math.Max(f1.Score, f2.Score);
            var minScore = Math.Min(f1.Score, f2.Score);
            if (minScore > 0 && maxScore > minScore*5) return false;*/
            
            if (f1.Score >= ScoreThreshold && f1.GoodEnougth
             && f2.Score >= ScoreThreshold && f2.GoodEnougth) return true;

            return false;
        }


        private IList<LcMsPeakCluster> RemoveOverlappedFeatures(SortedSet<LcMsPeakCluster> featureSet)
        {
            var outFeatures = new List<LcMsPeakCluster>();
            var tol = new Tolerance(5);
            while (true)
            {
                if (featureSet.Count < 1) break;

                var bestFeature = featureSet.First();
                featureSet.Remove(bestFeature);
                outFeatures.Add(bestFeature);
                var massTol = tol.GetToleranceAsTh(bestFeature.RepresentativeMass);

                var tempList = new List<LcMsPeakCluster>();
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

                bestFeature.InActivateMajorPeaks();
                foreach (var f in tempList)
                {
                    f.UpdateScore(_spectra);
                    f.Score = _scorer.GetScore(f);
                    if (f.Score > ScoreThreshold && f.GoodEnougth)
                    {
                        featureSet.Add(f);
                    }
                    else
                    {
                        //Console.WriteLine("{0}\t{1}\t{2} killed by {3}\t{4}\t{5}", f.Mass, f.MinScanNum, f.MaxScanNum, bestFeature.Mass, bestFeature.MinScanNum, bestFeature.MaxScanNum);
                    }
                }
            }

            return outFeatures;
        }

        public IList<SortedSet<LcMsPeakCluster>> GetAllConnectedFeatures()
        {
            var ret = new List<SortedSet<LcMsPeakCluster>>();
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

        private SortedSet<LcMsPeakCluster> GetConnectedFeatures(ref int startIndex)
        {
            var neighbors = new Queue<LcMsPeakCluster>();
            for (var i = startIndex; i < _featureList.Count; i++)
            {
                if (_featureList[i].Flag == 1) continue;
                neighbors.Enqueue(_featureList[i]);
                startIndex = i + 1;
                break;
            }

            var featureSet = new SortedSet<LcMsPeakCluster>(new LcMsFeatureComparer());
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

        public IEnumerable<LcMsPeakCluster> GetAllFeatureClusters()
        {
            return _featureList; 
        }

        private readonly LcMsFeatureLikelihood _scorer;
        private readonly List<LcMsPeakCluster> _featureList;
        private readonly List<Ms1Spectrum> _spectra;
        public int NumberOfFeatures { get { return _featureList.Count;  } }


        private class LcMsFeatureComparer : IComparer<LcMsPeakCluster>
        {
            public int Compare(LcMsPeakCluster x, LcMsPeakCluster y)
            {
                return (y.Score.CompareTo(x.Score));
            }
        }

    }
}
