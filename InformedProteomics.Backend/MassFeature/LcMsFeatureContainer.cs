using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.Backend.MassFeature
{
    
    public class LcMsFeatureContainer
    {
        internal class PostFeatureComparer : INodeComparer<LcMsPeakCluster>
        {
            private readonly Tolerance _tolerance = new Tolerance(10);
            public bool SameCluster(LcMsPeakCluster f1, LcMsPeakCluster f2)
            {
                var massTh = _tolerance.GetToleranceAsTh(f1.RepresentativeMass);
                var massDiff = Math.Abs(f1.RepresentativeMass - f2.RepresentativeMass);
                if (massDiff > massTh) return false;
                var coeLen = f1.CoElutionLength(f2);
                
                if (coeLen > f2.ElutionLength*0.5 || coeLen > f1.ElutionLength*0.5) return true; // significantly overlapped features

                if (f1.CoElutedByNet(f2, 0.003)) // close in elution time
                {
                    
                    // two differnt hills
                    if (Math.Abs(f1.ApexElutionTime - f2.ApexElutionTime) > 1.0d &&
                        f1.ApexIntensity/f1.MedianIntensity > 2 &&
                        f2.ApexIntensity/f2.MedianIntensity > 2) return false;
                    
                    // otherwise, they are fragmentized features, which
                    return true;
                }
                /*
                if (f1.MaxElutionTime < f2.MinElutionTime && f2.MinElutionTime - f1.MaxElutionTime < 0.5)
                    return true;

                if (f2.MaxElutionTime < f1.MinElutionTime && f1.MinElutionTime - f2.MaxElutionTime < 0.5)
                    return true;
                */
                return false;
            }
        }


        public LcMsFeatureContainer(List<Ms1Spectrum> ms1Spectra, LcMsFeatureLikelihood scorer)
        {
            _featureList = new List<LcMsPeakCluster>();
            _spectra = ms1Spectra;
            _scorer = scorer;
        }
        
        public bool Add(LcMsPeakCluster newFeature)
        {
            if (newFeature.Score < _scorer.ScoreThreshold) return false;
            if (!newFeature.GoodEnougth) return false;
            //var tolerance = new Tolerance(4);
            //var binNum = _comparer.GetBinNumber(newFeature.RepresentativeMass);
            //var massTol = tolerance.GetToleranceAsTh(newFeature.RepresentativeMass);
            for (var i = _featureList.Count - 1; i >= 0; i--)
            {
                var massDiff = Math.Abs(_featureList[i].RepresentativeMass - newFeature.RepresentativeMass);
                if (massDiff > 1.0d) break;
                if (massDiff < 1e-4)
                {
                    var coeLen = _featureList[i].CoElutionLength(newFeature);
                    if (coeLen > _featureList[i].ElutionLength*0.8 || coeLen > newFeature.ElutionLength*0.8) return false;
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
                filteredFeatures.AddRange(newList);
            }

            return filteredFeatures;
            //return PostFiltering(filteredFeatures, new PostFeatureComparer());
        }

        private IEnumerable<LcMsPeakCluster> MergeFeatures(List<LcMsPeakCluster> features, INodeComparer<LcMsPeakCluster> featureComparer, LcMsPeakMatrix featureFinder)
        {
            var featureSet = new NodeSet<LcMsPeakCluster>();
            featureSet.AddRange(features);
            var connectedFeatureSet = featureSet.ConnnectedComponents(featureComparer);
            //return connectedFeatureSet;
            
            var postFilteredSet = new List<LcMsPeakCluster>();
            foreach (var fSet in connectedFeatureSet)
            {
                if (fSet.Count == 1)
                {
                    postFilteredSet.Add(fSet[0]);
                }
                else
                {
                    var maxScan = fSet.Max(f => f.MaxScanNum);
                    var minScan = fSet.Min(f => f.MinScanNum);
                    var maxCharge = fSet.Max(f => f.MaxCharge);
                    var minCharge = fSet.Min(f => f.MinCharge);
                    
                    var mass = fSet.Select(f => f.Mass).Median();
                    var score = fSet.Max(f => f.Score);
                    var newFeature = featureFinder.GetLcMsPeakCluster(mass, minCharge, maxCharge, minScan, maxScan);
                    if (newFeature.Score < score) newFeature.Score = score;
                    postFilteredSet.Add(newFeature);
                }
                //if (selectedFeature != null) postFilteredSet.Add(selectedFeature);
            }
            return postFilteredSet.OrderBy(f => f.RepresentativeMass);
        }

        /*
        private IEnumerable<LcMsPeakCluster> PostFiltering(List<LcMsPeakCluster> features)
        {
            var orderedFeatures = features.OrderBy(f => f.RepresentativeMass).ToList();
            var postFiltered = new List<LcMsPeakCluster>();
            var tolerance = new Tolerance(16);

            foreach (var f in orderedFeatures) f.Flag = 0;

            for (var i = 0; i < orderedFeatures.Count; i++)
            {
                var f1 = orderedFeatures[i];

                if (f1.Flag == 2) continue;

                var massTh = tolerance.GetToleranceAsTh(f1.RepresentativeMass);

                //var isOkay = true;
                for (var j = i + 1; j < orderedFeatures.Count; j++)
                {
                    var f2 = orderedFeatures[j];

                    if (f2.Flag == 2) continue;

                    var massDiff = Math.Abs(f1.RepresentativeMass - f2.RepresentativeMass);
                    if (massDiff > massTh) break;

                    var coeLen = f1.CoElutionLength(f2);
                    if (coeLen > f2.ElutionLength*0.7 || coeLen > f1.ElutionLength*0.7)
                    {
                        if (f1.Score < f2.Score)
                        {
                            f1.Flag = 2;
                            break;
                        }
                        else
                        {
                            f1.Flag = 1;
                            f2.Flag = 2;
                        }
                    }
                }

                if (f1.Flag == 0) f1.Flag = 1;
            }
        }*/
        
        private bool SimilarScore(LcMsPeakCluster f1, LcMsPeakCluster f2)
        {
            /*var maxScore = Math.Max(f1.Score, f2.Score);
            var minScore = Math.Min(f1.Score, f2.Score);
            if (minScore > 0 && maxScore > minScore*5) return false;*/

            if (f1.Score >= _scorer.ScoreThreshold && f1.GoodEnougth
             && f2.Score >= _scorer.ScoreThreshold && f2.GoodEnougth) return true;

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
                    if (f.Score > _scorer.ScoreThreshold && f.GoodEnougth)
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
