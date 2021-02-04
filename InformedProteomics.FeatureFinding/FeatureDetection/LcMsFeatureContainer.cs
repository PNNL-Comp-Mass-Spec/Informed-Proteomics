using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.FeatureFinding.Clustering;
using InformedProteomics.FeatureFinding.Data;
using InformedProteomics.FeatureFinding.Scoring;
using InformedProteomics.FeatureFinding.Util;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.FeatureFinding.FeatureDetection
{
    public class LcMsFeatureContainer
    {
        public LcMsFeatureContainer(List<Ms1Spectrum> ms1Spectra, LcMsFeatureLikelihood scorer, INodeComparer<LcMsPeakCluster> mergeComparer)
        {
            _featureList = new List<LcMsPeakCluster>();
            _spectra = ms1Spectra;
            _scorer = scorer;
            _mergeComparer = mergeComparer;
        }

        public bool Add(LcMsPeakCluster newFeature)
        {
            if (newFeature.Score < _scorer.ScoreThreshold)
            {
                return false;
            }

            if (!newFeature.GoodEnough)
            {
                return false;
            }

            for (var i = _featureList.Count - 1; i >= 0; i--)
            {
                var massDiff = Math.Abs(_featureList[i].RepresentativeMass - newFeature.RepresentativeMass);
                if (massDiff > 1.0d)
                {
                    break;
                }

                if (massDiff < 1e-4)
                {
                    var coeLen = _featureList[i].CoElutionLength(newFeature);
                    if (coeLen > _featureList[i].ElutionLength * 0.7 || coeLen > newFeature.ElutionLength * 0.7)
                    {
                        return false;
                    }
                }
            }
            /*
            foreach (var peak in newFeature.GetMajorPeaks())
            {
                peak.TagMajorPeakOf(newFeature);
            }

            foreach (var peak in newFeature.GetMinorPeaks())
            {
                peak.TagMinorPeakOf(newFeature);
            }
            */
            _featureList.Add(newFeature);

            return true;
        }

        public void Add(IEnumerable<LcMsPeakCluster> features)
        {
            foreach(var f in features)
            {
                Add(f);
            }
        }

        public IEnumerable<LcMsPeakCluster> GetFilteredFeatures(LcMsPeakMatrix featureFinder)
        {
            var mergedFeatures = MergeFeatures(featureFinder, _featureList);
            SetFeatureList(mergedFeatures);
            var connectedFeatures = GetAllConnectedFeatures();
            var filteredFeatures = GetFilteredFeatures(connectedFeatures);
            return filteredFeatures.OrderBy(f => f.RepresentativeMass);

            /*
            var connectedFeatures = GetAllConnectedFeatures();
            var filteredFeatures = GetFilteredFeatures(connectedFeatures);
            var mergedFeatures = MergeFeatures(featureFinder, filteredFeatures);
            return mergedFeatures.OrderBy(f => f.RepresentativeMass);
            */
        }

        private void SetFeatureList(IEnumerable<LcMsPeakCluster> features)
        {
            _featureList.Clear();
            _featureList.AddRange(features);

            foreach (var newFeature in _featureList)
            {
                foreach (var peak in newFeature.GetMajorPeaks())
                {
                    peak.TagMajorPeakOf(newFeature);
                }
                foreach (var peak in newFeature.GetMinorPeaks())
                {
                    peak.TagMinorPeakOf(newFeature);
                }
            }
        }

        private List<LcMsPeakCluster> GetFilteredFeatures(IList<SortedSet<LcMsPeakCluster>> connectedFeatureList)
        {
            var filteredFeatures = new List<LcMsPeakCluster>();
            foreach (var featureSet in connectedFeatureList)
            {
                var newList = RemoveOverlappedFeatures(featureSet);
                filteredFeatures.AddRange(newList);
            }
            return filteredFeatures;
        }

        private IList<SortedSet<LcMsPeakCluster>> GetAllConnectedFeatures()
        {
            var ret = new List<SortedSet<LcMsPeakCluster>>();
            var startIndex = 0;
            while (true)
            {
                var featureSet = GetConnectedFeatures(ref startIndex);
                if (featureSet.Count < 1)
                {
                    break;
                }

                ret.Add(featureSet);
            }
            //stopwatch.Stop(); stopwatch.Reset();
            //var elapsed = (stopwatch.ElapsedMilliseconds) / 1000.0d;
            //Console.WriteLine("# of generated connected components = {0};  Elapsed Time = {1:0.000} sec", ret.Count, elapsed);
            return ret;
        }

        private List<LcMsPeakCluster> MergeFeatures(LcMsPeakMatrix featureFinder, List<LcMsPeakCluster> features)
        {
            //foreach (var f in _featureList) f.ActivateAllPeaks();
            var featureSet = new NodeSet<LcMsPeakCluster>();
            featureSet.AddRange(features);

            var connectedFeatureSet = featureSet.ConnectedComponents(_mergeComparer);
            var mergedFeatures = new List<LcMsPeakCluster>();

            foreach (var fSet in connectedFeatureSet)
            {
                if (fSet.Count == 1)
                {
                    mergedFeatures.Add(fSet[0]);
                }
                else
                {
                    var maxScan = fSet.Max(f => f.MaxScanNum);
                    var minScan = fSet.Min(f => f.MinScanNum);
                    var maxCharge = fSet.Max(f => f.MaxCharge);
                    var minCharge = fSet.Min(f => f.MinCharge);
                    var maxScore = double.MinValue;//fSet.Max(f => f.Score);
                    LcMsPeakCluster maxScoredClusterOriginal = null;
                    LcMsPeakCluster maxScoredCluster = null;
                    foreach (var f in fSet)
                    {
                        var newFeature = featureFinder.GetLcMsPeakCluster(f.RepresentativeMass, minCharge, maxCharge, minScan, maxScan);
                        if (newFeature != null && (maxScoredCluster == null || newFeature.Score > maxScoredCluster.Score))
                        {
                            maxScoredCluster = newFeature;
                        }

                        if (f.Score > maxScore)
                        {
                            maxScoredClusterOriginal = f;
                            maxScore = f.Score;
                        }
                    }
                    var feature = featureFinder.GetLcMsPeakCluster(fSet.Select(f => f.Mass).Mean(), minCharge, maxCharge, minScan, maxScan);
                    if (feature != null && (maxScoredCluster == null || feature.Score > maxScoredCluster.Score))
                    {
                        maxScoredCluster = feature;
                    }
                    //Console.WriteLine("------------- Merge -----------------");
                    //foreach (var f in fSet) Console.WriteLine("*\t{0}\t{1}\t{2}\t{3}", f.RepresentativeMass, f.MinScanNum, f.MaxScanNum, f.Score);
                    //Console.WriteLine("**\t{0}\t{1}\t{2}\t{3}", maxScoredCluster.RepresentativeMass, maxScoredCluster.MinScanNum, maxScoredCluster.MaxScanNum, maxScoredCluster.Score);
                    if (maxScoredCluster == null)
                    {
                        maxScoredCluster = maxScoredClusterOriginal;
                    }

                    if (maxScoredCluster != null && maxScoredCluster.Score < maxScore)
                    {
                        maxScoredCluster.Score = maxScore;
                    }

                    mergedFeatures.Add(maxScoredCluster);
                }
                //if (selectedFeature != null) postFilteredSet.Add(selectedFeature);
            }
            //return postFilteredSet.OrderBy(f => f.RepresentativeMass);

            return mergedFeatures;
        }

        private bool SimilarScore(LcMsPeakCluster f1, LcMsPeakCluster f2)
        {
            /*var maxScore = Math.Max(f1.Score, f2.Score);
            var minScore = Math.Min(f1.Score, f2.Score);
            if (minScore > 0 && maxScore > minScore*5) return false;*/

            if (f1.Score >= _scorer.ScoreThreshold && f1.GoodEnough
             && f2.Score >= _scorer.ScoreThreshold && f2.GoodEnough)
            {
                return true;
            }

            return false;
        }

        private IEnumerable<LcMsPeakCluster> RemoveOverlappedFeatures(ISet<LcMsPeakCluster> featureSet)
        {
            var outFeatures = new List<LcMsPeakCluster>();
            var tol = new Tolerance(5);
            while (true)
            {
                if (featureSet.Count < 1)
                {
                    break;
                }

                var bestFeature = featureSet.First();
                featureSet.Remove(bestFeature);
                outFeatures.Add(bestFeature);
                var massTol = tol.GetToleranceAsMz(bestFeature.RepresentativeMass);

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
                    if (f.Score > _scorer.ScoreThreshold && f.GoodEnough)
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

        private SortedSet<LcMsPeakCluster> GetConnectedFeatures(ref int startIndex)
        {
            var neighbors = new Queue<LcMsPeakCluster>();
            for (var i = startIndex; i < _featureList.Count; i++)
            {
                if (_featureList[i].Flag == 1)
                {
                    continue;
                }

                neighbors.Enqueue(_featureList[i]);
                startIndex = i + 1;
                break;
            }

            var featureSet = new SortedSet<LcMsPeakCluster>(new LcMsFeatureScoreComparer());
            while (true)
            {
                if (neighbors.Count < 1)
                {
                    break;
                }

                var feature = neighbors.Dequeue();
                if (feature.Flag == 1)
                {
                    continue;
                }

                feature.Flag = 1;
                featureSet.Add(feature);

                foreach (var f in feature.OverlappedFeatures)
                {
                    if (f.Flag == 1)
                    {
                        continue;
                    }

                    neighbors.Enqueue(f);
                }
            }

            return featureSet;
        }
        public int NumberOfFeatures => _featureList.Count;

        private readonly LcMsFeatureLikelihood _scorer;

        private readonly List<LcMsPeakCluster> _featureList;

        private readonly List<Ms1Spectrum> _spectra;
        private readonly INodeComparer<LcMsPeakCluster> _mergeComparer;

        private class LcMsFeatureScoreComparer : IComparer<LcMsPeakCluster>
        {
            public int Compare(LcMsPeakCluster x, LcMsPeakCluster y)
            {
                if (x == null)
                {
                    return y == null ? 0 : -1;
                }

                if (y == null)
                {
                    return 1;
                }

                return y.Score.CompareTo(x.Score);
            }
        }
    }
}
