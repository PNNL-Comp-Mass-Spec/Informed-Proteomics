using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsFeatureAlignment
    {
        public LcMsFeatureAlignment(ILcMsFeatureComparer featureCompare)
        {
            _featureSetList = new Dictionary<int, List<LcMsFeature>>();
            _featureList = new List<LcMsFeature>();
            RawFileList = new List<string>();
            _featureCompare = featureCompare;
        }

        private readonly ILcMsFeatureComparer _featureCompare;
        public LcMsFeatureAlignment(IList<string> featureFileList, IList<string> rawFileList, ILcMsFeatureComparer featureCompare)
        {
            RawFileList = rawFileList;
            _featureSetList = new Dictionary<int, List<LcMsFeature>>();
            _featureList = new List<LcMsFeature>();

            for (var i = 0; i < featureFileList.Count; i++)
            {
                var features = LoadProMexResult(featureFileList[i], RawFileList[i], i);
                features.Sort(new FeatureMassComparer());
                _featureSetList.Add(i, features);
                _featureList.AddRange(features);
            }
            _featureList.Sort(new FeatureMassComparer());

            _featureCompare = featureCompare;
        }
        
        public void AddDataSet(int dataId, List<LcMsFeature> features, string rawFilePath)
        {
            RawFileList.Add(rawFilePath);
            features.Sort(new FeatureMassComparer());
            _featureSetList.Add(dataId, features);
            _featureList.AddRange(features);
            _featureList.Sort(new FeatureMassComparer());        
        }
     
        public static List<LcMsFeature> LoadProMexResult(string featureFilePath, string rawFilePath = null, int dataid = 0)
        {
            var featureList = new List<LcMsFeature>();
            var tsvReader = new TsvFileParser(featureFilePath);
            var run = (rawFilePath == null || !File.Exists(rawFilePath)) ? null : PbfLcMsRun.GetLcMsRun(rawFilePath);
            var featureIds = tsvReader.GetData("FeatureID");

            var minScans = tsvReader.GetData("MinScan");
            var maxScans = tsvReader.GetData("MaxScan");
            var abu = tsvReader.GetData("Abundance");

            var minCharges = tsvReader.GetData("MinCharge");
            var maxCharges = tsvReader.GetData("MaxCharge");
            var monoMass = tsvReader.GetData("MonoMass");

            var repCharges = tsvReader.GetData("RepCharge");
            var repScans = tsvReader.GetData("RepScan");
            var repMzs = tsvReader.GetData("RepMz");

            var scores = tsvReader.GetData("LikelihoodRatio");
            var minElutionTime = tsvReader.GetData("MinElutionTime");
            var maxElutionTime = tsvReader.GetData("MaxElutionTime");

            for (var i = 0; i < tsvReader.NumData; i++)
            {
                var abundance = double.Parse(abu[i]);
                var repMass = double.Parse(monoMass[i]);
                var minCharge = int.Parse(minCharges[i]);
                var maxCharge = int.Parse(maxCharges[i]);
                var minScan = int.Parse(minScans[i]);
                var maxScan = int.Parse(maxScans[i]);
                var fid = int.Parse(featureIds[i]);
                var repCharge = int.Parse(repCharges[i]);
                var repMz = double.Parse(repMzs[i]);
                var repScanNum = int.Parse(repScans[i]);
                var score = double.Parse(scores[i]);
                var minEt = double.Parse(minElutionTime[i]);
                var maxEt = double.Parse(maxElutionTime[i]);
                
                var feature = new LcMsFeature(repMass, repCharge, repMz, repScanNum, abundance, minCharge, maxCharge, minScan, maxScan, minEt, maxEt, run)
                {
                    FeatureId = fid,
                    DataSetId = dataid,
                    Score = score,
                };

                featureList.Add(feature);
            }

            return featureList;
        }

        public void AlignFeatures()
        {
            _alignedFeatures = GroupFeatures(_featureList);
        }

        public List<LcMsFeature[]> GetAlignedFeatures()
        {
            return _alignedFeatures;
        }
        
        public void RefineAbundance()
        {
            if (_alignedFeatures == null) return;

            for (var i = 0; i < CountDatasets; i++)
            {
                var run = PbfLcMsRun.GetLcMsRun(RawFileList[i]);
                var ms1ScanNums = run.GetMs1ScanVector();
                var featureFinder = new LcMsPeakMatrix(run, new LcMsFeatureLikelihood());

                for (var j = 0; j < CountAlignedFeatures; j++)
                {
                    //if (_alignedFeatures[j][i] != null) continue;

                    var mass = 0d;
                    var charge = 0;
                    var minScanNum = -1;
                    var maxScanNum = ms1ScanNums.Last();
                    var repFt = GetRepFeatureInfo(_alignedFeatures[j]);

                    if (_alignedFeatures[j][i] == null)
                    {
                        mass = repFt.Mass;
                        charge = repFt.Charge;
                        var minNet = repFt.MinNet;
                        var maxNet = repFt.MaxNet;

                        for (var k = 0; k < ms1ScanNums.Length; k++)
                        {
                            var net = run.GetElutionTime(ms1ScanNums[k]) / run.GetElutionTime(run.MaxLcScan);
                            if (net > minNet && minScanNum < 0)
                            {
                                minScanNum = (k == 0) ? ms1ScanNums[k] : ms1ScanNums[k - 1];
                            }

                            if (net > maxNet)
                            {
                                maxScanNum = ms1ScanNums[k];
                                break;
                            }
                        }
                        if (minScanNum < 0) minScanNum = 0;
                    }
                    else
                    {
                        mass = _alignedFeatures[j][i].Mass;
                        charge = _alignedFeatures[j][i].RepresentativeCharge;
                        minScanNum = _alignedFeatures[j][i].MinScanNum;
                        maxScanNum = _alignedFeatures[j][i].MaxScanNum;
                    }
                    
                    var feature = featureFinder.GetLcMsPeakCluster(mass, charge, minScanNum, maxScanNum);
                    if (feature == null || feature.Score < -10)
                    {
                        feature = featureFinder.CollectLcMsPeaksWithNoise(mass, charge, minScanNum, maxScanNum, repFt.MinCharge, repFt.MaxCharge);
                    }
                    
                    _alignedFeatures[j][i] = feature;    
                }
                Console.WriteLine("{0} has been processed...", RawFileList[i]);
            }
        }

        private class RepFeatureInfo
        {
            internal double Mass;
            internal int Charge;
            internal double MinNet;
            internal double MaxNet;
            internal int MinCharge;
            internal int MaxCharge;
        }

        private RepFeatureInfo GetRepFeatureInfo(LcMsFeature[] features)
        {
            //var minElutionTime = double.MaxValue;
            //var maxElutionTime = 0d;
            var minNet = 1.0d;
            var maxNet = 0d;
            var massList = new List<double>();
            var chargeList = new List<double>();

            var minCharge = int.MaxValue;
            var maxCharge = -1;
            
            foreach (var f in features)
            {
                if (f == null) continue;
                minNet = Math.Min(minNet, f.MinNet);
                maxNet = Math.Max(maxNet, f.MaxNet);
                minCharge = Math.Min(minCharge, f.MinCharge);
                maxCharge = Math.Max(maxCharge, f.MaxCharge);
                //minElutionTime = Math.Min(minElutionTime, f.MinElutionTime);
                //maxElutionTime = Math.Max(maxElutionTime, f.MaxElutionTime);
                massList.Add(f.Mass);
                chargeList.Add(f.RepresentativeCharge);
            }
            
            var mass = massList.Median();
            var charge = (int) chargeList.Median();

            var ret = new RepFeatureInfo()
            {
                Mass = mass,
                Charge = charge,
                MinNet = minNet,
                MaxNet = maxNet,
                MinCharge = minCharge,
                MaxCharge = maxCharge,
            };
            return ret;
        }


        private LcMsFeature GetRepresentativeFeature(IList<LcMsFeature> features)
        {
            foreach (var f in features)
            {
                if (f == null) continue;
                return f;
            }
            return null;
        }
        
        public int CountDatasets { get { return RawFileList.Count; } }

        public int CountAlignedFeatures { get { return (_alignedFeatures == null) ? 0 : _alignedFeatures.Count; } }

        public void TryFillMissingFeature(List<LcMsFeature[]> alignedFeatures)
        {
            var featureInfoList = new List<LcMsFeature>();
            foreach (var features in alignedFeatures)
            {
                var featureInfo = GetRepresentativeFeature(features);
                featureInfoList.Add(featureInfo);
            }
            
            for (var i = 0; i < alignedFeatures.Count; i++)
            {
                var features = alignedFeatures[i];
                for (var j = 0; j < features.Length; j++)
                {
                    if (features[j] == null)
                    {
                        LcMsFeature altFt = null;

                        for (var k = i - 1; altFt == null && k >= 0; k--)
                        {
                            if (Math.Abs(featureInfoList[k].Mass - featureInfoList[i].Mass) > 2.5) break;

                            if (alignedFeatures[k][j] != null && _featureCompare.Match(featureInfoList[i], alignedFeatures[k][j]))
                                altFt = alignedFeatures[k][j];
                        }

                        for (var k = i + 1; altFt == null && k < alignedFeatures.Count; k++)
                        {
                            if (Math.Abs(featureInfoList[k].Mass - featureInfoList[i].Mass) > 2.5) break;
                            if (alignedFeatures[k][j] != null && _featureCompare.Match(featureInfoList[i], alignedFeatures[k][j]))
                                altFt = alignedFeatures[k][j];
                        }

                        if (altFt != null)
                        {
                            var newOnew = new LcMsFeature(featureInfoList[i].Mass, altFt.RepresentativeCharge,
                                altFt.RepresentativeMz, altFt.RepresentativeScanNum,
                                altFt.Abundance, altFt.MinCharge, altFt.MaxCharge, altFt.MinScanNum, altFt.MaxScanNum,
                                altFt.MinElutionTime, altFt.MaxElutionTime,
                                altFt.Run);

                            alignedFeatures[i][j] = newOnew;
                        }
                    }
                }
            }
        }

        private List<LcMsFeature[]> _alignedFeatures = null;
        private List<LcMsFeature[]> GroupFeatures(List<LcMsFeature> features)
        {
            var adjList = new List<int>[features.Count];

            for (var i = 0; i < features.Count; i++) adjList[i] = new List<int>();
            for (var i = 0; i < features.Count; i++)
            {
                var fi = features[i];

                for (var j = i + 1; j < features.Count; j++)
                {
                    var fj = features[j];
                    if (fj.Mass - fi.Mass > 1.5) break;
                    if (_featureCompare.Match(fi, fj))
                    {
                        adjList[i].Add(j);
                        adjList[j].Add(i);
                    }
                }
            }

            var ret = new List<LcMsFeature[]>();
            var components = GetConnectedComponents(adjList);
            for (var i = 0; i < components.Count; i++)
            {
                var component = new HashSet<int>(components[i]);
                while (component.Count > 0)
                {
                    var featureGroup = new LcMsFeature[CountDatasets];
                    var featureSet = GetAlignedFeatures(ref component, adjList);
                    foreach (var f in featureSet) featureGroup[f.DataSetId] = f;
                    ret.Add(featureGroup);
                }
            }

            return ret;
        }

        private List<LcMsFeature> GetAlignedFeatures(ref HashSet<int> component, List<int>[] adjList)
        {
            var nDataSet = RawFileList.Count;
            
            // find seed node
            var minScoreNode = 0;
            var minScore = 0d;
            var score = new double[nDataSet];
            var linkedIdx = new int[nDataSet];

            // find a cluster center that minimizing the sum of distances to other members
            foreach(var idx1 in component)
            {
                var f1 = _featureList[idx1];

                for (var k = 0; k < nDataSet; k++)
                {
                    score[k] = TolNet;
                    linkedIdx[k] = -1;
                }

                linkedIdx[f1.DataSetId] = idx1;
                score[f1.DataSetId] = 0;

                foreach(var idx2 in component)
                {
                    if (idx1 == idx2 || !adjList[idx1].Contains(idx2)) continue;

                    var f2 = _featureList[idx2];
                    var netDiff = NetDiff(f1, f2);

                    if (netDiff < score[f2.DataSetId])
                    {
                        score[f2.DataSetId] = netDiff;
                        linkedIdx[f2.DataSetId] = idx2;
                    }
                }

                var s = score.Sum();

                if (!(minScore > 0) || s < minScore)
                {
                    minScore = s;
                    minScoreNode = idx1;
                }
            }

            var chkDataset = new bool[nDataSet];
            //var nodeSet = new HashSet<int>();
            var nodeSet = new List<int>();
            nodeSet.Add(minScoreNode);
            chkDataset[_featureList[minScoreNode].DataSetId] = true;

            // clustering
            while (true)
            {
                minScore = 0d;
                minScoreNode = -1;

                foreach (var idx1 in nodeSet)
                {
                    var f1 = _featureList[idx1];
                    
                    foreach(var idx2 in component)
                    {
                        if (idx1 == idx2 || !adjList[idx1].Contains(idx2)) continue;
                        var f2 = _featureList[idx2];
                        if (chkDataset[f2.DataSetId]) continue;
                        var netDiff = NetDiff(f1, f2);

                        if (minScoreNode < 0 || netDiff < minScore)
                        {
                            minScore = netDiff;
                            minScoreNode = idx2;
                        }
                    }
                }
                if (minScoreNode < 0) break;

                // pick closest node to current cluster
                nodeSet.Add(minScoreNode);
                
                chkDataset[_featureList[minScoreNode].DataSetId] = true;

                if (nodeSet.Count == component.Count) break;
            }

            var alignedFeatures = new List<LcMsFeature>(nodeSet.Select(idx => _featureList[idx]));
            component.ExceptWith(nodeSet);

            return alignedFeatures;
        }


        private List<List<int>> GetConnectedComponents(List<int>[] adjList)
        {
            var nNodes = adjList.Length;
            //var groupedFeaters = new List<AlignedMs1FeatureSet>();
            var nodeSetList = new List<List<int>>();
            var visited = new bool[nNodes];
            for (var i = 0; i < nNodes; i++)
            {
                if (visited[i]) continue;

                //var featureSet = new AlignedMs1FeatureSet();
                var nodeSet = new List<int>();
                var neighbors = new Queue<int>();
                neighbors.Enqueue(i);
                while (true)
                {
                    if (neighbors.Count < 1) break;
                    var j = neighbors.Dequeue();
                    if (visited[j]) continue;
                    visited[j] = true;

                    //featureSet.Add(features[j]);
                    nodeSet.Add(j);
                    foreach (var k in adjList[j])
                    {
                        if (visited[k]) continue;
                        neighbors.Enqueue(k);
                    }
                }
                //groupedFeaters.Add(featureSet);
                nodeSetList.Add(nodeSet);

            }
            return nodeSetList;
        }


        private double NetDiff(LcMsFeature f1, LcMsFeature f2)
        {
            var n1 = Math.Abs(f1.Net - f2.Net);
            var n2 = Math.Abs(f1.MinNet - f2.MinNet);
            var n3 = Math.Abs(f1.MaxNet - f2.MaxNet);
            return Math.Min(Math.Min(n1, n2), n3);
        }

        private const double TolNet = 0.003;
        private readonly Dictionary<int, List<LcMsFeature>> _featureSetList;
        private readonly List<LcMsFeature> _featureList;
        public readonly IList<string> RawFileList;
        private readonly Tolerance _tolerance;

        private class FeatureMassComparer : IComparer<LcMsFeature>
        {
            public int Compare(LcMsFeature x, LcMsFeature y)
            {
                return (x.Mass.CompareTo(y.Mass));
            }
        }


    }
}
