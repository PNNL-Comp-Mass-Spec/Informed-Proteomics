using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsFeatureAlignment
    {
        public LcMsFeatureAlignment(List<string> featureFileList, List<string> rawFileList, double massTolerancePpm = 5)
        {
            _tolerance = new Tolerance(massTolerancePpm);

            FeatureFileList = featureFileList;
            RawFileList = rawFileList;

            _featureSetList = new Dictionary<int, List<LcMsFeature>>();
            _featureList = new List<LcMsFeature>();

            for (var i = 0; i < FeatureFileList.Count; i++)
            {
                var features = LoadProMexResult(RawFileList[i], FeatureFileList[i], i);
                features.Sort(new FeatureMassComparer());
                _featureSetList.Add(i, features);
                _featureList.AddRange(features);
            }
            _featureList.Sort(new FeatureMassComparer());
            
        }
        
        public static List<LcMsFeature> LoadProMexResult(string rawFilePath, string featureFilePath, int dataid = 0)
        {
            var featureList = new List<LcMsFeature>();
            var tsvReader = new TsvFileParser(featureFilePath);
            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);
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

                //var feature = new LcMsFeature(run, minChg, maxChg, minScan, maxScan, abu, mass, repChg, repMz, repScanNum)
                var feature = new LcMsFeature(repMass, repCharge, repMz, repScanNum, abundance, minCharge, maxCharge, minScan, maxScan, run)
                {
                    FeatureId = fid,
                    DataSetId = dataid,
                };

                featureList.Add(feature);
            }

            return featureList;
        }

        private List<LcMsFeature[]> _alignedFeatures = null;
        public List<LcMsFeature[]> GroupFeatures()
        {
            return _alignedFeatures ?? (_alignedFeatures = GroupFeatures(_featureList));
        }

        //public List<List<LcMsFeature>> GroupFeatures(List<LcMsFeature> features)
        public List<LcMsFeature[]> GroupFeatures(List<LcMsFeature> features)
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
                    if (Alignable(fi, fj))
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

        private bool Alignable(LcMsFeature f1, LcMsFeature f2)
        {
            if (f1.DataSetId == f2.DataSetId) return false;

            var massTol = Math.Min(_tolerance.GetToleranceAsTh(f1.Mass), _tolerance.GetToleranceAsTh(f2.Mass));

            if (Math.Abs(f1.Mass - f2.Mass) > massTol) return false;
            if (f1.CoElutedByNet(f2, 0.001)) return true;

            if (NetDiff(f1, f2) < TolNet) return true;

            return false;
        }

        private List<LcMsFeature> GetAlignedFeatures(ref HashSet<int> component, List<int>[] adjList)
        {
            var nDataSet = FeatureFileList.Count;
            
            // find seed node
            var minScoreNode = 0;
            var minScore = 0d;
            var score = new double[nDataSet];
            var linkedIdx = new int[nDataSet];

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

                //foreach (var idx2 in adjList[idx1])
                foreach(var idx2 in component)
                {
                    //if (!component.Contains(idx2)) continue;
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
            var nodeSet = new HashSet<int>();
            nodeSet.Add(minScoreNode);
            chkDataset[_featureList[minScoreNode].DataSetId] = true;

            // clustering
            while (true)
            {
                minScore = 0d;
                minScoreNode = -1;
                //minScoreNodeIndex = -1;

                foreach (var idx1 in nodeSet)
                {
                    var f1 = _featureList[idx1];
                    //foreach (var idx2 in adjList[idx1])
                    foreach(var idx2 in component)
                    {
                        //if (!component.Contains(idx2)) continue;
                        if (idx1 == idx2 || !adjList[idx1].Contains(idx2)) continue;
                    
                        var f2 = _featureList[idx2];
                        if (chkDataset[f2.DataSetId]) continue;
                        var netDiff = NetDiff(f1, f2);

                        if (minScoreNode < 0 || netDiff < minScore)
                        {
                            minScore = netDiff;
                            minScoreNode = idx2;
                            //minScoreNodeIndex = 
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


        public int CountDatasets { get { return RawFileList.Count;  } }

        public int CountAlignedFeatures { get { return (_alignedFeatures == null) ? 0 : _alignedFeatures.Count; } }

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
        
        public readonly List<string> FeatureFileList;
        public readonly List<string> RawFileList;

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
