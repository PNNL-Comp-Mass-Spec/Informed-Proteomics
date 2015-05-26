using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Quantification;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.Backend.Data.Spectrometry
{
   
    public class Ms1FeatureAlign
    {
        public Ms1FeatureAlign(List<string> featureFileList, List<string> rawFileList)
        {
            _tolerance = new Tolerance(5);

            _featureFileList = featureFileList;
            _rawFileList = rawFileList;

            _featureSetList = new Dictionary<int, List<Ms1Feature>>();
            _featureList = new List<Ms1Feature>();

            for (var i = 0; i < _featureFileList.Count; i++)
            {
                var features = Ms1FeatureResult.LoadProMexResult(_rawFileList[i], _featureFileList[i], i);
                features.Sort(new Ms1FeatureMassComparer());
                _featureSetList.Add(i, features);
                _featureList.AddRange(features);
            }
            _featureList.Sort(new Ms1FeatureMassComparer());
        }

        
        
        private const double TolNet = 0.003;

        public bool Alignable(Ms1Feature f1, Ms1Feature f2)
        {
            
            if (f1.DataSetId == f2.DataSetId) return false;

            var massTol = Math.Min(_tolerance.GetToleranceAsTh(f1.Mass), _tolerance.GetToleranceAsTh(f2.Mass));

            if (Math.Abs(f1.Mass - f2.Mass) > massTol) return false;
            if (f1.CoEluted(f2, 0.001)) return true;

            if (NetDiff(f1, f2) < TolNet) return true;

            return false;
        }

        private double NetDiff(Ms1Feature f1, Ms1Feature f2)
        {
            var n1 = Math.Abs(f1.Net - f2.Net);
            var n2 = Math.Abs(f1.MinNet - f2.MinNet);
            var n3 = Math.Abs(f1.MaxNet - f2.MaxNet);
            return Math.Min(Math.Min(n1, n2), n3);
        }

        public List<List<Ms1Feature>> GroupFeatures()
        {
            return GroupFeatures(_featureList);
        }

        public List<List<Ms1Feature>> GroupFeatures(List<Ms1Feature> features)
        {
            var adjList = new List<int>[features.Count];

            for (var i = 0; i < features.Count; i++) adjList[i] = new List<int>();

            for (var i = 0; i < features.Count; i++)
            {
                var fi = features[i];

                

                for (var j = i + 1; j < features.Count; j++)
                {
                    var fj = features[j];
                    /*
                    if (fi.DataSetId == 12-1 && fj.DataSetId == 11-1) 
                    {
                        if (fi.FeatureId == 18 && fj.FeatureId == 13)
                        {
                            var a = 0;
                        }
                    }*/
                    if (fj.Mass - fi.Mass > 1.5) break;
                    if (Alignable(fi, fj))
                    {
                        adjList[i].Add(j);
                        adjList[j].Add(i);
                    }
                }
            }

            var ret = new List<List<Ms1Feature>>();
            
            var components = GetConnectedComponents(adjList);
            for (var i = 0; i < components.Count; i++)
            {
                var component = new HashSet<int>(components[i]);
                
                while (component.Count > 0)
                {
                    //var featureSet = GetAlignedFeatureSet(component, adjList);    
                    var featureSet = GetAlignedFeatures(ref component, adjList);    

                    ret.Add(featureSet);
                }
            }
            return ret;

            /*
            var groupedFeaters = new List<AlignedMs1FeatureSet>();
            var visited = new bool[features.Count];
            for (var i = 0; i < features.Count; i++)
            {
                if (visited[i]) continue;

                var featureSet = new AlignedMs1FeatureSet();
                var neighbors = new Queue<int>();
                neighbors.Enqueue(i);
                while (true)
                {
                    if (neighbors.Count < 1) break;
                    var j = neighbors.Dequeue();
                    if (visited[j]) continue;
                    visited[j] = true;

                    featureSet.Add(features[j]);
                    foreach (var k in adjList[j])
                    {
                        if (visited[k]) continue;
                        neighbors.Enqueue(k);
                    }
                }
                groupedFeaters.Add(featureSet);
            }
            return groupedFeaters;*/
        }


        private List<Ms1Feature> GetAlignedFeatures(ref HashSet<int> component, List<int>[] adjList)
        {
            var nDataSet = _featureFileList.Count;
            
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

            var alignedFeatures = new List<Ms1Feature>(nodeSet.Select(idx => _featureList[idx]));

            component.ExceptWith(nodeSet);

            return alignedFeatures;
        }

        private List<Ms1Feature> GetAlignedFeatureSet(List<int> component, List<int>[] adjList)
        {
            var nDataSet = _featureFileList.Count;
            var minScore = 0d;
            int[] minScoreNodes = null;

            for (var i = 0; i < component.Count; i++)
            {
                var idx1 = component[i];
                var f1 = _featureList[idx1];

                var score = new double[nDataSet];
                var linkedIdx = new int[nDataSet];
                
                for (var k = 0; k < nDataSet; k++)
                {
                    score[k] = TolNet;
                    linkedIdx[k] = -1;
                }
                
                linkedIdx[f1.DataSetId] = idx1;
                score[f1.DataSetId] = 0;

                foreach (var idx2 in adjList[idx1])
                {
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
                    minScoreNodes = linkedIdx;
                }
            }

            
            var ret = new List<Ms1Feature>();
            var temp = new HashSet<int>();
            var temp1 = new HashSet<int>(component);

            foreach (var idx in minScoreNodes)
            {
                if (idx >= 0)
                {
                    ret.Add(_featureList[idx]);
                    temp.Add(idx);
                }
            }

            temp1.ExceptWith(temp);

            component.Clear();
            component.AddRange(temp1);

            return ret;
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

        /*
        private string GetHashString(AlignedMs1FeatureSet set)
        {
            var featureIdList = new int[_featureFileList.Count];
            foreach (var f in set)
            {
                featureIdList[f.DataSetId] = f.FeatureId;
            }
            var ret = new StringBuilder();

            for(var i = 0; i < featureIdList.Length; i++)
            {
                if (i != 0) ret.Append("-");
                ret.Append(featureIdList[i]);
            }
            //return BitConverter.ToString(featureIdList);
            return ret.ToString();
        }

        private bool SameFeature(AlignedMs1FeatureSet fset, Ms1FeatureResult f2)
        {
            const double tolNet = 0.001;
            //if (f1.DataSetIdIndex == f2.DataSetIdIndex) return false;
            //if (!f1.MassMatch(f2)) return false;
            //if (f1.CoEluted(f2, tolNet)) return true;

            foreach (var f in fset)
            {
                if (f.DataSetId == f2.DataSetId) continue;
                if (f.MassMatch(f2) && f.CoEluted(f2, tolNet)) return true;
            }

            return false;
        }

        public List<AlignedMs1FeatureSet> GroupFeatures()
        {
            for (var i = 0; i < _featureFileList.Count; i++)
            {
                for (var j = i + 1; j < _featureFileList.Count; j++)
                {
                    MatchPairDataset(_featureSetList[i], _featureSetList[j]);
                }
            }

            var alignments = new HashSet<string>();
            var groupedFeaters = new List<AlignedMs1FeatureSet>();
            foreach (var feature in _featureList)
            {
                var alignmentId = GetHashString(feature.AlignedSet);
                if (alignments.Contains(alignmentId)) continue;
                alignments.Add(alignmentId);
                groupedFeaters.Add(feature.AlignedSet);
            }
            
            return groupedFeaters;
        }

        private void MatchPairDataset(List<Ms1FeatureResult> seti, List<Ms1FeatureResult> setj)
        {
            var tol = new Tolerance(10);

            for (var i = 0; i < seti.Count; i++)
            {
                var fi = seti[i];
                var k = setj.BinarySearch(fi);
                if (k < 0) k = ~k;
                var massTol = tol.GetToleranceAsTh(fi.Mass);

                var minNetDiff = 1.0;
                Ms1FeatureResult alignee = null;

                for (var j = k; j < setj.Count; j++)
                {
                    var fj = setj[j];

                    if (Math.Abs(fj.Mass - fi.Mass) > massTol) break;
                    var netDiff = Math.Abs(fi.Net - fj.Net);
                    var lenDiff = Math.Abs(fi.NetLength - fj.NetLength);

                    if (netDiff < 0.12 && netDiff < minNetDiff && lenDiff < 0.05)
                    {
                        minNetDiff = netDiff;
                        alignee = fj;
                    }
                }
                for (var j = k - 1; j >= 0; j--)
                {
                    var fj = setj[j];
                    if (Math.Abs(fj.Mass - fi.Mass) > massTol) break;
                    var netDiff = Math.Abs(fi.Net - fj.Net);
                    var lenDiff = Math.Abs(fi.NetLength - fj.NetLength);

                    if (netDiff < 0.12 && netDiff < minNetDiff && lenDiff < 0.05)
                    {
                        minNetDiff = netDiff;
                        alignee = fj;
                    }
                }

                if (alignee != null)
                {
                    fi.AlignedSet.Add(alignee);
                    alignee.AlignedSet.Add(fi);
                }
            }
            
        }*/


        /*
       public List<AlignedMs1FeatureSet> ClusterFeatures(out double[,] alignedAbundance)
       {
           //var groupedFeatres = GroupFeatures(_featureList);
           var groupedFeatres = GroupFeatures();
           alignedAbundance = new double[groupedFeatres.Count, _featureFileList.Count];
            
           var featureID = 0;
           foreach (var featureSet in groupedFeatres)
           {
               foreach (var feature in featureSet)
               {
                   alignedAbundance[featureID, feature.DataSetIdIndex] += feature.Abundance;
               }
               featureID++;
           }
           return groupedFeatres;
       }

       
        public void AlignNetAll(int refDatasetIndex = 0)
        {
            for (var i = 0; i < _featureSetList.Count; i++)
            {
                if (i == refDatasetIndex) continue;

                AlignNet(_featureSetList[refDatasetIndex], _featureSetList[i]);
            }
        }
        
        private void GetInitialAlignPoint(List<Ms1FeatureResult> referenceSet, List<Ms1FeatureResult> targetSet, out double[] x, out double[] y)
        {
            referenceSet.Sort();
            targetSet.Sort();
            var tempX = new List<double>();
            var tempY = new List<double>();

            foreach (var refFt in referenceSet)
            {
                var index = targetSet.BinarySearch(refFt);
                if (index < 0) index = ~index;

                // go up
                var minNetDiff = 1.0d;
                var minNetDiffIdx = -1;

                var i = index;
                while (true)
                {
                    if (i >= targetSet.Count) break;
                    if (Math.Abs(refFt.Mass - targetSet[i].Mass) > 1.5) break;

                    var lenDiff = Math.Abs(refFt.NetLength - targetSet[i].NetLength);
                    var lenTh = Math.Max(Math.Max(refFt.NetLength, targetSet[i].NetLength)*0.4, 0.003);
                    var netDiff = Math.Abs(refFt.Net - targetSet[i].Net);

                    if (netDiff < 0.01 && lenDiff < lenTh && refFt.MassMatch(targetSet[i], true))
                    {
                        if (netDiff < minNetDiff)
                        {
                            minNetDiff = netDiff;
                            minNetDiffIdx = i;
                        }
                    }

                    i++;
                }

                // go down
                i = index - 1;
                while (true)
                {
                    if (i < 0) break;
                    if (Math.Abs(refFt.Mass - targetSet[i].Mass) > 1.5) break;

                    var lenDiff = Math.Abs(refFt.NetLength - targetSet[i].NetLength);
                    var lenTh = Math.Max(Math.Max(refFt.NetLength, targetSet[i].NetLength) * 0.4, 0.003);
                    var netDiff = Math.Abs(refFt.Net - targetSet[i].Net);

                    if (netDiff < 0.01 && lenDiff < lenTh && refFt.MassMatch(targetSet[i], true))
                    {
                        if (netDiff < minNetDiff)
                        {
                            minNetDiff = netDiff;
                            minNetDiffIdx = i;
                        }
                    }
                    i--;
                }

                if (minNetDiffIdx >= 0)
                {
                    tempX.Add(targetSet[minNetDiffIdx].ElutionTime);
                    tempY.Add(refFt.ElutionTime);
                }
            }

            x = tempX.ToArray();
            y = tempY.ToArray();
        }

        private void AlignNet(List<Ms1Feature> referenceSet, List<Ms1Feature> targetSet)
        {
            var allset = new List<Ms1Feature>();

            allset.AddRange(referenceSet);
            allset.AddRange(targetSet);
            allset.Sort();

            double[] x;
            double[] y;
            GetInitialAlignPoint(referenceSet, targetSet, out x, out y);
            var param = Fit.Line(x, y);
            foreach (var f in targetSet) f.SetElutionTimeFitParam(param);
        }
        */


        private readonly Dictionary<int, List<Ms1Feature>> _featureSetList;
        private readonly List<Ms1Feature> _featureList;
        
        private List<string> _featureFileList;
        private List<string> _rawFileList;

        private readonly Tolerance _tolerance;

        private class Ms1FeatureMassComparer : IComparer<Ms1Feature>
        {
            public int Compare(Ms1Feature x, Ms1Feature y)
            {
                return (x.Mass.CompareTo(y.Mass));
            }
        }       
        /*
        public static List<TargetFeature> MathMs1Ms2Result(int dataid, string featureFilePath, string rawFilePath, string idFilePath)
        {
            var ms2List = new List<TargetFeature>();
            var featureList = new List<Ms1FeatureResult>();
            var tsvReader = new TsvFileParser(featureFilePath);
            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);

            var featureIds = tsvReader.GetData("FeatureID");
            var minScans = tsvReader.GetData("MinScan");
            var maxScans = tsvReader.GetData("MaxScan");
            var minCharges = tsvReader.GetData("MinCharge");
            var maxCharges = tsvReader.GetData("MaxCharge");
            var monoMass = tsvReader.GetData("MonoMass");
            var abundance = tsvReader.GetData("Abundance");

            for (var i = 0; i < tsvReader.NumData; i++)
            {
                var abu = double.Parse(abundance[i]);
                var mass = double.Parse(monoMass[i]);
                var minChg = int.Parse(minCharges[i]);
                var maxChg = int.Parse(maxCharges[i]);
                var minScan = int.Parse(minScans[i]);
                var maxScan = int.Parse(maxScans[i]);
                var fid = int.Parse(featureIds[i]);
                var feature = new Ms1FeatureResult(run, dataid, fid, mass, minChg, maxChg, minScan, maxScan, abu);

                featureList.Add(feature);
            }
            featureList.Sort();

            tsvReader = new TsvFileParser(idFilePath);

            var scans = tsvReader.GetData("Scan");
            var mass2 = tsvReader.GetData("Mass");
            var charges = tsvReader.GetData("Charge");
            var qv = tsvReader.GetData("QValue");
            var comps = tsvReader.GetData("Composition");
            var pnames = tsvReader.GetData("ProteinName");


            var tol = new Tolerance(10);
            var tol2 = new Tolerance(15);
            for (var i = 0; i < tsvReader.NumData; i++)
            {
                var qvalue = double.Parse(qv[i]);
                var mass = double.Parse(mass2[i]);
                var charge = int.Parse(charges[i]);
                var scan = int.Parse(scans[i]);
                if (qvalue > 0.01) break;

                var ms2 = new TargetFeature(mass, charge, scan);
                ms2.ProteinName = pnames[i];
                ms2.Composition = comps[i];

                var k = featureList.BinarySearch(new Ms1FeatureResult(run, 0, 0, mass, 2, 2, scan, scan, 0));
                if (k < 0) k = ~k;

                for (var j = k; j < featureList.Count; j++)
                {
                    if (Math.Abs(featureList[j].Mass - mass) > tol.GetToleranceAsTh(mass)) break;
                    if (featureList[k].MinScanNum <= scan && scan <= featureList[k].MaxScanNum &&
                        featureList[k].MinCharge <= charge && charge <= featureList[k].MaxCharge)
                    {
                        ms2.SetMs1Feature(featureList[j]);
                        break;
                    }
                }

                if (ms2.LinkedMs1Feature == null)
                {
                    for (var j = k - 1; j >= 0; j--)
                    {
                        if (Math.Abs(featureList[j].Mass - mass) > tol.GetToleranceAsTh(mass)) break;
                        if (featureList[k].MinScanNum <= scan && scan <= featureList[k].MaxScanNum &&
                            featureList[k].MinCharge <= charge && charge <= featureList[k].MaxCharge)
                        {
                            ms2.SetMs1Feature(featureList[j]);
                            break;
                        }
                    }
                }

                if (ms2.LinkedMs1Feature == null)
                {
                    for (var j = k; j < featureList.Count; j++)
                    {
                        if (Math.Abs(featureList[j].Mass - mass) > tol2.GetToleranceAsTh(mass)) break;
                        if (featureList[k].MinScanNum <= scan && scan <= featureList[k].MaxScanNum &&
                            featureList[k].MinCharge <= charge && charge <= featureList[k].MaxCharge)
                        {
                            ms2.SetMs1Feature(featureList[j]);
                            break;
                        }
                    }

                    if (ms2.LinkedMs1Feature == null)
                    {
                        for (var j = k - 1; j >= 0; j--)
                        {
                            if (Math.Abs(featureList[j].Mass - mass) > tol2.GetToleranceAsTh(mass)) break;
                            if (featureList[k].MinScanNum <= scan && scan <= featureList[k].MaxScanNum &&
                                featureList[k].MinCharge <= charge && charge <= featureList[k].MaxCharge)
                            {
                                ms2.SetMs1Feature(featureList[j]);
                                break;
                            }
                        }
                    }
                }

                if (ms2.LinkedMs1Feature != null)
                {
                    var exist = false;
                    foreach (var ms2Ret in ms2List)
                    {
                        if (ms2Ret.LinkedMs1Feature == ms2.LinkedMs1Feature)
                        {
                            exist = true;
                            break;
                        }
                    }
                    if (!exist) ms2List.Add(ms2);
                }
            }
            return ms2List;
        }

        public static List<List<TargetFeature>> GroupTargetFeatures(List<TargetFeature> features)
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

                    if (fi.ProteinName == fj.ProteinName && fi.Composition == fj.Composition)
                    {
                        adjList[i].Add(j);
                        adjList[j].Add(i);
                    }
                }
            }

            var groupedFeaters = new List<List<TargetFeature>>();
            var visited = new bool[features.Count];
            for (var i = 0; i < features.Count; i++)
            {
                if (visited[i]) continue;

                var featureSet = new List<TargetFeature>();
                var neighbors = new Queue<int>();
                neighbors.Enqueue(i);
                while (true)
                {
                    if (neighbors.Count < 1) break;
                    var j = neighbors.Dequeue();
                    if (visited[j]) continue;
                    visited[j] = true;

                    featureSet.Add(features[j]);
                    foreach (var k in adjList[j])
                    {
                        if (visited[k]) continue;
                        neighbors.Enqueue(k);
                    }
                }
                groupedFeaters.Add(featureSet);
            }
            return groupedFeaters;
        }
        */
    }
 
}
