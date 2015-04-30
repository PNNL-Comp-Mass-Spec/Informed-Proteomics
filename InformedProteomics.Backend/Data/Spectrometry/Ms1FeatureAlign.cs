using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /*
    public class AlignedMs1FeatureSet : List<Ms1Feature>
    {
        public AlignedMs1FeatureSet() : base()
        {
        }
        
        public double MaxNet
        {
            get { return this.Max(f => f.MaxNet); }
        }

        public double MinNet
        {
            get { return this.Min(f => f.MinNet); }
        }

        public IList<double> GetMasses()
        {
            //var ret = this.Select(f => f.Mass).ToList().Distinct().ToArray();
            var ret = new List<double>();
            foreach (var f in this) ret.AddRange(f.GetMasses());
            return ret;
        }

        public double RepresentativeMass
        {
            get
            {
                var ret = GetMasses();
                return ret.Median();
            }
        }

        public int[] GetFeatureId(int dataSetId)
        {
            return this.Where(f => f.DataSetId == dataSetId).Select(f => f.FeatureId).ToArray();
        }
    }
    
    public class Ms1FeatureAlign
    {
        public Ms1FeatureAlign(List<string> featureFileList, List<string> rawFileList)
        {
            _tolerance = new Tolerance(5);

            _featureFileList = featureFileList;
            _rawFileList = rawFileList;
            _featureSetList = new Dictionary<int, List<Ms1Feature>>();
            _featureList = new List<Ms1Feature>();
            _nFeatures = new int[featureFileList.Count];

            ElutionLengths = new List<Tuple<int, double>>();

            for (var i = 0; i < featureFileList.Count; i++)
            {
                var features = LoadFeatureData(i, featureFileList[i], rawFileList[i]);
                _nFeatures[i] = features.Count;

                // merge features in the same dataset
                var groupList = GroupFeatures(features);
                var featuresMerged = new List<Ms1Feature>();
                foreach (var featureGroup in groupList)
                {
                    var repFeature = featureGroup[0];
                    for (var j = 1; j < featureGroup.Count; j++) repFeature.Merge(featureGroup[j]);
                    featuresMerged.Add(repFeature);
                }

                _featureSetList.Add(i, featuresMerged);
                _featureList.AddRange(featuresMerged);
            }
            _featureList.Sort();
        }
        
        public List<AlignedMs1FeatureSet> ClusterFeatures(out double[,] alignedAbundance)
        {
            var groupedFeatres = GroupFeatures(_featureList);
            alignedAbundance = new double[groupedFeatres.Count, _featureFileList.Count];
            
            var featureID = 0;
            foreach (var featureSet in groupedFeatres)
            {
                foreach (var feature in featureSet)
                {
                    alignedAbundance[featureID, feature.DataSetId] += feature.Abundance;
                }
                featureID++;
            }
            return groupedFeatres;
        }

        public List<AlignedMs1FeatureSet> GroupFeatures(List<Ms1Feature> features)
        {
            var adjList = new List<int>[features.Count];

            for (var i = 0; i < features.Count; i++) adjList[i] = new List<int>();

            for (var i = 0; i < features.Count; i++)
            {
                for (var j = i + 1; j < features.Count; j++)
                {
                    if (features[j].Mass - features[i].Mass > 2.3) break;
                    
                    //if (SameFeature(features[i], features[j]))
                    if (features[i].Equals(features[j]))
                    {
                        adjList[i].Add(j);
                        adjList[j].Add(i);
                    }
                }
            }

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
            return groupedFeaters;
        }
        
        public void AlignNetAll()
        {
            for(var i = 1; i < _featureSetList.Count; i++)
            {
                AlignNet(_featureSetList[0], _featureSetList[i]);
            }
        }

        private void GetInitialAlignPoint(List<Ms1Feature> referenceSet, List<Ms1Feature> targetSet, out double[] x, out double[] y)
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

        public List<Tuple<int, double>> ElutionLengths;
        private List<Ms1Feature> LoadFeatureData(int dataid, string featureFilePath, string rawFilePath)
        {
            var featureList = new List<Ms1Feature>();
            var tsvReader = new TsvFileParser(featureFilePath);
            var run = PbfLcMsRun.GetLcMsRun(rawFilePath);

            var maxScanNum = run.MaxLcScan;
            //var minScanNum = run.MinLcScan;

            ElutionLengths.Add(new Tuple<int, double>(maxScanNum, run.GetElutionTime(maxScanNum)));

            //Console.WriteLine(featureFilePath);
            //Console.WriteLine(tsvReader.NumData);

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
                //if (!(abu > 0)) continue;
                var mass = double.Parse(monoMass[i]);
                var minChg = int.Parse(minCharges[i]);
                var maxChg = int.Parse(maxCharges[i]);
                var minScan = int.Parse(minScans[i]);
                var maxScan = int.Parse(maxScans[i]);
                //var minNet = (double)(minScan - minScanNum) / (double)(maxScanNum - minScanNum);
                //var maxNet = (double)(maxScan - minScanNum) / (double)(maxScanNum - minScanNum);
                var fid = int.Parse(featureIds[i]);
                var feature = new Ms1Feature(run, dataid, fid, mass, minChg, maxChg, minScan, maxScan, abu);
                //feature.SetNet(minNet, maxNet);
                featureList.Add(feature);
            }
            return featureList;
        }

        public int GetFeatureCount(int index)
        {
            return _nFeatures[index];
        }

        private readonly Dictionary<int, List<Ms1Feature>> _featureSetList;
        private readonly List<Ms1Feature> _featureList;
        private List<string> _featureFileList;
        private List<string> _rawFileList;
        private readonly Tolerance _tolerance;
        private readonly int[] _nFeatures;

        private class Ms1FeatureAbundanceComparer : IComparer<Ms1Feature>
        {
            public int Compare(Ms1Feature x, Ms1Feature y)
            {
                return (y.Abundance.CompareTo(x.Abundance));
            }
        }
    }
     */
}
