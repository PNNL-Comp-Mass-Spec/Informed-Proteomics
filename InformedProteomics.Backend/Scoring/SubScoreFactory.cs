using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace InformedProteomics.Backend.Scoring
{
    public static class SubScoreFactory
    {
        class Node 
        {
            public IonType Ion { get; private set; }
            public float AttachementCost { get; set; }

            public Node(IonType ion)
            {
                Ion = ion;
            }

            public override int GetHashCode()
            {
 	             return Ion.GetHashCode();
            }

            public override bool Equals(object obj)
            {
                if (!(obj is Node)) return false;
                var o = obj as Node;
                return Ion.Equals(o.Ion);
            }

        }

        private static IEnumerable<float> GetmST(Dictionary<Node, Dictionary<Node, float>> edgeSet)
        {
            var queue = new PriorityQueue<Node, float>();
            var mST = new List<float>();
            var nodesNotInTree = new List<Node>();
            var connectedNode = new Dictionary<Node, Node>();
            
            foreach (var node in edgeSet.Keys)
            {
                node.AttachementCost = node.Ion.IsPrecursor ? 0f : -float.MaxValue;
                nodesNotInTree.Add(node);
                queue.Enqueue(node, node.AttachementCost);
            }

            while (nodesNotInTree.Count>0)
            {
                var node = queue.Dequeue();
                Node cnode;
                if (connectedNode.TryGetValue(node, out cnode))
                {
                    var e = edgeSet[cnode][node];
                    mST.Add(e);
                }
                //Console.WriteLine(node.Ion + "\t" + nodesNotInTree.Count);
                nodesNotInTree.Remove(node);

                foreach (var n in nodesNotInTree)
                {
                    //Console.WriteLine(node.Ion + "\t" + n.Ion);
                    var weight = edgeSet[node][n];
                    if (weight < n.AttachementCost) continue;
                    queue.Enqueue(n, weight);
                    connectedNode[n] = node;
                }
            }
            return mST;
        }


        public static void Read(string fileName)
        {
            _ionProbabilityScore = new Dictionary<FragmentSpectrumParameter, Dictionary<IonType, float>>();
            _ratioProbabilityScore = new Dictionary<FragmentSpectrumParameter, Dictionary<IonType, Dictionary<IonType, Dictionary<int, float>>>>();
            var ions = new List<IonType>();
            var sr = new StreamReader(fileName);
            string s;
            var phase = 0;
            FragmentSpectrumParameter par = null;
            int ionIndex1 = 0, ionIndex2 = 0;

            while ((s = sr.ReadLine()) != null)
            {
                var t = s.Split('\t');
                if (s.StartsWith("##"))
                {
                    phase = int.Parse(t[1]);
                    continue;
                }
                if (phase == 1)
                {
                    if (s.StartsWith("ION"))
                    {
                        var ion = new IonType(t[2], int.Parse(t[1]));
                        ions.Add(ion);
                    }
                    else
                    {
                        par = FragmentSpectrumParameter.Parse(t[0]);
                        if(!_ionProbabilityScore.ContainsKey(par))
                            _ionProbabilityScore[par] = new Dictionary<IonType, float>();
                        for (var i = 0; i < ions.Count; i++)
                            _ionProbabilityScore[par][ions[i]] = float.Parse(t[i + 1]);
                    }
                }else if (phase == 2)
                {
                    if (s.StartsWith("PARA"))
                    {
                        par = FragmentSpectrumParameter.Parse(t[1]);
                    }else if (s.StartsWith("ION"))
                    {
                        ionIndex1 = int.Parse(t[1]);
                        ionIndex2 = int.Parse(t[3]);
                    }
                    else if (par!=null)
                    {
                        var ion1 = ions[ionIndex1];
                        var ion2 = ions[ionIndex2];
                        if(!_ratioProbabilityScore.ContainsKey(par))
                            _ratioProbabilityScore[par] = new Dictionary<IonType, Dictionary<IonType, Dictionary<int, float>>>();
                        if(!_ratioProbabilityScore[par].ContainsKey(ion1))
                            _ratioProbabilityScore[par][ion1] = new Dictionary<IonType, Dictionary<int, float>>();
                        if(!_ratioProbabilityScore[par][ion1].ContainsKey(ion2))
                            _ratioProbabilityScore[par][ion1][ion2] = new Dictionary<int, float>();

                        foreach (var t1 in t)
                        {
                            if (t1.Length == 0) continue;
                            var st = t1.Split(':');
                            var r = int.Parse(st[0]);
                            _ratioProbabilityScore[par][ion1][ion2][r] = float.Parse(st[1]);
                        }
                    }
                }
            }
            sr.Close();
        }

        

        private static float GetNormalLogRatio(float x, float numVar)
        {
            x += (numVar-1)*.1f;
            var sigma = (float)Math.Sqrt(0.2);
            const float mean1 = 0.9f;
            const float mean2 = 0.2f;
            return (float)(-.5f * (Math.Pow((x - mean1) / sigma, 2) - Math.Pow((x - mean2) / sigma, 2)))*10;
        }
         
        internal static float GetPrecursorIonLikelihoodRatioScore(float rawScore, float numVar) // spectrum para
        {
            return GetNormalLogRatio(rawScore, numVar);
        }

        internal static float GetProductIonXICLikelihoodRatioScore(float rawScore, float numVar, FragmentParameter par) // spectrum para
        {
            return GetNormalLogRatio(rawScore, numVar);
        }

        internal static float GetPrecursorIonCorrelationCoefficient(int c1, int c2) // spectrum para
        {
            if (c1 == c2) return 1;
            return 0.8f;
        }

        internal static float GetIonXICCorrelationCoefficient(IonType ion1, IonType ion2, FragmentParameter par) // spectrum para, peak para (of ion1)
        {
            if (ion1.Equals(ion2)) return 1;
            if (ion1.Charge != ion2.Charge) return 0.4f;
            if (ion1.IsPrefix == ion2.IsPrefix) return 0.7f;
           
            return 0.7f;
            
        }

        internal static float GetIonXICCorrelationCoefficient(IonType ion, FragmentParameter par) // spectrum para, peak para
        {
            return 0.7f;
        }

        public static float GetProductIonSpectrumScore(FragmentSpectrum spectrum, FragmentSpectrumParameter par)
        {
            var score = 0f;
            var maxCharge = (from ion in spectrum.Keys where !ion.IsPrecursor select ion.Charge).Concat(new[] {0}).Max();
            /* foreach(var ion in spectrum.Keys)
            {
                if (ion.IsPrecursor) continue;
                if (maxCharge < ion.Charge) maxCharge = ion.Charge;
            }*/
            for (var charge = 1; charge <= maxCharge; charge++)
            {
                score += GetProductIonSpectrumScore(spectrum, par, charge);
            }
            return score;
        }

        private static float GetProductIonSpectrumScore(FragmentSpectrum spectrum,
                                                       FragmentSpectrumParameter par, int charge)
        {
            var edgeSet = new Dictionary<Node, Dictionary<Node, float>>();
            var s = new Node(spectrum.PrecursorIon);
            edgeSet[s] = new Dictionary<Node, float>();
            foreach (var ion1 in _ionProbabilityScore[par].Keys)
            {
                if (ion1.Charge != charge) continue;
                var l = new Node(ion1);
                edgeSet[s][l] = GetRatioScore(spectrum, spectrum.PrecursorIon, ion1, par);
                if (!edgeSet.ContainsKey(l)) edgeSet[l] = new Dictionary<Node, float>();
                foreach (var ion2 in _ionProbabilityScore[par].Keys)
                {
                    if (ion1.Equals(ion2) || ion2.Charge != charge) continue;
                    var r = new Node(ion2);
                    edgeSet[l][r] = GetRatioScore(spectrum, ion1, ion2, par);
                }
            }

            var edgeList = GetmST(edgeSet);
            //foreach (var e in edgeList) Console.Write(e + " ");
            //Console.WriteLine(spectrum.Keys.Count);
            return edgeList.Sum(e => e);
        }

        private static
            Dictionary<FragmentSpectrumParameter, Dictionary<IonType, Dictionary<IonType, Dictionary<int, float>>>>
            _ratioProbabilityScore;

        private static Dictionary<FragmentSpectrumParameter, Dictionary<IonType, float>> _ionProbabilityScore;

        private static int GetRatio(double v1, double v2)
        {
            if (v1 <= 0)
            {
                if (v2 <= 0) return -11;
                return -12;
            }
            if (v2 <= 0) return 11;
            double r;
            var f = 1;
            if (v1 > v2) r = v1/v2;
            else
            {
                r = v2/v1;
                f = -1;
            }
            return (int) (Math.Min(10, r)*f);
        }

        private static float GetRatioScore(FragmentSpectrum spectrum, IonType ion1, IonType ion2,
                                           FragmentSpectrumParameter par)
        {
            if (ion1.IsPrecursor)
                return _ionProbabilityScore[par][ion2];
            var r = GetRatio(spectrum.GetIntensity(ion1), spectrum.GetIntensity(ion2));
            //Console.WriteLine(par.ToFileString() + "\t" + ion1 + "\t" + ion2 + "\t" + r);
            return _ratioProbabilityScore[par][ion1][ion2][r];
        }

       /* public static float GetIonXICMean(FragmentIonClassBase FragmentIonClassBase, GroupParameter par)
        {
            return 0.6f;
        }

        public static float GetIonXICVariance(FragmentIonClassBase FragmentIonClassBase, GroupParameter par)
        {
            return 1;
        }*/
    }
}
