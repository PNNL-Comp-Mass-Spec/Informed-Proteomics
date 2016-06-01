using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using System.Diagnostics;
using System.IO;

namespace InformedProteomics.TopDown.Quantification
{
    public class MSDeconvClusterer
    {

        public MSDeconvClusterer(LcMsRun run)
        {
            _run = run;
        }

        public List<List<MSDeconvNode>> GetClustersList(Tolerance tol, double elutionInterval,List<MSDeconvNode> nodeList)
        {
            var edgeDict = GetEdges(nodeList,tol,elutionInterval);

            var connectedComponents = new List<List<MSDeconvNode>>();
            var nodesInComponenets = new Dictionary<MSDeconvNode,int>();

            for (int i = 0; i < nodeList.Count; i++)
            {
                if (!nodesInComponenets.ContainsKey(nodeList[i]))
                {
                    connectedComponents.Add(GetConnectedComponenet(nodeList[i],edgeDict,nodesInComponenets));
                }
            }

            /**foreach (var c in connectedComponents)
            {
                Console.WriteLine("Connected Componenet Elements:");
                foreach (var n in c)
                {
                    Console.WriteLine("{0} {1}", n.ScanNumber,n.RealMonoMass);
                }
                Console.WriteLine();
            }*/

            return connectedComponents;
        }

        public void SaveClusterInfoToFile(String filename, List<List<MSDeconvNode>> clusters)
        {
            var tsv = new StringBuilder();
            tsv.Append("FeatureID \t MinScan \t MaxScan \t MinCharge \t MaxCharge \t MonoMass \t Abundance \t MinElutionTime \t MaxElutionTime \t \n");
            var id = 1;
            foreach (var c in clusters)
            {
                var maxScan = c.Max(x => x.ScanNumber);
                var minScan = c.Min(x => x.ScanNumber);
                var minCharge = c.Min(x => x.Charge);
                var maxCharge = c.Max(x => x.Charge);
                var abundance = c.Sum(x => x.RealIntensitySum);
                var maxScanNum = c.Max(x => x.ScanNumber);
                var minScanNum = c.Min(x => x.ScanNumber);

                var sortedMono = c.OrderBy(x => x.RealMonoMass);
                var halfIndex = c.Count() / 2;
                var medianMass = c[0].RealMonoMass;
                if (c.Count()%2 == 0)
                {
                    medianMass = (sortedMono.ElementAt(halfIndex).RealMonoMass +
                              sortedMono.ElementAt(halfIndex - 1).RealMonoMass)/2;
                }
                else medianMass = sortedMono.ElementAt(halfIndex).RealMonoMass;

                var info = String.Format("{0} \t {1} \t {2} \t {3} \t {4} \t {5} \t {6} \t {7} \t {8} \t \n", id, minScan, maxScan, minCharge, maxCharge, medianMass, abundance, _run.GetElutionTime(minScanNum), _run.GetElutionTime(maxScanNum));
                tsv.Append(info);
                id++;
            }
            File.WriteAllText(filename+".tsv", tsv.ToString());
        }

        private Dictionary<MSDeconvNode, List<MSDeconvNode>> GetEdges(List<MSDeconvNode> nodeList, Tolerance tol, double elutionInterval)
        {
            var edgeDict = new Dictionary<MSDeconvNode, List<MSDeconvNode>>();
            var currentNode = nodeList[0];
            var currentIndex = GetIndexRange(nodeList, elutionInterval, currentNode);
            foreach (var n1 in nodeList)
            {
                if (n1.ScanNumber != currentNode.ScanNumber)
                {
                    currentNode = n1;
                    currentIndex = GetIndexRange(nodeList, elutionInterval, currentNode);
                    if (currentIndex.Item2 == -1) break;
                }
                edgeDict.Add(n1 ,new List<MSDeconvNode>());
                for (var i = currentIndex.Item1 ; i < currentIndex.Item2; i++)
                {
                    var deltaMass = Math.Abs(nodeList[i].RealMonoMass - n1.RealMonoMass);
                    if(deltaMass <= tol.GetValue()) edgeDict[n1].Add(nodeList[i]);
                }

            }
            return edgeDict;
        }

        private List<MSDeconvNode> GetConnectedComponenet(MSDeconvNode node, Dictionary<MSDeconvNode, List<MSDeconvNode>> edges, Dictionary<MSDeconvNode,int> connectedNodes)
        {
            var connectedComp = new List<MSDeconvNode>();
            var toVisit = new List<MSDeconvNode>();
            toVisit.Add(node);
            while (toVisit.Count != 0)
            {
                var current = toVisit[0];
                toVisit.Remove(current);
                if (!connectedComp.Contains(current))
                {
                    connectedComp.Add(current);
                    if (connectedNodes.ContainsKey(current)) break;
                    connectedNodes.Add(current,0);
                    if (!edges.ContainsKey(current)) break;
                    foreach (var n in edges[current])
                    {
                        toVisit.Add(n);
                    }
                }
            }
            return connectedComp;
        }


        private Tuple<int, int> GetIndexRange(List<MSDeconvNode> nodeList, double elutionInterval, MSDeconvNode node)
        {

            var startIndex = -1;
            var endIndex = -1;
            var setStart = false;
            for (int i = nodeList.IndexOf(node); i < nodeList.Count-1; i++)
            {
                if(nodeList[i]. ScanNumber == node.ScanNumber && !setStart)
                {
                    setStart = true;
                    startIndex = i + 1;
                }

                var deltaElution = _run.GetElutionTime(nodeList[i+1].ScanNumber) - _run.GetElutionTime(node.ScanNumber);
                if (deltaElution > elutionInterval)
                {
                    endIndex = i + 1;
                    break;
                }
            }

            return new Tuple<int, int>(startIndex,endIndex);
        }

        private readonly LcMsRun _run;
    }
}
