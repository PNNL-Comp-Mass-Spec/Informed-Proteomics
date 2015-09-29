using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Remoting.Messaging;

namespace InformedProteomics.Backend.SequenceTag
{
    public class SequenceTagGraph<T> where T : GraphEdge
    {
        private List<T>[] _adjList;
        private bool[] _hasInEdge;

        public SequenceTagGraph(int maxTagLen = 8)
        {
            MaxTagLen = maxTagLen;
        }

        public SequenceTagGraph(int nodeCount, int maxTagLen = 8)
        {
            SetNodeCount(nodeCount);
            MaxTagLen = maxTagLen;
        }

        public void SetNodeCount(int nodeCount)
        {
            _adjList = new List<T>[nodeCount];
            _hasInEdge = new bool[nodeCount];
            for (var i = 0; i < nodeCount; i++)
            {
                _adjList[i] = new List<T>();
            }            
        }

        public void AddEdge(T edge)
        {
            _adjList[edge.Node1].Add(edge);
            _hasInEdge[edge.Node2] = true;
        }

        public bool HasInEdge(int node)
        {
            return _hasInEdge[node];
        }

        public int GetEdgeCount()
        {
            return _adjList.Sum(edges => edges.Count);
        }

        public int GetNodeCount()
        {
            return _adjList.Length;
        }

        public IEnumerable<int> Nodes()
        {
            return Enumerable.Range(0, _adjList.Length);
        }

        public IEnumerable<int> Roots()
        {
            for (var i = 0; i < _hasInEdge.Length; i++)
            {
                if (_hasInEdge[i] == false) yield return i;
            }
        }

        public bool HasEdge(int node1, int node2)
        {
            return _adjList[node1].Any(edge => edge.Node2 == node2);
        }

        public void RemoveEdge(int node1, int node2)
        {
            for(var i = 0; i < _adjList[node1].Count; i++)
            {
                if (_adjList[node1][i].Node2 == node2)
                {
                    _adjList[node1].RemoveAt(i);
                    return;
                }
            }
        }

        public List<T> OutEdges(int node)
        {
            return _adjList[node];
        }

        public List<T> InEdges(int node)
        {
            return (from t in _adjList from edge in t where edge.Node2 == node select edge).ToList();
        }


        public List<List<int>> ConnnectedComponents()
        {
            /*
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
            }*/

            var componentSet = new List<List<int>>();
            var visited = new bool[GetNodeCount()];
            for (var i = 0; i < GetNodeCount(); i++)
            {
                if (visited[i]) continue;

                var component = new List<int>();
                var neighbors = new Queue<int>();
                neighbors.Enqueue(i);
                while (true)
                {
                    if (neighbors.Count < 1) break;
                    var j = neighbors.Dequeue();
                    if (visited[j]) continue;
                    visited[j] = true;

                    component.Add(j);
                    foreach (var edge in _adjList[j])
                    {
                        if (visited[edge.Node2]) continue;
                        neighbors.Enqueue(edge.Node2);
                    }
                }
                componentSet.Add(component);
            }
            return componentSet;
        }


        protected virtual bool ProcessPath(IEnumerable<T> edges)
        {
            foreach (var edge in edges)
            {
                Console.Write(edge.Node1 + "-" + edge.Node2 + "\t");
            }
            Console.Write("\n");

            return true;
        }
        /*
        protected virtual bool AlreadyUsedPeak(int node)
        {
            return false;
        }*/

        protected bool StopFindPath;
        protected readonly Stack<T> EdgeList = new Stack<T>();
        
        protected bool[] NodeVisitFlag;
        protected int MaxTagLen;
        public void FindPaths(int node, bool firstCall = true, T e = null)
        {
            if (firstCall)
            {
                NodeVisitFlag = new bool[_adjList.Length];
                //EdgeList = new Stack<T>();
                EdgeList.Clear();
            }
            else
            {
                if (EdgeList.Count >= MaxTagLen)
                {
                    ProcessPath(EdgeList.Reverse());
                    EdgeList.Clear();
                    return;
                }

                //if (NumberOfProcessedPaths > MaxNumberOfProcessedPaths)
                if (StopFindPath)
                {
                    EdgeList.Clear();
                    return;
                }
            }

            if (e != null) EdgeList.Push(e);
            NodeVisitFlag[node] = true;
            
            var flag = false;
            foreach(var edge in OutEdges(node))
            {
                //if (!NodeVisitFlag[edge.Node2] && !AlreadyUsedPeak(edge.Node2))
                if (!NodeVisitFlag[edge.Node2])
                {
                    flag = true;
                    FindPaths(edge.Node2, false, edge);
                }
            }
            
            if(!flag)
            {
                var t = ProcessPath(EdgeList.Reverse());
                if (t == false) return;
            }
            NodeVisitFlag[node] = false;
            
            if (EdgeList.Count > 0) EdgeList.Pop();
        }
    }
}
