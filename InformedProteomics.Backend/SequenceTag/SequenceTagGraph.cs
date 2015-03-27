using System;
using System.Collections.Generic;
using System.Linq;

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

        protected virtual void ProcessPath(IEnumerable<T> edges)
        {
            foreach (var edge in edges)
            {
                Console.Write(edge.Node1 + "-" + edge.Node2 + "\t");
            }
            Console.Write("\n");
        }

        //private Stack<int> _nodeList;
        private Stack<T> _edgeList;
        private bool[] _nodeVisitFlag;
        protected int MaxTagLen;
        public void FindPaths(int node, bool firstCall = true, T e = null)
        {
            if (firstCall)
            {
                _nodeVisitFlag = new bool[_adjList.Length];
                _edgeList = new Stack<T>();
            }
            else
            {
                if (_edgeList.Count >= MaxTagLen)
                {
                    ProcessPath(_edgeList.Reverse());
                    _edgeList.Clear();
                    return;
                }
            }

            if (e != null) _edgeList.Push(e);

            _nodeVisitFlag[node] = true;
            
            bool flag = false;
            foreach(var edge in OutEdges(node))
            {
                if (!_nodeVisitFlag[edge.Node2])
                {
                    flag = true;
                    FindPaths(edge.Node2, false, edge);
                }
            }
            
            if(!flag)
            {
                ProcessPath(_edgeList.Reverse());
            }
            _nodeVisitFlag[node] = false;
            
            if (_edgeList.Count > 0) _edgeList.Pop();
        }
    }
}
