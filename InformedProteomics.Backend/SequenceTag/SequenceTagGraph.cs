using System;
using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Backend.SequenceTag
{
    public class GraphEdge : IEquatable<GraphEdge>
    {
        public int Node1 { get; private set; }
        public int Node2 { get; private set; }

        public GraphEdge(int node1, int node2)
        {
            Node1 = node1;
            Node2 = node2;
        }

        public bool Equals(GraphEdge other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return Node1 == other.Node1 && Node2 == other.Node2;
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((GraphEdge)obj);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                return (Node1 * 397) ^ Node2;
            }
        }

        public static bool operator ==(GraphEdge edge1, GraphEdge edge2)
        {
            if ((object)edge1 == null || ((object)edge2) == null)
                return Object.Equals(edge1, edge2);

            return edge1.Equals(edge2);
        }

        public static bool operator !=(GraphEdge edge1, GraphEdge edge2)
        {
            if (edge1 == null || edge2 == null)
                return !Object.Equals(edge1, edge2);

            return !(edge1.Equals(edge2));
        }
    }

    public class SequenceTagGraph<T> where T : GraphEdge
    {
        private List<T>[] _adjList;
        private bool[] _hasInEdge;

        public SequenceTagGraph() { }

        public SequenceTagGraph(int nodeCount)
        {
            SetNodeCount(nodeCount);
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

        public void FindPaths(int node, bool firstCall = true, T e = null)
        {
            if (firstCall) 
            {
                _nodeVisitFlag = new bool[_adjList.Length];
                //_nodeList = new Stack<int>();
                _edgeList = new Stack<T>();
            }

            //_nodeList.Push(node);
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
            //_nodeList.Pop();
            if (_edgeList.Count > 0) _edgeList.Pop();
        }
    }
}
