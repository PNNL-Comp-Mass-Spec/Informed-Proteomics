using System.Collections.Generic;


namespace InformedProteomics.Backend.MassFeature
{
    public class NodeSet<T> : List<T>
    {
        public IList<NodeSet<T>> ConnnectedComponents(INodeComparer<T> comparer)
        {
            SetAdjacentList(comparer);

            var componentSet = new List<NodeSet<T>>();
            var visited = new bool[Count];
            for (var i = 0; i < Count; i++)
            {
                if (visited[i]) continue;

                var component = new NodeSet<T>();
                var neighbors = new Queue<int>();
                neighbors.Enqueue(i);
                while (true)
                {
                    if (neighbors.Count < 1) break;
                    var j = neighbors.Dequeue();
                    if (visited[j]) continue;
                    visited[j] = true;

                    component.Add(this[j]);
                    foreach (var nd in _adjList[j])
                    {
                        if (visited[nd]) continue;
                        neighbors.Enqueue(nd);
                    }
                }
                componentSet.Add(component);
            }
            return componentSet;
        }

        private void SetAdjacentList(INodeComparer<T> comparer)
        {
            _adjList = new List<int>[Count];
            for (var i = 0; i < Count; i++)
            {
                _adjList[i] = new List<int>();
            }

            for (var i = 0; i < Count; i++)
            {
                for (var j = i + 1; j < Count; j++)
                {
                    if (comparer.SameCluster(this[i], this[j]))
                    {
                        _adjList[i].Add(j);
                        _adjList[j].Add(i);
                    }
                }
            }
        }


        private List<int>[] _adjList;
     }
}
