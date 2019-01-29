using System;

namespace InformedProteomics.TopDown.SequenceTag
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
                return object.Equals(edge1, edge2);

            return edge1.Equals(edge2);
        }

        public static bool operator !=(GraphEdge edge1, GraphEdge edge2)
        {
            if (edge1 == null || edge2 == null)
                return !object.Equals(edge1, edge2);

            return !(edge1.Equals(edge2));
        }
    }
}
