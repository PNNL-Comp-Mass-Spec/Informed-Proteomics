namespace InformedProteomics.FeatureFinding.Util
{
    public interface INodeComparer<T>
    {
        bool SameCluster(T node1, T node2);
    }
}
