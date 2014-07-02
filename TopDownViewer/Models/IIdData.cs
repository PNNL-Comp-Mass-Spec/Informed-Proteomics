namespace InformedProteomics.TopDownViewer.Models
{
    public interface IIdData
    {
        void Add(PrSm data);
        PrSm GetHighestScoringPrSm();
    }
}
