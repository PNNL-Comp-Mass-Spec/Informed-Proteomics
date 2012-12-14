namespace InformedProteomics.Backend.Utils
{
    public interface IMolecule : IMatter
    {
        Composition GetComposition();
    }
}
