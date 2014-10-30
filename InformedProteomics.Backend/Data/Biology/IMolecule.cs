namespace InformedProteomics.Backend.Data.Biology
{
    public interface IMolecule : IMatter
    {
        Composition.Composition Composition { get; }
    }
}
