using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Science
{
    public interface IMolecule : IMatter
    {
        Composition GetComposition();
    }
}
