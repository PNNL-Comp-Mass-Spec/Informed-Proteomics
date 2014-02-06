
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Biology
{
    public interface IMolecule : IMatter
    {
        Composition GetComposition();
    }
}
