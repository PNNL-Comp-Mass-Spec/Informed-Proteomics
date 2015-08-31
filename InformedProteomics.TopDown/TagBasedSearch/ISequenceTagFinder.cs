using System.Collections.Generic;
using InformedProteomics.Backend.SequenceTag;

namespace InformedProteomics.TopDown.TagBasedSearch
{
    public interface ISequenceTagFinder
    {
        IList<SequenceTag> GetAllSequenceTagString(int ms2ScanNum);
        long NumberOfGeneratedTags();

        //IEnumerable<int> GetMs2ScanNumsContainingTags();
    }
}
