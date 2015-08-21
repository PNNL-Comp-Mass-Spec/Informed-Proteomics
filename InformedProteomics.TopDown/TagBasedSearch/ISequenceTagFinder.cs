using System.Collections.Generic;
using InformedProteomics.Backend.SequenceTag;

namespace InformedProteomics.TopDown.TagBasedSearch
{
    public interface ISequenceTagFinder
    {
        IList<SequenceTagString> GetAllSequenceTagString(int ms2ScanNum);
        long NumberOfGeneratedTags();

        //IEnumerable<int> GetMs2ScanNumsContainingTags();
    }
}
