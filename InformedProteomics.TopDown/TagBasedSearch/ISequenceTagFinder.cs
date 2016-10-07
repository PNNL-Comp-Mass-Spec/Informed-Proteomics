using System.Collections.Generic;

namespace InformedProteomics.TopDown.TagBasedSearch
{
    public interface ISequenceTagFinder
    {
        IList<SequenceTag.SequenceTag> GetAllSequenceTagString(int ms2ScanNum);
        long NumberOfGeneratedTags();

        //IEnumerable<int> GetMs2ScanNumsContainingTags();
    }
}
