using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.SequenceTag;

namespace InformedProteomics.TopDown.TagBasedSearch
{
    public class SequenceTagGenerator : ISequenceTagFinder
    {
        public SequenceTagGenerator(LcMsRun run, Tolerance tolerance, int minTagLength = 5, int maxTagLength = 8,
            AminoAcid[] aminoAcidsArray = null)
        {
            _run = run;
            _tolerance = tolerance;

            _minTagLen = minTagLength;
            _maxTagLen = maxTagLength;
            _aminoAcids = aminoAcidsArray ?? AminoAcid.StandardAminoAcidArr;
            _ms2ScanToTagMap = new Dictionary<int, IList<SequenceTag.SequenceTag>>();
        }

        public void SearchAllSpectra()
        {
            foreach (var ms2ScanNum in _run.GetScanNumbers(2))
            {
                Generate(ms2ScanNum);
            }
        }

        public void Generate(int ms2ScanNum)
        {
            var spec = _run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            if (spec == null) return;
            var tagFinder = new SequenceTagFinder(spec, _tolerance, _minTagLen, _maxTagLen, _aminoAcids);

            var tags = tagFinder.GetAllSequenceTagString();

            lock (_ms2ScanToTagMap)
            {
                _ms2ScanToTagMap[ms2ScanNum] = tags;
            }
        }

        public IList<SequenceTag.SequenceTag> GetAllSequenceTagString(int ms2ScanNum)
        {
            IList<SequenceTag.SequenceTag> tags;

            lock (_ms2ScanToTagMap)
            {
                if (_ms2ScanToTagMap.TryGetValue(ms2ScanNum, out tags))
                {
                    return tags;
                }
            }

            var spec = _run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            if (spec == null) return new List<SequenceTag.SequenceTag>();
            var tagFinder = new SequenceTagFinder(spec, _tolerance, _minTagLen, _maxTagLen, _aminoAcids);
            tags = tagFinder.GetAllSequenceTagString();

            lock (_ms2ScanToTagMap)
            {
                _ms2ScanToTagMap[ms2ScanNum] = tags;
            }

            return tags;
        }

        public long NumberOfGeneratedTags()
        {
            return _ms2ScanToTagMap.Values.Sum(tags => tags.Count);
        }

        public IEnumerable<int> GetMs2ScanNumsContainingTags()
        {
            return _ms2ScanToTagMap.Keys;
        }

        //private readonly List<int> _ms2ScansContainingTags;
        private readonly LcMsRun _run;
        private readonly AminoAcid[] _aminoAcids;
        private readonly Tolerance _tolerance;
        private readonly int _minTagLen;
        private readonly int _maxTagLen;
        private readonly Dictionary<int, IList<SequenceTag.SequenceTag>> _ms2ScanToTagMap;
    }
}
