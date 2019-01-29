using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.PostProcessing
{
    public class SequenceTagParser
    {
        public SequenceTagParser(string tagFilePath, int minTagLength = 0, int numTagsPerScan = int.MaxValue)
        {
            TagFilePath = tagFilePath;
            _scanToTags = new Dictionary<int, IList<SequenceTag.SequenceTag>>();
            _minTagLength = minTagLength;
            _numTagsPerScan = numTagsPerScan;
            Parse(tagFilePath);
        }

        public IEnumerable<SequenceTag.SequenceTag> GetSequenceTags(int scanNum)
        {
            if (_scanToTags.TryGetValue(scanNum, out var tags)) return tags;
            return new List<SequenceTag.SequenceTag>();
        }

        public IEnumerable<int> GetScanNums()
        {
            return _scanToTags.Keys.OrderBy(s => s);
        }

        public string TagFilePath { get; }
        private readonly int _minTagLength;
        private readonly int _numTagsPerScan;
        private readonly Dictionary<int, IList<SequenceTag.SequenceTag>> _scanToTags;

        private void Parse(string tagFilePath)
        {
            var tagParser = new TsvFileParser(tagFilePath);
            var scan = tagParser.GetData("ScanNum").Select(s => Convert.ToInt32(s)).ToArray();
            var sequence = tagParser.GetData("SequenceTag").ToArray();
            var isPrefix = tagParser.GetData("IsPrefix").Select(s => s.Equals("1")).ToArray();
            var flankingMass = tagParser.GetData("FlankingMass").Select(Convert.ToDouble).ToArray();
            for (var i = 0; i < tagParser.NumData; i++)
            {
                if (sequence[i].Length < _minTagLength) continue;
                var tag = new SequenceTag.SequenceTag(scan[i], sequence[i], isPrefix[i], flankingMass[i]);

                if (_scanToTags.TryGetValue(scan[i], out var tagList))
                {
                    if (tagList.Count < _numTagsPerScan) tagList.Add(tag);
                }
                else
                {
                    _scanToTags.Add(scan[i], new List<SequenceTag.SequenceTag> { tag });
                }
            }
        }
    }
    /*
    public class SequenceTag
    {
        public SequenceTag(int scanNum, string sequence, bool isPrefix, double flankingMass)
        {
            ScanNum = scanNum;
            Sequence = sequence;
            IsPrefix = isPrefix;
            FlankingMass = flankingMass;
        }

        public int ScanNum { get; private set; }
        public string Sequence { get; private set; }
        public bool IsPrefix { get; private set; }
        public double FlankingMass { get; private set; }

        public double TagMass
        {
            get
            {
                return _tagMass ??
                       (double)(_tagMass = AminoAcidSet.GetStandardAminoAcidSet().GetComposition(Sequence).Mass);
            }
        }

        public double? GetNTermFlankingMass(double? sequenceMass)
        {
            return GetNTermFlankingMass(sequenceMass, BaseIonType.B, BaseIonType.Y);
        }

        public double? GetNTermFlankingMass(double? sequenceMass, BaseIonType prefixIonType, BaseIonType suffixIonType)
        {
            if(IsPrefix) return FlankingMass - prefixIonType.OffsetComposition.Mass;
            if (sequenceMass == null) return null;
            return sequenceMass - (FlankingMass - suffixIonType.OffsetComposition.Mass) - TagMass;
        }

        public double? GetCTermFlankingMass(double? sequenceMass)
        {
            return GetCTermFlankingMass(sequenceMass, BaseIonType.B, BaseIonType.Y);
        }

        public double? GetCTermFlankingMass(double? sequenceMass, BaseIonType prefixIonType, BaseIonType suffixIonType)
        {
            if (!IsPrefix) return FlankingMass - suffixIonType.OffsetComposition.Mass;
            if (sequenceMass == null) return null;
            return sequenceMass - (FlankingMass - prefixIonType.OffsetComposition.Mass) - TagMass;
        }

        private double? _tagMass;
    }
    */
}
