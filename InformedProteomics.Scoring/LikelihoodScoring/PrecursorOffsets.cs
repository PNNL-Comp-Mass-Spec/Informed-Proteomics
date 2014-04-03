using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class PrecursorOffsets
    {
        public int Charge { get; private set; }

        private readonly Dictionary<int, List<double>> _offsets;
        private readonly double _probabilityThreshold;

        public PrecursorOffsets(int charge, double probabilityThreshold=0.0)
        {
            Charge = charge;
            _offsets = new Dictionary<int, List<double>>();
            _probabilityThreshold = probabilityThreshold;

            for (int i = 1; i <= charge; i++)
            {
                _offsets.Add(i, new List<double>());
            }
        }

        public void AddOffsets(int charge, IEnumerable<OffsetProbability> offsetProbabilities)
        {
            if (charge > Charge)
            {
                throw new Exception("charge must be less than or equal precursor charge.");
            }

            var offsets = (from offsetProb in offsetProbabilities
                           where offsetProb.Prob >= _probabilityThreshold
                           select offsetProb.Prob).ToList();

            _offsets.Add(charge, offsets);
            _offsets[charge] = _offsets[charge].Distinct().ToList(); // remove duplicates
        }

        public void AddOffsetsFromFile(string fileName)
        {
            var offsetFile = new TsvFileParser(fileName);
            var offsets = offsetFile.GetData("Offset");
            var offsetFileData = offsetFile.GetAllData();

            foreach (var chargeListLabel in offsetFileData.Keys)
            {
                if (chargeListLabel == "Offset") continue;
                var chargeStr = chargeListLabel.Remove(0, 3); // remove "Ch "
                var charge = Convert.ToInt32(chargeStr);
                var chargeList = offsetFileData[chargeListLabel];
                var offsetList = new List<double>();
                for (int i = 0; i < chargeList.Count; i++)
                {
                    var offset = Convert.ToDouble(offsets[i]);
                    var probability = Convert.ToDouble(chargeList[i]);
                    if (probability >= _probabilityThreshold)
                        offsetList.Add(offset);
                }
                _offsets[charge].AddRange(offsetList);
                _offsets[charge] = _offsets[charge].Distinct().ToList(); // remove duplicates
            }
        }

        public List<double> GetChargeOffsets(int charge)
        {
            return _offsets[charge];
        }
    }
}
