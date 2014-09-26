using System.Collections.Generic;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDown.Scoring
{
    class TopDownTrainer
    {
        private Dictionary<int, int> _missingXicCounterTarget; // charge, count
        private Dictionary<int, Dictionary<int, int>> _xicCorrScoreCounterTarget;    // charge, correlation raw score, count
        private Dictionary<int, int> _missingXicCounterDecoy; // charge, count
        private Dictionary<int, Dictionary<int, int>> _xicCorrScoreCounterDecoy;    // charge, correlation raw score, count
        private InMemoryLcMsRun _run;
        private Tolerance _tolerance;

        public void Train(string outFileName, InMemoryLcMsRun run, Tolerance tolerance, string annotationFileName)
        {
            // charge, scan number, protein
            _run = run;
            _tolerance = tolerance;
            var target = new Dictionary<int, Dictionary<int, Composition>>();
            var decoy = new Dictionary<int, Dictionary<int, Composition>>();

            _missingXicCounterTarget = new Dictionary<int, int>();
            _missingXicCounterDecoy = new Dictionary<int, int>();
            _xicCorrScoreCounterTarget = new Dictionary<int, Dictionary<int, int>>();
            _xicCorrScoreCounterDecoy = new Dictionary<int, Dictionary<int, int>>();
        }

        private void SubTrain(Dictionary<int, Dictionary<int, Composition>> dictionary, bool isTarget)
        {
            Dictionary<int, int> missingXicCounter;
            Dictionary<int, Dictionary<int, int>> corrScoreCounter;

            if (isTarget)
            {
                missingXicCounter = _missingXicCounterTarget;
                corrScoreCounter = _xicCorrScoreCounterTarget;
            }
            else
            {
                missingXicCounter = _missingXicCounterDecoy;
                corrScoreCounter = _xicCorrScoreCounterDecoy;
            }

            foreach (var charge in dictionary.Keys)
            {
                var subDictionary = dictionary[charge];
                foreach (var scanNumber in subDictionary.Keys)
                {
                    var composition = subDictionary[scanNumber];
                    var scorer = new TopDownScorer(composition, _run, _tolerance, null);
                    var areXicMissing = scorer.AreXicMissing(charge, scanNumber);
                    var corrScores = scorer.GetIsotopeCorrelationIntensityRawScores(charge, scanNumber);
                    for (var c = TopDownScorer.MinCharge; c <= TopDownScorer.MaxCharge;c++ )
                    {
                        if (areXicMissing[c - TopDownScorer.MinCharge])
                        {
                            var t = missingXicCounter[charge];
                            missingXicCounter[charge] = t + 1;
                        }
                        else
                        {
                            var t = corrScoreCounter[c] ?? new Dictionary<int, int>();
                            var score = corrScores[c - TopDownScorer.MinCharge];
                            t[score] = t[score] + 1;
                        }
                    }
                }
            }
        }
    }
}
 