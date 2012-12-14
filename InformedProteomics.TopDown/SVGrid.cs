using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend;
using InformedProteomics.Backend.Scoring;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics
{
    public class SequenceVariantSet
    {
        #region Private Members

        private readonly Sequence _sequence;
        private readonly Composition[] _prefixCompArr;

        private readonly ScoredSpectra _scorer;
        private readonly ModificationParams _modParams;
        private readonly int _maxModIndex;

        /// <summary>
        /// Peak score of modIndex i
        /// </summary>
        private float?[] precursorScore;

        /// <summary>
        /// THe score of the best-scoring variant of _sequence
        /// </summary>
        private float? _bestScore;

        /// <summary>
        /// the best score among SVs corresponding to modIndex i
        /// </summary>
        private float?[] _bestScoreOfModIndex;

        /// <summary>
        /// 
        /// </summary>
        private SequenceVariant[] _bestScoringSVs;

        #endregion


        public SequenceVariantSet(Sequence sequence, ModificationParams modParams, ScoredSpectra scorer)
        {
            _sequence = sequence;
            _prefixCompArr = new Composition[_sequence.Count - 1];
            _prefixCompArr[0] = Composition.Zero;
            for (int i = 1; i < _sequence.Count - 1; i++)
                _prefixCompArr[i] = _prefixCompArr[i - 1] + _sequence[i].Composition;

            _modParams = modParams;
            _scorer = scorer;

            _maxModIndex = _modParams.GetNumMassShifts();

            _bestScoreOfModIndex = new float?[_maxModIndex];

            // Compute best score
            for (int modIndex = 0; modIndex < _maxModIndex; modIndex++)
                ComputeBestScore(modIndex);
        }

        public SequenceVariantSet()
        {
            throw new NotImplementedException();
        }

        public float? GetBestScore()
        {
            return _bestScore;
        }

        public float? GetBestScore(int modIndex)
        {
            return _bestScoreOfModIndex[modIndex];
        }

        public SequenceVariant[] GetBestScoringSequenceVariants()
        {
            return _bestScoringSVs;
        }

        private float GetPrecursorScore(int modIndex)
        {

            return 0f;
        }

        /// <summary>
        /// Compute the best score among SVs corresponding to _sequence + modIndex
        /// </summary>
        /// <param name="modIndex">Modification index applied to the sequence</param>
        private void ComputeBestScore(int modIndex)
        {
            Composition seqComp = _sequence.Composition + _modParams.GetModificationComposition(modIndex);
            float?[,] _score = new float?[_sequence.Count - 1, modIndex + 1];

            float?[] maxSeqScore = new float?[modIndex + 1];

            maxSeqScore[0] = 0f;
            for (int i = 0; i <= modIndex; i++)
            {
                _score[0, i] = 0f;
                _score[_sequence.Count - 1, i] = 0f;
            }

            for(int i=0; i<_sequence.Count; i++)
            {
                SequenceLocation location = SequenceLocation.Inner;
                if (i == 1)
                {
                    if (_sequence[0] == AminoAcid.PeptideNTerm)
                        location = SequenceLocation.PeptideNTerm;
                    else if (_sequence[0] == AminoAcid.ProteinNTerm)
                        location = SequenceLocation.ProteinNTerm;
                }
                else if (i == _sequence.Count - 2)
                {
                    if (_sequence[_sequence.Count - 1] == AminoAcid.PeptideCTerm)
                        location = SequenceLocation.PeptideCTerm;
                    else if (_sequence[_sequence.Count - 1] == AminoAcid.ProteinCTerm)
                        location = SequenceLocation.ProteinCTerm;
                }

                float?[] newMaxSeqScore = new float?[modIndex + 1];

                for (int prevModIndex = 0; prevModIndex <= modIndex; prevModIndex++)
                {
                    if (maxSeqScore[prevModIndex] != null)
                    {
                        int[] curModIndices = _modParams.GetModificationIndices(prevModIndex, _sequence[i], location);
                        foreach (int curModIndex in curModIndices)
                        {
                            float? cutScore = _score[i,curModIndex];
                            if (cutScore == null)
                            {
                                cutScore = _scorer.GetProductIonScore(_prefixCompArr[i] + _modParams.GetModificationComposition(curModIndex), seqComp);
                                _score[i, curModIndex] = cutScore;
                                newMaxSeqScore[curModIndex] = maxSeqScore[prevModIndex] + cutScore;
                            }
                            else
                            {
                                float seqScore = maxSeqScore[prevModIndex].Value + cutScore.Value;
                                if (seqScore > newMaxSeqScore[curModIndex])
                                    newMaxSeqScore[curModIndex] = seqScore;
                            }
                        }
                    }
                }
                maxSeqScore = newMaxSeqScore;
            }

            _bestScoreOfModIndex[modIndex] = maxSeqScore[modIndex];
        }
    }
}