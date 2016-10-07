using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.SequenceTag;

namespace InformedProteomics.TopDown.PostProcessing
{
    public class MatchedTagSet
    {
        public MatchedTagSet(string sequence,
            AminoAcidSet aminoAcidSet, Tolerance tolerance, Tolerance relaxedTolerance)
        {
            _sequence = sequence;
            _aminoAcidSet = aminoAcidSet;
            _tolerance = tolerance;
            _relaxedTolerance = relaxedTolerance;
            _tags = new List<MatchedTag>();
        }

        public string Sequence { get { return _sequence; } }
        public List<MatchedTag> Tags { get { return _tags; } }

        /// <summary>
        /// Adds a tag to this tag set.
        /// </summary>
        /// <param name="tag">a matched tag to add</param>
        /// <returns>true if tag is merged to an existingTag tag. false otherwise</returns>
        public bool Add(MatchedTag tag)
        {
            if (_tags.Any(existingTag => TryMerge(existingTag, tag)))
            {
                return true;
            }
            _tags.Add(tag);
            return false;
        }

        private bool TryMerge(MatchedTag existingTag, MatchedTag newTag)
        {
            // N-term
            var newStartIndex = Math.Min(existingTag.StartIndex, newTag.StartIndex);
            double? adjustedNTermFlankingMass = null;
            if (existingTag.NTermFlankingMass != null && newTag.NTermFlankingMass != null)
            {
                var newNTermFlankingMassFromExistingTag = (double)existingTag.NTermFlankingMass -
                                           GetSequenceMass(newTag.StartIndex, existingTag.StartIndex);
                var newNTermFlankingMassFromNewTag = (double)newTag.NTermFlankingMass -
                                           GetSequenceMass(existingTag.StartIndex, newTag.StartIndex);

                var toleranceNTerm = existingTag.IsNTermFlankingMassReliable == newTag.IsNTermFlankingMassReliable
                    ? _tolerance
                    : _relaxedTolerance;
                if (!toleranceNTerm.IsWithin(newNTermFlankingMassFromExistingTag, newNTermFlankingMassFromNewTag)) return false;

                if (existingTag.IsNTermFlankingMassReliable || newTag.IsNTermFlankingMassReliable)    // flanking mass from this tag is reliable
                {
                    adjustedNTermFlankingMass =
                        (existingTag.NumReliableNTermFlankingMasses * newNTermFlankingMassFromExistingTag
                        + newTag.NumReliableNTermFlankingMasses * newNTermFlankingMassFromNewTag) /
                        (existingTag.NumReliableNTermFlankingMasses + newTag.NumReliableNTermFlankingMasses);
                }
                else
                {
                    adjustedNTermFlankingMass =
                        (existingTag.NumMergedSequenceTags * newNTermFlankingMassFromExistingTag
                        + newTag.NumMergedSequenceTags * newNTermFlankingMassFromNewTag) /
                        (existingTag.NumMergedSequenceTags + newTag.NumMergedSequenceTags);
                }
            }
            else if (existingTag.NTermFlankingMass != newTag.NTermFlankingMass) return false;

            // C-term
            var newEndIndex = Math.Max(existingTag.EndIndex, newTag.EndIndex);
            double? adjustedCTermFlankingMass = null;
            if (existingTag.CTermFlankingMass != null && newTag.CTermFlankingMass != null)
            {
                var newCTermFlankingMassFromExistingTag = (double) existingTag.CTermFlankingMass -
                                                          GetSequenceMass(existingTag.EndIndex, newTag.EndIndex);
                var newCTermFlankingMassFromNewTag = (double) newTag.CTermFlankingMass -
                                                     GetSequenceMass(newTag.EndIndex, existingTag.EndIndex);

                var toleranceCTerm = existingTag.IsCTermFlankingMassReliable == newTag.IsCTermFlankingMassReliable
                    ? _tolerance
                    : _relaxedTolerance;
                if (!toleranceCTerm.IsWithin(newCTermFlankingMassFromExistingTag, newCTermFlankingMassFromNewTag)) return false;

                if (existingTag.IsCTermFlankingMassReliable || newTag.IsCTermFlankingMassReliable)    // flanking mass from this tag is reliable
                {
                    adjustedCTermFlankingMass =
                        (existingTag.NumReliableCTermFlankingMasses * newCTermFlankingMassFromExistingTag
                        + newTag.NumReliableCTermFlankingMasses * newCTermFlankingMassFromNewTag) /
                        (existingTag.NumReliableCTermFlankingMasses + newTag.NumReliableCTermFlankingMasses);
                }
                else
                {
                    adjustedCTermFlankingMass =
                        (existingTag.NumMergedSequenceTags * newCTermFlankingMassFromExistingTag
                        + newTag.NumMergedSequenceTags * newCTermFlankingMassFromNewTag) /
                        (existingTag.NumMergedSequenceTags + newTag.NumMergedSequenceTags);
                }
            }
            else if (existingTag.CTermFlankingMass != newTag.CTermFlankingMass) return false;

            existingTag.Mass += GetSequenceMass(newStartIndex, existingTag.StartIndex) + GetSequenceMass(existingTag.EndIndex, newEndIndex);
            existingTag.StartIndex = newStartIndex;
            existingTag.EndIndex = newEndIndex;
            existingTag.NTermFlankingMass = adjustedNTermFlankingMass;
            existingTag.CTermFlankingMass = adjustedCTermFlankingMass;
            existingTag.NumMergedSequenceTags += newTag.NumMergedSequenceTags;
            existingTag.NumReliableNTermFlankingMasses += newTag.NumReliableNTermFlankingMasses;
            existingTag.NumReliableCTermFlankingMasses += newTag.NumReliableCTermFlankingMasses;

            return true;
        }

        // startIndex: inclusive, endIndex: exclusive
        // This method is inefficient
        private double GetSequenceMass(int startIndex, int endIndex)
        {
            var mass = 0.0;
            for (var i = startIndex; i < endIndex; i++) mass += (_aminoAcidSet.GetAminoAcid(_sequence[i]) ?? AminoAcid.Empty).Mass;
            return mass;
        }

        private readonly string _sequence;
        private readonly AminoAcidSet _aminoAcidSet;
        private readonly Tolerance _tolerance;
        private readonly Tolerance _relaxedTolerance;
        private readonly List<MatchedTag> _tags;
    }

    public class MatchedTag
    {
        public MatchedTag(SequenceTag tag, int startIndex): this(tag, startIndex, null)
        {
        }

        public MatchedTag(SequenceTag tag, int startIndex, double? featureMass)
        {
            StartIndex = startIndex;
            EndIndex = startIndex + tag.Sequence.Length;
            NumMergedSequenceTags = 1;
            NTermFlankingMass = tag.GetNTermFlankingMass(featureMass - Composition.H2O.Mass);
            CTermFlankingMass = tag.GetCTermFlankingMass(featureMass - Composition.H2O.Mass);

            NumReliableNTermFlankingMasses = tag.IsPrefix ? 1 : 0;
            NumReliableCTermFlankingMasses = tag.IsPrefix ? 0 : 1;
        }

        public MatchedTag(int startIndex, int endIndex, double nTermFlankingMass, double cTermFlankingMass,
            int numMergedSequenceTags, int numReliableNTermFlankingMasses, int numReliableCTermFlankingMasses)
        {
            StartIndex = startIndex;
            EndIndex = endIndex;
            NumMergedSequenceTags = numMergedSequenceTags;
            NTermFlankingMass = nTermFlankingMass;
            CTermFlankingMass = cTermFlankingMass;

            NumReliableNTermFlankingMasses = numReliableNTermFlankingMasses;
            NumReliableCTermFlankingMasses = numReliableCTermFlankingMasses;
        }

        public int StartIndex { get; internal set; }
        public int EndIndex { get; internal set; }    // exclusive
        public int Length { get { return EndIndex - StartIndex; } }
        public double Mass { get; set; }

        public int NumMergedSequenceTags { get; internal set; }
        public double? NTermFlankingMass { get; internal set; }
        public double? CTermFlankingMass { get; internal set; }

        public bool IsNTermFlankingMassReliable { get { return NumReliableNTermFlankingMasses > 0; } }
        public bool IsCTermFlankingMassReliable { get { return NumReliableCTermFlankingMasses > 0; } }

        public int NumReliableNTermFlankingMasses { get; internal set; }
        public int NumReliableCTermFlankingMasses { get; internal set; }

        public MatchedTag Add(MatchedTag tag)
        {
            // N-term
            var newStartIndex = Math.Min(StartIndex, tag.StartIndex);
            if (tag.IsNTermFlankingMassReliable)    // flanking mass from this tag is reliable
            {
                NTermFlankingMass = (NumReliableNTermFlankingMasses * NTermFlankingMass + tag.NTermFlankingMass) /
                                    (NumReliableNTermFlankingMasses + 1);
                ++NumReliableNTermFlankingMasses;
            }
            else  // flanking mass is not reliable
            {
                if (IsNTermFlankingMassReliable)
                {
                    // do nothing
                }
                else
                {
                    NTermFlankingMass = (NumMergedSequenceTags * NTermFlankingMass + tag.NTermFlankingMass) /
                                        (NumMergedSequenceTags + 1);
                }
            }
            StartIndex = newStartIndex;

            // C-term
            var newEndIndex = Math.Max(EndIndex, tag.EndIndex);

            if (tag.IsCTermFlankingMassReliable)    // flanking mass is reliable
            {
                CTermFlankingMass = (NumReliableCTermFlankingMasses * CTermFlankingMass + tag.CTermFlankingMass) /
                                    (NumReliableCTermFlankingMasses + 1);
                ++NumReliableCTermFlankingMasses;
            }
            else  // flanking mass is not reliable
            {
                if (IsCTermFlankingMassReliable)
                {
                    // do nothing
                }
                else
                {
                    CTermFlankingMass = (NumMergedSequenceTags * CTermFlankingMass + tag.CTermFlankingMass) /
                                        (NumMergedSequenceTags + 1);
                }
            }
            EndIndex = newEndIndex;

            ++NumMergedSequenceTags;

            return this;
        }
    }
}
