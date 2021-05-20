using System;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Database;

namespace InformedProteomics.TopDown.Execution
{
    /// <summary>
    /// This class tracks a single proteoform match to a mass spectrum
    /// </summary>
    public class DatabaseSequenceSpectrumMatch : IComparable<DatabaseSequenceSpectrumMatch>
    {
        public DatabaseSequenceSpectrumMatch(string sequence, char pre, char post, int scanNum, long offset,
            int numNTermCleavages, ModificationCombination modifications, Ion ion, double score,
            bool isDecoy,
            double specEvalue = 0.0, int featureId = 0)
        {
            Sequence = sequence;
            Pre = pre == FastaDatabaseConstants.Delimiter ? '-' : pre;
            Post = post == FastaDatabaseConstants.Delimiter ? '-' : post;
            ScanNum = scanNum;
            Offset = offset;
            NumNTermCleavages = numNTermCleavages;
            Modifications = modifications;
            Ion = ion;
            Score = score;
            SpecEvalue = specEvalue;
            IsDecoy = isDecoy;
            FeatureId = featureId;
        }

        public Sequence IpSequence { get; set; }
        public string Sequence { get; }
        public char Pre { get; }
        public char Post { get; }
        public int ScanNum { get; }
        public long Offset { get; }
        public int NumNTermCleavages { get; }
        public ModificationCombination Modifications { get; }
        public Ion Ion { get; }
        public readonly bool IsDecoy;

        public double Score { get; internal set; }
        public double SpecEvalue { get; internal set; }
        public string ModificationText { get; internal set; }

        public int NumMatchedFragments { get; internal set; }

        public int FeatureId { get; private set; }

        public AminoAcid NTerm => Pre == '-' ? AminoAcid.ProteinNTerm : AminoAcid.PeptideNTerm;

        public AminoAcid CTerm => Post == '-' ? AminoAcid.ProteinCTerm : AminoAcid.PeptideCTerm;

        public int CompareTo(DatabaseSequenceSpectrumMatch other)
        {
            return Score.CompareTo(other.Score);
        }

        public void UpdateFeatureId(int ms1FeatureId)
        {
            FeatureId = ms1FeatureId;
        }
    }
}
