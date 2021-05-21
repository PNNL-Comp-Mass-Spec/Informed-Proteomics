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
        // Ignore spelling: proteoform

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="sequence">Proteoform sequence</param>
        /// <param name="pre">Residue in the protein before the sequence</param>
        /// <param name="post">Residue in the protein after the sequence</param>
        /// <param name="scanNum">Scan number</param>
        /// <param name="offset">Offset in the concatenated protein sequence where the protein for this PSM starts</param>
        /// <param name="numNTermCleavages"></param>
        /// <param name="modifications"></param>
        /// <param name="ion"></param>
        /// <param name="score"></param>
        /// <param name="isDecoy"></param>
        /// <param name="specEvalue"></param>
        /// <param name="featureId"></param>
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

        /// <summary>
        /// Proteoform sequence
        /// </summary>
        public string Sequence { get; }

        /// <summary>
        /// Residue in the protein before the sequence
        /// </summary>
        public char Pre { get; }

        /// <summary>
        /// Residue in the protein after the sequence
        /// </summary>
        public char Post { get; }

        /// <summary>
        /// Scan number
        /// </summary>
        public int ScanNum { get; }

        /// <summary>
        /// Offset in the concatenated protein sequence where the protein for this PSM starts
        /// </summary>
        public long Offset { get; }

        public int NumNTermCleavages { get; }
        public ModificationCombination Modifications { get; }
        public Ion Ion { get; }
        public readonly bool IsDecoy;

        public double Score { get; internal set; }
        public double SpecEvalue { get; internal set; }
        public string ModificationText { get; internal set; }

        public int NumMatchedFragments { get; internal set; }

        /// <summary>
        /// MS1 feature ID
        /// </summary>
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
