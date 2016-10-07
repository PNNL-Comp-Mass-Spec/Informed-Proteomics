using System;
using System.Collections.Generic;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.Utils;
using PNNLOmics.Utilities;

namespace InformedProteomics.TopDown.Execution
{
    public class DatabaseSequenceSpectrumMatch: IComparable<DatabaseSequenceSpectrumMatch>
    {
        public DatabaseSequenceSpectrumMatch(string sequence, char pre, char post, int scanNum, long offset,
            int numNTermCleavages, ModificationCombination modifications, Ion ion, double score,
            bool isDecoy,
            double specEvalue = 0.0)
        {
            Sequence = sequence;
            Pre = pre == FastaDatabase.Delimiter ? '-' : pre;
            Post = post == FastaDatabase.Delimiter ? '-' : post;
            ScanNum = scanNum;
            Offset = offset;
            NumNTermCleavages = numNTermCleavages;
            Modifications = modifications;
            Ion = ion;
            Score = score;
            SpecEvalue = specEvalue;
            IsDecoy = isDecoy;
        }

        public string Sequence { get; private set; }
        public char Pre { get; private set; }
        public char Post { get; private set; }
        public int ScanNum { get; private set; }
        public long Offset { get; private set; }
        public int NumNTermCleavages { get; private set; }
        public ModificationCombination Modifications { get; private set; }
        public Ion Ion { get; private set; }
        public readonly bool IsDecoy;

        public double Score { get; internal set; }
        public double SpecEvalue { get; internal set; }
        public string ModificationText { get; internal set; }

        public int NumMatchedFragments { get; internal set; }

        public AminoAcid NTerm
        {
            get { return Pre == '-' ? AminoAcid.ProteinNTerm : AminoAcid.PeptideNTerm; }
        }

        public AminoAcid CTerm
        {
            get { return Post == '-' ? AminoAcid.ProteinCTerm : AminoAcid.PeptideCTerm; }
        }

        public int CompareTo(DatabaseSequenceSpectrumMatch other)
        {
            return Score.CompareTo(other.Score);
        }
    }
}
