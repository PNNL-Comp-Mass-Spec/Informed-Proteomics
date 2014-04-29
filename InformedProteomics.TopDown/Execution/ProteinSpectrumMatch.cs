using System;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.TopDown.Execution
{
    public class ProteinSpectrumMatch: IComparable<ProteinSpectrumMatch>
    {
        public ProteinSpectrumMatch(string sequence, int scanNum, long offset, int numNTermCleavages, ModificationCombination modifications, Ion ion, double score)
        {
            Sequence = sequence;
            ScanNum = scanNum;
            Offset = offset;
            NumNTermCleavages = numNTermCleavages;
            Modifications = modifications;
            Ion = ion;
            Score = score;
        }

        public string Sequence { get; private set; }
        public int ScanNum { get; private set; }
        public long Offset { get; private set; }
        public int NumNTermCleavages { get; private set; }
        public ModificationCombination Modifications { get; private set; }
        public Ion Ion { get; private set; }
        public double Score { get; private set; }

        public int CompareTo(ProteinSpectrumMatch other)
        {
            return Score.CompareTo(other.Score);
        }
    }
}
