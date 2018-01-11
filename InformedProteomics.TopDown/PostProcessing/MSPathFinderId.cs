using InformedProteomics.Backend.Data.Composition;

namespace InformedProteomics.TopDown.PostProcessing
{
    public class MsPathFinderId
    {
        public MsPathFinderId(int scan, char pre, string sequence, char post, string modifications,
            Composition composition, string proteinName, string proteinDesc, int proteinLength,
            int start, int end, int charge, double mostAbundantIsotopeMz, double mass,
            int numMatchedFragments, double qValue, double pepQValue)
        {
            Scan = scan;
            Pre = pre;
            Sequence = sequence;
            Post = post;
            Modifications = modifications;
            Composition = composition;
            ProteinName = proteinName;
            ProteinDesc = proteinDesc;
            ProteinLength = proteinLength;
            Start = start;
            End = end;
            Charge = charge;
            MostAbundantIsotopeMz = mostAbundantIsotopeMz;
            Mass = mass;
            NumMatchedFragments = numMatchedFragments;
            QValue = qValue;
            PepQValue = pepQValue;
        }

        public int Scan { get; }
        public char Pre { get; }
        public string Sequence { get; }
        public char Post { get; }
        public string Modifications { get; }
        public Composition Composition { get; }
        public string ProteinName { get; }
        public string ProteinDesc { get; }
        public int ProteinLength { get; }
        public int Start { get; }
        public int End { get; }
        public int Charge { get; }
        public double MostAbundantIsotopeMz { get; }
        public double Mass { get; }
        public int NumMatchedFragments { get; }
        public double QValue { get; }
        public double PepQValue { get; }
    }
}
