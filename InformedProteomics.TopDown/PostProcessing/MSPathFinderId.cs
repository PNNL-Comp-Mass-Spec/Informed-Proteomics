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

        public int Scan { get; private set; }
        public char Pre { get; private set; }
        public string Sequence { get; private set; }
        public char Post { get; private set; }
        public string Modifications { get; private set; }
        public Composition Composition { get; private set; }
        public string ProteinName { get; private set; }
        public string ProteinDesc { get; private set; }
        public int ProteinLength { get; private set; }
        public int Start { get; private set; }
        public int End { get; private set; }
        public int Charge { get; private set; }
        public double MostAbundantIsotopeMz { get; private set; }
        public double Mass { get; private set; }
        public int NumMatchedFragments { get; private set; }
        public double QValue { get; private set; }
        public double PepQValue { get; private set; }
    }
}
