using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.SequenceTag
{
    public class SequenceTagGraphEdge : GraphEdge
    {
        public double Mass { get; private set; }
        public byte NominalMass { get; private set; }
        public List<AminoAcid> AminoAcidList { get; private set; }

        public SequenceTagGraphEdge(int peak1Index, int peak2Index, double peakGapDistance)
            : base(peak1Index, peak2Index)
        {
            Mass = peakGapDistance;
            NominalMass = (byte)Math.Round(Mass * Constants.RescalingConstant);
            AminoAcidList = new List<AminoAcid>();
        }

        public void AddMatchedAminoAcid(AminoAcid aa, double massError)
        {
            AminoAcidList.Add(aa);
        }
    }
}
