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
        public List<KeyValuePair<AminoAcid, double>> AminoAcidList { get; private set; }

        public SequenceTagGraphEdge(int peak1Index, int peak2Index, double peakGapDistance)
            : base(peak1Index, peak2Index)
        {
            Mass = peakGapDistance;
            NominalMass = (byte)Math.Round(Mass * Constants.RescalingConstant);
            AminoAcidList = new List<KeyValuePair<AminoAcid, double>>();
        }

        public void AddMatchedAminoAcid(AminoAcid aa, double massError)
        {
            AminoAcidList.Add(new KeyValuePair<AminoAcid, double>(aa, massError));
        }
    }
}
