using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.FeatureFinding.SpectrumMatching
{
    public class ProteinSpectrumMatchSet : List<ProteinSpectrumMatch>
    {
        public ProteinSpectrumMatchSet(int dataid)
        {
            DataId = dataid;
        }
        public ProteinSpectrumMatchSet(int dataid, IEnumerable<ProteinSpectrumMatch> groupedMatches)
        {
            DataId = dataid;
            AddRange(groupedMatches);
        }

        public void Merge(ProteinSpectrumMatchSet other)
        {
            AddRange(other);
        }

        public int MinScanNum { get { return this.Min(item => item.ScanNum); } }
        public int MaxScanNum { get { return this.Max(item => item.ScanNum); } }

        public int MinCharge { get { return this.Min(item => item.Charge); } }
        public int MaxCharge { get { return this.Max(item => item.Charge); } }

        public double Mass { get { return this.Select(item => item.Mass).Median(); } }

        public void SetDataId(int dataid)
        {
            DataId = dataid;
        }

        public int DataId { get; private set; }

        public bool ShareProteinId(ProteinSpectrumMatchSet other)
        {
            foreach (var prsm1 in this)
            {
                foreach (var prsm2 in other)
                {
                    if (prsm1.ProteinId.Equals(prsm2.ProteinId)) return true;
                }
            }
            return false;
        }
    }

    public class ProteinSpectrumMatch : IEquatable<ProteinSpectrumMatch>, IComparable<ProteinSpectrumMatch>
    {
        public ProteinSpectrumMatch(string sequence, int scanNum, double mass, int charge, string protName, string protDesc, int firstResidue, int lastResidue, double score, SearchTool searchTool = SearchTool.Unknown)
        {
            Sequence = sequence;
            ScanNum = scanNum;
            Mass = mass;
            Charge = charge;

            ProteinName = protName;

            FirstResidue = firstResidue;
            LastResidue = lastResidue;
            SearchToolType = searchTool;
            Score = score;
            ProteinDesc = protDesc;
        }

        public double SpectralEvalue { get; internal set; }

        public string Sequence { get; private set; }
        public int ScanNum { get; private set; }
        public double Mass { get; private set; }
        public int Charge { get; private set; }

        public string ProteinName { get; private set; }
        public string ProteinDesc { get; private set; }

        public int ProteinLength { get; internal set; }

        public string Pre { get; internal set; }
        public string Post { get; internal set; }

        public int FirstResidue { get; private set; }
        public int LastResidue { get; private set; }
        public double Score { get; private set; }
        public string SequenceText { get; internal set; }

        public string Modifications { get; internal set; }

        public SearchTool SearchToolType { get; private set; }
        public string ProteinId { get; set; }

        public bool Equals(ProteinSpectrumMatch other)
        {
            if (SearchToolType == other.SearchToolType)
            {
                return SequenceText.Equals(other.SequenceText);
            }

            var massDiff = Math.Abs(Mass - other.Mass);
            var tol = new Tolerance(10);
            if (massDiff < tol.GetToleranceAsMz(Mass) && FirstResidue == other.FirstResidue && LastResidue == other.LastResidue) return true;

            return false;
        }

        public enum SearchTool
        {
            MsAlign,
            MsPathFinder,
            MsGfPlus,
            Unknown,
        }

        public int CompareTo(ProteinSpectrumMatch other)
        {
            return other.Score.CompareTo(Score);
        }

        public Sequence GetSequence()
        {
            if (SearchToolType == SearchTool.MsGfPlus)
            {
                return Backend.Data.Sequence.Sequence.GetSequenceFromMsGfPlusPeptideStr(Sequence);
            }
            else if (SearchToolType == SearchTool.MsPathFinder)
            {
                return Backend.Data.Sequence.Sequence.CreateSequence(Sequence, Modifications, new AminoAcidSet());
            }
            // todo : MsAlign

            return null;
        }
    }
}
