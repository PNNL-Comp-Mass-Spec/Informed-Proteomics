using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.MassFeature
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
        public readonly int DataId;

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
        public ProteinSpectrumMatch(string sequence, int scanNum, double mass, int charge, string protName, int firstResidue, int lastResidue, double score = 0.0)
        {
            Sequence = sequence;
            ScanNum = scanNum;
            Mass = mass;
            Charge = charge;

            ProteinName = protName;

            FirstResidue = firstResidue;
            LastResidue = lastResidue;
            Score = score;
        }

        public string Sequence { get; private set; }
        public int ScanNum { get; private set; }
        public double Mass { get; private set; }
        public int Charge { get; private set; }

        public string ProteinName { get; private set; }
        public int FirstResidue { get; private set; }
        public int LastResidue { get; private set; }
        public double Score { get; private set; }
        public string SequenceText { get; internal set; }
        public SearchTool SearchToolType { get; internal set; }

        public string ProteinId { get; set; }

        //public LcMsFeature LcMsFeature { get; set; }

        public bool Equals(ProteinSpectrumMatch other)
        {
            if (SearchToolType == other.SearchToolType)
            {
                return SequenceText.Equals(other.SequenceText);
            }
            
            var massDiff = Math.Abs(Mass - other.Mass);
            var tol = new Tolerance(10);                
            if (massDiff < tol.GetToleranceAsTh(Mass) && FirstResidue == other.FirstResidue && LastResidue == other.LastResidue) return true;
            
            return false;
        }

        public enum SearchTool
        {
            MsAlign,
            MsPathFinder,
        }

        public int CompareTo(ProteinSpectrumMatch other)
        {
            return other.Score.CompareTo(Score);
        }
    }
}
