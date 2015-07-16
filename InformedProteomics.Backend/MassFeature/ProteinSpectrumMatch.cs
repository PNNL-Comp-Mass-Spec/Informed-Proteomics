using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.MassFeature
{
    public class ProteinSpectrumMatch : IEquatable<ProteinSpectrumMatch>, IComparable<ProteinSpectrumMatch>
    {
        public ProteinSpectrumMatch(string sequence, int scanNum, double mass, 
            string protName, 
            int firstResidue, int lastResidue, double score = 0.0)
        {
            Sequence = sequence;
            ScanNum = scanNum;
            Mass = mass;

            ProteinName = protName;

            FirstResidue = firstResidue;
            LastResidue = lastResidue;
            Score = score;
        }

        public string Sequence { get; private set; }
        public int ScanNum { get; private set; }
        public double Mass { get; private set; }

        public string ProteinName { get; private set; }
        public int FirstResidue { get; private set; }
        public int LastResidue { get; private set; }
        public double Score { get; private set; }
        public string SequenceText { get; internal set; }
        public SearchTool SearchToolType { get; internal set; }

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
