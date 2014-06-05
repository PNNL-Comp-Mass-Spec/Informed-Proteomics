using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDownViewer
{
    public class ChargeStateData
    {
        public int Scan { get; private set; }
        public ProductSpectrum Ms2Spectrum { get; private set; }
        public string Sequence { get; private set; }
        public string Pre { get; private set; }
        public string Post { get; private set; }
        public string Modifications { get; private set; }
        public string Composition { get; private set; }
        public string ProteinName { get; private set; }
        public string ProteinDesc { get; private set; }
        public int ProteinLength { get; private set; }
        public int Start { get; private set; }
        public int End { get; private set; }
        public int Charge { get; private set; }
        public double MostAbundantIsotopeMz { get; private set; }
        public double Mass { get; private set; }
        public int MatchedFragments { get; private set; }
        public double IsotopeCorrPrevMs1 { get; private set; }
        public double IsotopeCorrNextMs1 { get; private set; }
        public double CorrMostAbundantPlusOneIsoptope { get; private set; }
        public double ChargeCorrMinusOne { get; private set; }
        public double ChargeCorrPlusOne { get; private set; }
        public double QValue { get; private set; }
        public double PepQValue { get; private set; }

        public ChargeStateData(string line, LcMsRun lcms)
        {
            var parts = line.Split('\t');
            Scan = Convert.ToInt32(parts[0]);
            Ms2Spectrum = lcms.GetSpectrum(Scan) as ProductSpectrum;
            Pre = parts[1];
            Sequence = parts[2];
            Post = parts[3];
            Modifications = parts[4];
            Composition = parts[5];
            ProteinName = parts[6];
            ProteinDesc = parts[7];
            ProteinLength = Convert.ToInt32(parts[8]);
            Start = Convert.ToInt32(parts[9]);
            End = Convert.ToInt32(parts[10]);
            Charge = Convert.ToInt32(parts[11]);
            MostAbundantIsotopeMz = Convert.ToDouble(parts[12]);
            Mass = Convert.ToDouble(parts[13]);
            MatchedFragments = Convert.ToInt32(parts[14]);
            IsotopeCorrPrevMs1 = Convert.ToDouble(parts[15]);
            IsotopeCorrNextMs1 = Convert.ToDouble(parts[16]);
            CorrMostAbundantPlusOneIsoptope = Convert.ToDouble(parts[17]);
            ChargeCorrMinusOne = Convert.ToDouble(parts[18]);
            ChargeCorrPlusOne = Convert.ToDouble(parts[19]);
            QValue = Convert.ToDouble(parts[20]);
            PepQValue = Convert.ToDouble(parts[21]);
        }
    }
}
