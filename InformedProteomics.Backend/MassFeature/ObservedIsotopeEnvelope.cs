using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassFeature
{
    public class ObservedIsotopeEnvelope : IsotopeEnvelope
    {
        public ObservedIsotopeEnvelope(double monoMass, int charge, int scanNum, Ms1Peak[] peaks, TheoreticalIsotopeEnvelope theoreticalEnvelope)
        {
            MonoMass = monoMass;
            Charge = charge;
            ScanNum = scanNum;
            TheoreticalEnvelope = theoreticalEnvelope;
            Peaks = new Ms1Peak[theoreticalEnvelope.Size];

            Array.Copy(peaks, Peaks, theoreticalEnvelope.Size);
            GoodEnough = false;
        }

        public double BhattacharyyaDistance
        {
            get { return TheoreticalEnvelope.GetBhattacharyyaDistance(Peaks); }
        }

        public double PearsonCorrelation
        {
            get { return TheoreticalEnvelope.GetPearsonCorrelation(Peaks); }
        }

        public override Isotope GetMostAbundantIsotope()
        {
            return TheoreticalEnvelope.GetMostAbundantIsotope();
        }

        public override double[] Probability
        {
            get
            {
                var s = Abundance;
                var prob = new double[Size];
                for(var i = 0; i < Size; i++) if (Peaks[i] != null && Peaks[i].Active) prob[i] = Peaks[i].Intensity/s;
                return prob;
            }
        }

        public Ms1Peak[] Peaks { get; private set; }
        public Ms1Peak MinMzPeak { get { return Peaks.FirstOrDefault(t => t != null && t.Active); } }
        public Ms1Peak MaxMzPeak { get { return Peaks.LastOrDefault(t => t != null && t.Active); } }
        public int NumberOfPeaks { get { return Peaks.Count(x => x != null && x.Active); } }
        public double Abundance { get { return Peaks.Where(p => p != null && p.Active).Sum(p => p.Intensity); } }

        public Ms1Peak RepresentativePeak
        {
            get
            {
                foreach (var i in TheoreticalEnvelope.IndexOrderByRanking) if (Peaks[i] != null && Peaks[i].Active) return Peaks[i];
                return null;
            }
        }

        public double HighestIntensity
        {
            get
            {
                var ret = 0d;
                foreach (var peak in Peaks.Where(p => p != null && p.Active)) if (peak.Intensity > ret) ret = peak.Intensity;
                return ret;
            }
        }

        public void SumEnvelopeTo(double[] targetEnvelope)
        {
            Peaks.SumEnvelopeTo(targetEnvelope);
        }

        public readonly int Charge;
        public readonly int ScanNum;
        public readonly TheoreticalIsotopeEnvelope TheoreticalEnvelope;

        public bool GoodEnough { get; internal set; }
    }
}
