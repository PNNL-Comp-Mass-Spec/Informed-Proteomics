using System;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.FeatureFinding.Data;

namespace InformedProteomics.FeatureFinding.IsotopicEnvelope
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

        public double BhattacharyyaDistance => TheoreticalEnvelope.GetBhattacharyyaDistance(Peaks);

        public double PearsonCorrelation => TheoreticalEnvelope.GetPearsonCorrelation(Peaks);

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
                for (var i = 0; i < Size; i++)
                {
                    if (Peaks[i]?.Active == true)
                    {
                        prob[i] = Peaks[i].Intensity / s;
                    }
                }

                return prob;
            }
        }

        public Ms1Peak[] Peaks { get; }
        public Ms1Peak MinMzPeak { get { return Peaks.FirstOrDefault(t => t?.Active == true); } }
        public Ms1Peak MaxMzPeak { get { return Peaks.LastOrDefault(t => t?.Active == true); } }
        public int NumberOfPeaks { get { return Peaks.Count(x => x?.Active == true); } }
        public double Abundance { get { return Peaks.Where(p => p?.Active == true).Sum(p => p.Intensity); } }

        public Ms1Peak RepresentativePeak
        {
            get
            {
                foreach (var i in TheoreticalEnvelope.IndexOrderByRanking)
                {
                    if (Peaks[i]?.Active == true)
                    {
                        return Peaks[i];
                    }
                }

                return null;
            }
        }

        public double HighestIntensity
        {
            get
            {
                var ret = 0d;
                foreach (var peak in Peaks.Where(p => p?.Active == true))
                {
                    if (peak.Intensity > ret)
                    {
                        ret = peak.Intensity;
                    }
                }

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
