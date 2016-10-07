using System;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.FeatureFinding.MassFeature
{
    public abstract class IsotopeEnvelope
    {
        public abstract Isotope GetMostAbundantIsotope();
        public abstract double[] Probability { get; }

        public double MonoMass { get; protected set;  }
        public int Size { get { return Probability.Length; } }

        public const double MaxBhattacharyyaDistance = 10.0d;

        public double GetBhattacharyyaDistance(IsotopeEnvelope other)
        {
            var bc = 0d;
            for (var i = 0; i < Size; i++)
            {
                var p = Probability[i];
                var q = other.Probability[i];
                bc += Math.Sqrt(p * q);
            }

            if (!(bc > 0)) return MaxBhattacharyyaDistance;

            return -Math.Log(bc);
        }

        public double GetBhattacharyyaDistance(Ms1Peak[] isotopePeaks)
        {
            var s2 = 0d;

            for (var i = 0; i < Size; i++)
            {
                if (isotopePeaks[i] != null && isotopePeaks[i].Active)
                    s2 += isotopePeaks[i].Intensity;
            }

            if (!(s2 > 0)) return MaxBhattacharyyaDistance;

            var bc = 0d;
            for (var i = 0; i < Size; i++)
            {
                var p = Probability[i];
                var q = (isotopePeaks[i] != null && isotopePeaks[i].Active) ? isotopePeaks[i].Intensity / s2 : 0;
                bc += Math.Sqrt(p * q);
            }

            if (!(bc > 0)) return MaxBhattacharyyaDistance;

            return -Math.Log(bc);
        }

        public double GetBhattacharyyaDistance(double[] intensities)
        {
            var s2 = 0d;

            for (var i = 0; i < Size; i++)
            {
                s2 += intensities[i];
            }

            if (!(s2 > 0)) return MaxBhattacharyyaDistance;

            var bc = 0d;
            for (var i = 0; i < Size; i++)
            {
                var p = Probability[i];
                var q = intensities[i] / s2;
                bc += Math.Sqrt(p * q);
            }
            if (!(bc > 0)) return MaxBhattacharyyaDistance;

            return -Math.Log(bc);
        }

        public double GetPearsonCorrelation(IsotopeEnvelope other)
        {
            var m = 1 / Size;
            var cov = 0.0;
            var s1 = 0.0;
            var s2 = 0.0;

            for (var i = 0; i < Size; i++)
            {
                var d1 = Probability[i] - m;
                var d2 = other.Probability[i] - m;
                cov += d1 * d2;
                s1 += d1 * d1;
                s2 += d2 * d2;
            }

            if (s1 <= 0 || s2 <= 0) return 0;

            return cov < 0 ? 0d : cov / Math.Sqrt(s1 * s2);
        }

        public double GetPearsonCorrelation(Ms1Peak[] isotopePeaks)
        {
            var m1 = 1.0/Size;
            var m2 = 0.0;

            for (var i = 0; i < Size; i++)
            {
                if (isotopePeaks[i] != null && isotopePeaks[i].Active) m2 += isotopePeaks[i].Intensity;
            }
            m2 /= Size;

            // compute Pearson correlation
            var cov = 0.0;
            var s1 = 0.0;
            var s2 = 0.0;

            for (var i = 0; i < Size; i++)
            {
                var d1 = Probability[i] - m1;
                var d2 = (isotopePeaks[i] != null && isotopePeaks[i].Active) ? isotopePeaks[i].Intensity - m2 : -m2;
                cov += d1 * d2;
                s1 += d1 * d1;
                s2 += d2 * d2;
            }

            if (s1 <= 0 || s2 <= 0) return 0;

            return cov < 0 ? 0d : cov / Math.Sqrt(s1 * s2);
        }

        public double GetPearsonCorrelation(double[] intensities)
        {
            var m1 = 1.0 / Size;
            var m2 = 0.0;

            for (var i = 0; i < Size; i++)
            {
                m2 += intensities[i];
            }
            m2 /= Size;

            // compute Pearson correlation
            var cov = 0.0;
            var s1 = 0.0;
            var s2 = 0.0;

            for (var i = 0; i < Size; i++)
            {
                var d1 = Probability[i] - m1;
                var d2 = intensities[i] - m2;
                cov += d1 * d2;
                s1 += d1 * d1;
                s2 += d2 * d2;
            }

            if (s1 <= 0 || s2 <= 0) return 0;

            return cov < 0 ? 0d : cov / Math.Sqrt(s1 * s2);
        }
    }
}
