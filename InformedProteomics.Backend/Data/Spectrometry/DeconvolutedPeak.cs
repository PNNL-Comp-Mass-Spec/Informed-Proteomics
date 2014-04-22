using System;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class DeconvolutedPeak: IComparable<DeconvolutedPeak>
    {
        public DeconvolutedPeak(double mass, double intensity, int charge)
        {
            Mass = mass;
            Intensity = intensity;
            Charge = charge;
        }

        public double Mass { get; private set; }
        public double Intensity { get; private set; }
        public int Charge { get; private set; }

        public int CompareTo(DeconvolutedPeak other)
        {
            return Mass.CompareTo(other.Mass);
        }
    }

}
