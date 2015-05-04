using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class DeconvolutedPeak: IComparable<DeconvolutedPeak>
    {
        public DeconvolutedPeak(Peak mzPeak, double mass, int charge)
        {
            Mass = mass;
            //Intensity = intensity;
            Charge = charge;
            MzPeak = mzPeak;
        }

        public double Mass { get; private set; }
        //public double Intensity { get; private set; }
        public double Intensity { get { return MzPeak != null ? MzPeak.Intensity : 0; } }
        public double Mz { get { return MzPeak.Mz;  } }
        public int Charge { get; private set; }
        public Peak MzPeak { get; private set; }

        public int CompareTo(DeconvolutedPeak other)
        {
            return Mass.CompareTo(other.Mass);
        }

        public bool Equals(DeconvolutedPeak other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return Math.Abs(Mass - other.Mass) < 1e-9;
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != typeof(DeconvolutedPeak)) return false;
            return Equals((DeconvolutedPeak)obj);
        }

        public override int GetHashCode()
        {
            return Mass.GetHashCode();
        }

        public static bool operator ==(DeconvolutedPeak left, DeconvolutedPeak right)
        {
            return Equals(left, right);
        }

        public static bool operator !=(DeconvolutedPeak left, DeconvolutedPeak right)
        {
            return !Equals(left, right);
        }
    }


}
