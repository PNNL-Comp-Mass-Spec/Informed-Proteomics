using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class DeconvolutedPeak: IComparable<DeconvolutedPeak>
    {
        public DeconvolutedPeak(double mass, int charge, Peak[] isotopePeaks = null, double[] envelop = null, double relativeThreshold = 0.3)
        {
            Mass = mass;
            Charge = charge;
            Intensity = 0d;
            IsotopePeaks = new HashSet<Peak>();

            if (isotopePeaks != null)
            {
                for(var i = 0; i < isotopePeaks.Length; i++)
                {
                    if (isotopePeaks[i] == null) continue;
                    Intensity += isotopePeaks[i].Intensity;
                    if (envelop[i] > relativeThreshold) IsotopePeaks.Add(isotopePeaks[i]);
                }
            }
        }

        public double Mass { get; private set; }
        
        //public double Intensity { get { return MzPeak != null ? MzPeak.Intensity : 0; } }
        //public double Mz { get { return MzPeak.Mz;  } }
        //public Peak MzPeak { get; private set; }
        public double Intensity { get; private set; }

        public bool PeakShare(DeconvolutedPeak other)
        {
            return IsotopePeaks.Any(peak => other.IsotopePeaks.Contains(peak));
        }

        internal readonly HashSet<Peak> IsotopePeaks;
        public int Charge { get; private set; }

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
