using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Deconvoluted peak class
    /// </summary>
    public class DeconvolutedPeak : Peak
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="mass"></param>
        /// <param name="intensity"></param>
        /// <param name="charge"></param>
        /// <param name="corr"></param>
        /// <param name="dist"></param>
        /// <param name="isotopePeaks"></param>
        public DeconvolutedPeak(double mass, double intensity, int charge, double corr, double dist, Peak[] isotopePeaks = null) : base(mass, intensity)
        {
            Charge = charge;

            Corr = corr;
            Dist = dist;

            ObservedPeaks = isotopePeaks;
            //_isotopePeaks = new HashSet<Peak>();
            /*
            if (isotopePeaks == null)
            {
                SummedIntensity = intensity;
                return;
            }*/

            //SummedIntensity = 0;
            /*
            foreach (var peak in isotopePeaks.Where(peak => peak != null))
            {
                _isotopePeaks.Add(peak);
                SummedIntensity += peak.Intensity;
            }*/
            ObservedPeakIndices = new List<int>();
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="mzPeak"></param>
        /// <param name="charge"></param>
        /// <param name="corr"></param>
        /// <param name="dist"></param>
        public DeconvolutedPeak(Peak mzPeak, int charge, double corr = 0, double dist = 0)
            : this(charge * (mzPeak.Mz - Constants.Proton), mzPeak.Intensity, charge, corr, dist)
        {
        }

        /// <summary>
        /// Get the summed intensity of the deconvoluted peak
        /// </summary>
        public double SummedIntensity { get { return ObservedPeaks?.Where(p => p != null).Sum(p => p.Intensity) ?? Intensity; } }

        /// <summary>
        /// Get the m/z without the adduct ion mass
        /// </summary>
        public double MzWithoutAdductIonMass => Mass / Charge;

        /// <summary>
        /// Get the mass (m/z) of the deconvoluted peak
        /// </summary>
        public double Mass => Mz;

        /// <summary>
        /// Get the charge of the deconvoluted peak
        /// </summary>
        public int Charge { get; }

        /// <summary>
        /// Get the Pearson Correlation of the deconvoluted peak
        /// </summary>
        public double Corr { get; }

        /// <summary>
        /// Get the Bhattacharyya distance of the deconvoluted peak
        /// </summary>
        public double Dist { get; }

        /*
        public bool PeakShare(DeconvolutedPeak other)
        {
            return ObservedPeaks.Any(peak => other.ObservedPeaks.Contains(peak));
        }*/

        /// <summary>
        /// Observed peaks of this deconvoluted peak
        /// </summary>
        public Peak[] ObservedPeaks { get; private set; }
        //private readonly HashSet<Peak> _isotopePeaks;

        /// <summary>
        /// Used for retrieving Observed peaks when reading from a .dpbf file, in conjunction with a .pbf file.
        /// </summary>
        public List<int> ObservedPeakIndices { get; }

        /// <summary>
        /// Uses the peaks in spec to find and set the <see cref="ObservedPeaks"/> according to the indices in <see cref="ObservedPeakIndices"/>
        /// </summary>
        /// <param name="spec"></param>
        public void SetObservedPeaksFromSpectrum(Spectrum spec)
        {
            var peaks = new List<Peak>();
            foreach (var index in ObservedPeakIndices)
            {
                peaks.Add(spec.Peaks[index]);
            }
            ObservedPeaks = peaks.ToArray();
        }
    }

    /*
    public class DeconvolutedPeak: IComparable<DeconvolutedPeak>
    {
        public DeconvolutedPeak(double mass, int charge, Peak[] isotopePeaks = null, double[] envelope = null, double relativeThreshold = 0.3)
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
                    if (envelope[i] > relativeThreshold) IsotopePeaks.Add(isotopePeaks[i]);
                }
            }
        }

        public double Mass { get; private set; }
        public double Corr { get; internal set; }
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
            if (other == null) return false;
            if (ReferenceEquals(this, other)) return true;
            return Math.Abs(Mass - other.Mass) < 1e-9;
        }

        public override bool Equals(object obj)
        {
            if (obj == null) return false;
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
    }*/
}
