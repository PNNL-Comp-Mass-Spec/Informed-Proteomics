namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Like <see cref="Peak"/>, but with additional per-peak data, like the Noise value available from Thermo Label data.
    /// </summary>
    public class PeakDetailed : Peak
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="mz"></param>
        /// <param name="intensity"></param>
        /// <param name="noise"></param>
        public PeakDetailed(double mz, double intensity, double noise) : base(mz, intensity)
        {
            Noise = noise;
        }

        /// <summary>
        /// Noise value for the peak
        /// </summary>
        public double Noise { get; }

        /// <summary>
        /// Calculated SignalToNoise value
        /// </summary>
        public double SignalToNoise => Intensity / (!Noise.Equals(0) ? Noise : double.Epsilon);
    }
}
