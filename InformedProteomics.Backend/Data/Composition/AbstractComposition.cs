using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Composition
{
    /// <summary>
    /// Composition base class
    /// </summary>
    public abstract class AbstractComposition
    {
        /// <summary>
        /// Composition mass
        /// </summary>
        public abstract double Mass { get; }

        /// <summary>
        /// Composition nominal mass
        /// </summary>
        public abstract int NominalMass { get; }

        #region Methods to get masses

        /// <summary>
        /// Gets the mass of nth isotope
        /// </summary>
        /// <param name="isotopeIndex">isotope index. 0 means mono-isotope, 1 means 2nd isotope, etc.</param>
        /// <returns>m/z</returns>
        public double GetIsotopeMass(int isotopeIndex)
        {
            return Mass + isotopeIndex * Constants.C13MinusC12;
        }

        /// <summary>
        /// Gets the m/z of nth isotope
        /// </summary>
        /// <param name="isotopeIndexInRealNumber">isotope index in real number. 0 means mono-isotope, 0.5 means the center of mono and 2nd isotopes</param>
        /// <returns>m/z</returns>
        public double GetIsotopeMass(double isotopeIndexInRealNumber)
        {
            return Mass + isotopeIndexInRealNumber * Constants.C13MinusC12;
        }

        #endregion

        #region Methods to get Isotopomer Envelopes

        /// <summary>
        /// Get the <see cref="IsotopomerEnvelope"/> for this composition
        /// </summary>
        /// <returns>Isotopomer envelope</returns>
        public IsotopomerEnvelope GetIsotopomerEnvelope()
        {
            return _isotopomerEnvelope ??
                   (_isotopomerEnvelope = Averagine.GetIsotopomerEnvelopeFromNominalMass(NominalMass));
        }

        /// <summary>
        /// Get the relative intensities of the Isotopomer envelope for this composition
        /// </summary>
        /// <returns>Array of intensities</returns>
        public double[] GetIsotopomerEnvelopeRelativeIntensities()
        {
            if (_isotopomerEnvelope == null)
            {
                _isotopomerEnvelope = Averagine.GetIsotopomerEnvelopeFromNominalMass(NominalMass);
            }

            return _isotopomerEnvelope.Envelope;
        }

        /// <summary>
        /// Get the zero-based index of the most abundant isotope, according to the isotopomer envelope
        /// </summary>
        /// <returns>Index of the most abundant isotope</returns>
        public int GetMostAbundantIsotopeZeroBasedIndex()
        {
            if (_isotopomerEnvelope == null)
            {
                _isotopomerEnvelope = Averagine.GetIsotopomerEnvelopeFromNominalMass(NominalMass);
            }

            return _isotopomerEnvelope.MostAbundantIsotopeIndex;
        }

        private IsotopomerEnvelope _isotopomerEnvelope;

        #endregion
    }
}
