using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Composition
{
    using System;

    [Serializable]
    public abstract class AbstractComposition
    {
        public abstract double Mass { get; }
        public abstract int NominalMass { get; }

        #region Methods to get masses

        /// <summary>
        /// Gets the mass of ith isotope
        /// </summary>
        /// <param name="isotopeIndex">isotope index. 0 means mono-isotope, 1 means 2nd isotope, etc.</param>
        /// <returns></returns>
        public double GetIsotopeMass(int isotopeIndex)
        {
            return Mass + isotopeIndex * Constants.C13MinusC12;
        }

        /// <summary>
        /// Gets the m/z of ith isotope
        /// </summary>
        /// <param name="isotopeIndexInRealNumber">isotope index in real number. 0 means mono-isotope, 0.5 means the center of mono and 2nd isotopes.</param>
        /// <returns></returns>
        public double GetIsotopeMass(double isotopeIndexInRealNumber)
        {
            return Mass + isotopeIndexInRealNumber * Constants.C13MinusC12;
        }

        #endregion

        #region Methods to get Isotopomer Envelpops

        public IsotopomerEnvelope GetIsotopomerEnvelope()
        {
            return _isotopomerEnvelop ??
                   (_isotopomerEnvelop = Averagine.GetIsotopomerEnvelopeFromNominalMass(NominalMass));
        }

        public double[] GetIsotopomerEnvelopeRelativeIntensities()
        {
            if (_isotopomerEnvelop == null) _isotopomerEnvelop = Averagine.GetIsotopomerEnvelopeFromNominalMass(NominalMass);
            return _isotopomerEnvelop.Envolope;
        }

        public int GetMostAbundantIsotopeZeroBasedIndex()
        {
            if (_isotopomerEnvelop == null) _isotopomerEnvelop = Averagine.GetIsotopomerEnvelopeFromNominalMass(NominalMass);
            return _isotopomerEnvelop.MostAbundantIsotopeIndex;
        }

        private IsotopomerEnvelope _isotopomerEnvelop;

        #endregion
    }
}
