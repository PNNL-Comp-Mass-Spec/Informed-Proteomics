using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Biology
{
    public class Ion
    {
        public Composition Composition { get; private set; }
        public int Charge { get; private set; }

        public Ion(Composition composition, int charge)
        {
            Composition = composition;
            Charge = charge;
        }

        public double GetMz()
        {
            return (Composition.GetMass() + Charge * Constants.H) / Charge;
        }

        /// <summary>
        /// Gets the m/z of ith isotope
        /// </summary>
        /// <param name="isotopeIndex">isotope index. 0 means mono-isotope, 1 means 2nd isotope, etc.</param>
        /// <returns></returns>
        public double GetIsotopeMz(int isotopeIndex)
        {
            return (Composition.GetIsotopeMass(isotopeIndex) + Charge * Constants.H) / Charge;
        }
    }
}
