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
            return (Composition.GetMass() + Constants.H2O + Charge * Constants.H) / Charge;
        }
    }
}
