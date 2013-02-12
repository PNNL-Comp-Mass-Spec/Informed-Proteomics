using System.Collections.Generic;

namespace InformedProteomics.Backend.Scoring
{
    internal class FragmentXICSet : Dictionary<IonType, double[]>
    {
        public double[] GetIonXIC(IonType ion)
        {
            return !ContainsKey(ion) ? null : this[ion];
        }

        public void AddIonXIC(IonType ion, double[] ionXIC)
        {
            var prevIonXIC = ContainsKey(ion) ? this[ion] : new double[ionXIC.Length];
            var updatedIonXIC = new double[ionXIC.Length];
            for (int i = 0; i < ionXIC.Length; i++)
            {
                updatedIonXIC[i] = ionXIC[i] + prevIonXIC[i];
            }
            this[ion] = updatedIonXIC;
        }
    }
}
