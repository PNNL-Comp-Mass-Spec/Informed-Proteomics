using System.Collections.Generic;

namespace InformedProteomics.Backend.Scoring
{
    public class FragmentSpectrum : Dictionary<IonType, double>
    {
        public double GetIntensity(IonType ion)
        {
            if (ion == null) return 0;
            return !ContainsKey(ion) ? 0 : this[ion];
        }

        public void AddIntensity(IonType ion, double intensity)
        {
            var prevIntensity = 0.0;
            if (ContainsKey(ion)) prevIntensity = this[ion];
            this[ion] = prevIntensity + intensity;
        }
    }
}
