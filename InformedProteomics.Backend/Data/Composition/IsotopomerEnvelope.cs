namespace InformedProteomics.Backend.Data.Composition
{
    using System;

    [Serializable]
    public class IsotopomerEnvelope
    {
        public IsotopomerEnvelope(double[] envolope, int mostAbundantIsotopeIndex)
        {
            Envolope = envolope;
            MostAbundantIsotopeIndex = mostAbundantIsotopeIndex;
        }

        public double[] Envolope { get; private set; }
        public int MostAbundantIsotopeIndex { get; private set; }
    }
}
