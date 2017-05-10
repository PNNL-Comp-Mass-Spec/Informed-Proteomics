namespace InformedProteomics.Backend.Data.Composition
{
    using System;

    [Serializable]
    public class IsotopomerEnvelope
    {
        public IsotopomerEnvelope(double[] envelope, int mostAbundantIsotopeIndex)
        {
            Envelope = envelope;
            MostAbundantIsotopeIndex = mostAbundantIsotopeIndex;
        }

        public double[] Envelope { get; private set; }
        public int MostAbundantIsotopeIndex { get; private set; }
    }
}
