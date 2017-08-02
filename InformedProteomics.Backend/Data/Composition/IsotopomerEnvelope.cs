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

        public double[] Envelope { get; }
        public int MostAbundantIsotopeIndex { get; }
    }
}
