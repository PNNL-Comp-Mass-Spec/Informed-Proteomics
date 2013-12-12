using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.DIA.Scoring
{
    public class FragmentOffsetFrequency
    {
        public FragmentOffsetFrequency(IonType ionType, float frequency)
        {
            IonType = ionType;
            Frequency = frequency;
        }

        public IonType IonType { get; private set; }
        public float Frequency { get; private set; }
    }
}
