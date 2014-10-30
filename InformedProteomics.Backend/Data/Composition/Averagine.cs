using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Composition
{
    public class Averagine
    {
        public static IsotopomerEnvelope GetIsotopomerEnvelope(double monoIsotopeMass)
        {
            var nominalMass = (int) Math.Round(monoIsotopeMass*Constants.RescalingConstant);
            return GetIsotopomerEnvelopeFromNominalMass(nominalMass);
        }

        public static List<Peak> GetTheoreticalIsotopeProfile(double monoIsotopeMass, int charge, double relativeIntensityThreshold = 0.1)
        {
            var peakList = new List<Peak>();
            var envelope = GetIsotopomerEnvelope(monoIsotopeMass);
            for (var isotopeIndex = 0; isotopeIndex < envelope.Envolope.Length; isotopeIndex++)
            {
                var intensity = envelope.Envolope[isotopeIndex];
                if (intensity < relativeIntensityThreshold) continue;
                var mz = Ion.GetIsotopeMz(monoIsotopeMass, charge, isotopeIndex);
                peakList.Add(new Peak(mz, intensity));
            }
            return peakList;
        }

        public static IsotopomerEnvelope GetIsotopomerEnvelopeFromNominalMass(int nominalMass)
        {
            IsotopomerEnvelope envelope;
            if (IsotopeEnvelopMap.TryGetValue(nominalMass, out envelope)) return envelope;

            var mass = nominalMass/Constants.RescalingConstant;
            envelope = ComputeIsotopomerEnvelope(mass);
            IsotopeEnvelopMap[nominalMass] = envelope;

            return envelope;
        }

        private const double C = 4.9384;
        private const double H = 7.7583;
        private const double N = 1.3577;
        private const double O = 1.4773;
        private const double S = 0.0417;
        private const double AveragineMass = C * Atom.C + H * Atom.H + N * Atom.N + O * Atom.O + S * Atom.S;

        private static readonly Dictionary<int, IsotopomerEnvelope> IsotopeEnvelopMap; // NominalMass -> Isotope Envelop

        static Averagine()
        {
            IsotopeEnvelopMap = new Dictionary<int, IsotopomerEnvelope>();
        }

        private static IsotopomerEnvelope ComputeIsotopomerEnvelope(double mass)
        {
            var numAveragines = mass / AveragineMass;
            var numC = (int)Math.Round(C * numAveragines);
            var numH = (int)Math.Round(H * numAveragines);
            var numN = (int)Math.Round(N * numAveragines);
            var numO = (int)Math.Round(O * numAveragines);
            var numS = (int)Math.Round(S * numAveragines);

            if (numH == 0) numH = 1;
            return IsoProfilePredictor.GetIsotopomerEnvelop(numC, numH, numN, numO, numS);
        }
    }
}
