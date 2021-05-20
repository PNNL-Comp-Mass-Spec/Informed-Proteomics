using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Composition
{
    /// <summary>
    /// Averagine algorithm - creates isotopic envelopes based on an averaged elemental composition
    /// </summary>
    public class Averagine
    {
        // Ignore Spelling: Averagine, Isotopomer

        /// <summary>
        /// Map of NominalMass to Isotope Envelope
        /// </summary>
        private readonly ConcurrentDictionary<int, IsotopomerEnvelope> IsotopeEnvelopeMap;

        /// <summary>
        /// Constructor
        /// </summary>
        public Averagine()
        {
            IsotopeEnvelopeMap = new ConcurrentDictionary<int, IsotopomerEnvelope>();
        }

        /// <summary>
        /// Get the Isotopomer envelope for <paramref name="monoIsotopeMass"/>
        /// </summary>
        /// <param name="monoIsotopeMass"></param>
        /// <returns></returns>
        public static IsotopomerEnvelope GetIsotopomerEnvelope(double monoIsotopeMass)
        {
            return DefaultAveragine.GetIsotopomerEnvelopeInst(monoIsotopeMass);
        }

        /// <summary>
        /// Get the Isotopomer envelope for <paramref name="monoIsotopeMass"/> using <paramref name="isoProfilePredictor"/>
        /// </summary>
        /// <param name="monoIsotopeMass"></param>
        /// <param name="isoProfilePredictor"></param>
        /// <returns></returns>
        public IsotopomerEnvelope GetIsotopomerEnvelopeInst(double monoIsotopeMass, IsoProfilePredictor isoProfilePredictor = null)
        {
            var nominalMass = (int)Math.Round(monoIsotopeMass * Constants.RescalingConstant);
            return GetIsotopomerEnvelopeFromNominalMassInst(nominalMass, isoProfilePredictor);
        }

        /// <summary>
        /// Get the theoretical Isotope profile for <paramref name="monoIsotopeMass"/> at charge <paramref name="charge"/>
        /// </summary>
        /// <param name="monoIsotopeMass"></param>
        /// <param name="charge"></param>
        /// <param name="relativeIntensityThreshold"></param>
        /// <returns></returns>
        public static List<Peak> GetTheoreticalIsotopeProfile(double monoIsotopeMass, int charge, double relativeIntensityThreshold = 0.1)
        {
            return DefaultAveragine.GetTheoreticalIsotopeProfileInst(monoIsotopeMass, charge, relativeIntensityThreshold);
        }

        /// <summary>
        /// Get the theoretical Isotope profile for <paramref name="monoIsotopeMass"/> at charge <paramref name="charge"/> using <paramref name="isoProfilePredictor"/>
        /// </summary>
        /// <param name="monoIsotopeMass"></param>
        /// <param name="charge"></param>
        /// <param name="relativeIntensityThreshold"></param>
        /// <param name="isoProfilePredictor"></param>
        /// <returns></returns>
        public List<Peak> GetTheoreticalIsotopeProfileInst(double monoIsotopeMass, int charge, double relativeIntensityThreshold = 0.1, IsoProfilePredictor isoProfilePredictor = null)
        {
            var peakList = new List<Peak>();
            var envelope = GetIsotopomerEnvelopeInst(monoIsotopeMass, isoProfilePredictor);
            for (var isotopeIndex = 0; isotopeIndex < envelope.Envelope.Length; isotopeIndex++)
            {
                var intensity = envelope.Envelope[isotopeIndex];
                if (intensity < relativeIntensityThreshold)
                {
                    continue;
                }

                var mz = Ion.GetIsotopeMz(monoIsotopeMass, charge, isotopeIndex);
                peakList.Add(new Peak(mz, intensity));
            }
            return peakList;
        }

        /// <summary>
        /// Get the Isotopomer envelope for the nominal mass <paramref name="nominalMass"/>
        /// </summary>
        /// <param name="nominalMass"></param>
        /// <returns></returns>
        public static IsotopomerEnvelope GetIsotopomerEnvelopeFromNominalMass(int nominalMass)
        {
            return DefaultAveragine.GetIsotopomerEnvelopeFromNominalMassInst(nominalMass);
        }

        /// <summary>
        /// Get the Isotopomer envelope for the nominal mass <paramref name="nominalMass"/> using <paramref name="isoProfilePredictor"/>
        /// </summary>
        /// <param name="nominalMass"></param>
        /// <param name="isoProfilePredictor"></param>
        /// <returns></returns>
        public IsotopomerEnvelope GetIsotopomerEnvelopeFromNominalMassInst(int nominalMass, IsoProfilePredictor isoProfilePredictor = null)
        {
            var nominalMassFound = IsotopeEnvelopeMap.TryGetValue(nominalMass, out var envelope);
            if (nominalMassFound)
            {
                return envelope;
            }

            var mass = nominalMass / Constants.RescalingConstant;
            envelope = ComputeIsotopomerEnvelope(mass, isoProfilePredictor);
            IsotopeEnvelopeMap.AddOrUpdate(nominalMass, envelope, (key, value) => value);

            return envelope;
        }

        private const double C = 4.9384;
        private const double H = 7.7583;
        private const double N = 1.3577;
        private const double O = 1.4773;
        private const double S = 0.0417;
        private const double AveragineMass = C * Atom.C + H * Atom.H + N * Atom.N + O * Atom.O + S * Atom.S;

        /// <summary>
        /// Default averagine formula
        /// </summary>
        public static Averagine DefaultAveragine;

        static Averagine()
        {
            DefaultAveragine = new Averagine();
        }

        private IsotopomerEnvelope ComputeIsotopomerEnvelope(double mass, IsoProfilePredictor isoProfilePredictor = null)
        {
            var averagineCount = mass / AveragineMass;
            var numC = (int)Math.Round(C * averagineCount);
            var numH = (int)Math.Round(H * averagineCount);
            var numN = (int)Math.Round(N * averagineCount);
            var numO = (int)Math.Round(O * averagineCount);
            var numS = (int)Math.Round(S * averagineCount);

            if (numH == 0)
            {
                numH = 1;
            }

            isoProfilePredictor = isoProfilePredictor ?? IsoProfilePredictor.Predictor;
            return isoProfilePredictor.GetIsotopomerEnvelope(numC, numH, numN, numO, numS);
        }
    }
}
