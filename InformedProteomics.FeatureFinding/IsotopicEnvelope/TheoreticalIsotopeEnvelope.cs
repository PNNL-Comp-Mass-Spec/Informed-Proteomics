using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.FeatureFinding.IsotopicEnvelope
{
    public class TheoreticalIsotopeEnvelope : IsotopeEnvelope
    {
        public readonly int[] IndexOrderByRanking;
        public readonly int[] Ranking;

        public TheoreticalIsotopeEnvelope(double monoMass, int maxNumOfIsotopes, double relativeIntensityThreshold = 0.1)
        {
            MonoMass = monoMass;
            var isoEnv = Averagine.GetIsotopomerEnvelope(monoMass);
            var isotopeRankings = ArrayUtil.GetRankings(isoEnv.Envolope);

            Isotopes = new List<Isotope>(maxNumOfIsotopes);
            var ratioSum = 0d;
            for (var i = 0; i < isoEnv.Envolope.Length; i++)
            {
                if (isoEnv.Envolope[i] < relativeIntensityThreshold || isotopeRankings[i] > maxNumOfIsotopes) continue;
                ratioSum += isoEnv.Envolope[i];
                Isotopes.Add(new Isotope(i, isoEnv.Envolope[i]));
            }

            if (!(ratioSum > 0)) throw new Exception("Abnormal Theoretical Envelope");

            _probability            = new double[Isotopes.Count];
            Ranking                 = new int[Isotopes.Count];
            IndexOrderByRanking     = new int[Isotopes.Count];
            for(var i = 0; i < Isotopes.Count; i++)
            {
                _probability[i] = Isotopes[i].Ratio / ratioSum;
                Ranking[i] = isotopeRankings[Isotopes[i].Index];
                IndexOrderByRanking[isotopeRankings[Isotopes[i].Index] - 1] = i;
            }
        }

        public override double[] Probability
        {
            get { return _probability; }
        }

        public Isotope GetIsotopeRankedAt(int ranking)
        {
            return Isotopes[IndexOrderByRanking[ranking - 1]];
        }

        public override Isotope GetMostAbundantIsotope() { return GetIsotopeRankedAt(1); }

        // zero-based internal index
        public double GetIsotopeMz(int charge, int internalIndex)
        {
            return Ion.GetIsotopeMz(MonoMass, charge, Isotopes[internalIndex].Index);
        }

        public readonly List<Isotope> Isotopes;
        private readonly double[] _probability;
    }
}
