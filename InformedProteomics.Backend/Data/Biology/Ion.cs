using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Biology
{
    public class Ion
    {
        public Composition.Composition Composition { get; private set; }
        public int Charge { get; private set; }

        public Ion(Composition.Composition composition, int charge)
        {
            Composition = composition;
            Charge = charge;
        }

        public double GetMonoIsotopicMz()
        {
            return (Composition.Mass + Charge * Constants.Proton) / Charge;
        }

        /// <summary>
        /// Gets the m/z of ith isotope
        /// </summary>
        /// <param name="isotopeIndex">isotope index. 0 means mono-isotope, 1 means 2nd isotope, etc.</param>
        /// <returns></returns>
        public double GetIsotopeMz(int isotopeIndex)
        {
            return (Composition.GetIsotopeMass(isotopeIndex) + Charge * Constants.Proton) / Charge;
        }

        /// <summary>
        /// Gets the m/z of the most abundant isotope peak
        /// </summary>
        /// <returns>m/z of the most abundant isotope peak</returns>
        public double GetMostAbundantIsotopeMz()
        {
            return GetIsotopeMz(Composition.GetMostAbundantIsotopeZeroBasedIndex());
        }

        /// <summary>
        /// Gets the m/z of ith isotope
        /// </summary>
        /// <param name="isotopeIndexInRealNumber">isotope index in real number. 0 means mono-isotope, 0.5 means the center of mono and 2nd isotopes.</param>
        /// <returns></returns>
        public double GetIsotopeMz(double isotopeIndexInRealNumber)
        {
            return (Composition.GetIsotopeMass(isotopeIndexInRealNumber) + Charge * Constants.Proton) / Charge;
        }

        /// <summary>
        /// Gets theoretical isotope peaks whose intensities are relative to the most intense isotope
        /// </summary>
        /// <param name="relativeIntensityThreshold">relative isotope intensity threshold</param>
        /// <returns>Enumerable of isotope peaks</returns>
        public IEnumerable<Isotope> GetIsotopes(double relativeIntensityThreshold)
        {
            var isotopeIndex = -1;
            foreach (var isotopeRatio in Composition.GetIsotopomerEnvelopeRelativeIntensities())
            {
                ++isotopeIndex;
                if (isotopeRatio > relativeIntensityThreshold)
                {
                    yield return new Isotope(isotopeIndex, isotopeRatio);
                }
            }
        }

        /// <summary>
        /// Gets top n (numIsotopes) theoretical isotope peaks ordered by the ratios of isotopes (higher first)
        /// </summary>
        /// <param name="numIsotopes">number of isotopes</param>
        /// <returns>Enumerable of isotope peaks</returns>
        public IEnumerable<Isotope> GetIsotopes(int numIsotopes)
        {
            var isotopes = Composition.GetIsotopomerEnvelopeRelativeIntensities();
            var index = Enumerable.Range(0, isotopes.Length).ToArray();

            Array.Sort(index, (i,j) => isotopes[j].CompareTo(isotopes[i]));

            for (var i = 0; i < numIsotopes && i<index.Length; i++)
            {
                yield return new Isotope(index[i], isotopes[index[i]]);
            }
        }

        public IList<Isotope> GetTop3Isotopes()
        {
            var isotopes = Composition.GetIsotopomerEnvelopeRelativeIntensities();

            var top3 = new List<Isotope>();
            var indexOfMostAbundantIsotope = Composition.GetMostAbundantIsotopeZeroBasedIndex();
            if (indexOfMostAbundantIsotope == 0)
            {
                for (var i = 0; i < 3 && i < isotopes.Length; i++)
                {
                    top3.Add(new Isotope(i, isotopes[i]));
                }
            }
            else
            {
                for (var i = indexOfMostAbundantIsotope - 1; i <= indexOfMostAbundantIsotope+1 && i < isotopes.Length; i++)
                {
                    top3.Add(new Isotope(i, isotopes[i]));
                }
            }

            return top3;
        }

        public static double GetIsotopeMz(double monoIsotopicMass, int charge, int isotopeIndex)
        {
            var isotopeMass = monoIsotopicMass + isotopeIndex*Constants.C13MinusC12;
            return isotopeMass/charge + Constants.Proton;
        }

        public static double GetMonoIsotopicMass(double isotopeMz, int charge, int isotopeIndex)
        {
            var isotopeMass = (isotopeMz - Constants.Proton) * charge;
            var monoIsotopeMass = isotopeMass - isotopeIndex * Constants.C13MinusC12;
            return monoIsotopeMass;
        }

    }
}
