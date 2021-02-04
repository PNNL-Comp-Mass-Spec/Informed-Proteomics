using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

// ReSharper disable UnusedMember.Global
namespace InformedProteomics.Backend.Data.Biology
{
    /// <summary>
    /// Ion class: a biological component with a charge
    /// </summary>
    public class Ion
    {
        /// <summary>
        /// Elemental composition of the ion
        /// </summary>
        public Composition.Composition Composition { get; }

        /// <summary>
        /// Electrical charge of the ion
        /// </summary>
        public int Charge { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="composition"></param>
        /// <param name="charge"></param>
        public Ion(Composition.Composition composition, int charge)
        {
            Composition = composition;
            Charge = charge;
        }

        /// <summary>
        /// Get the monoisotopic m/z of the ion
        /// </summary>
        /// <returns></returns>
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
        /// <param name="minimumRelativeIntensity">Minimum intensity threshold for including the isotope in the results</param>
        /// <returns>Enumerable of isotope peaks</returns>
        public IEnumerable<Isotope> GetIsotopes(double minimumRelativeIntensity)
        {
            var isotopeIndex = -1;
            foreach (var isotopeRatio in Composition.GetIsotopomerEnvelopeRelativeIntensities())
            {
                ++isotopeIndex;
                if (isotopeRatio >= minimumRelativeIntensity)
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

            Array.Sort(index, (i, j) => isotopes[j].CompareTo(isotopes[i]));

            for (var i = 0; i < numIsotopes && i < index.Length; i++)
            {
                yield return new Isotope(index[i], isotopes[index[i]]);
            }
        }

        /// <summary>
        /// Get the top 3 isotopes for this ion
        /// </summary>
        /// <returns></returns>
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
                for (var i = indexOfMostAbundantIsotope - 1; i <= indexOfMostAbundantIsotope + 1 && i < isotopes.Length; i++)
                {
                    top3.Add(new Isotope(i, isotopes[i]));
                }
            }

            return top3;
        }

        /// <summary>
        /// Get the m/z of the specified isotope
        /// </summary>
        /// <param name="monoIsotopicMass"></param>
        /// <param name="charge"></param>
        /// <param name="isotopeIndex"></param>
        /// <returns></returns>
        public static double GetIsotopeMz(double monoIsotopicMass, int charge, int isotopeIndex)
        {
            var isotopeMass = monoIsotopicMass + isotopeIndex * Constants.C13MinusC12;
            return isotopeMass / charge + Constants.Proton;
        }

        /// <summary>
        /// Get the monoisotopic mass of the specified isotope
        /// </summary>
        /// <param name="isotopeMz"></param>
        /// <param name="charge"></param>
        /// <param name="isotopeIndex"></param>
        /// <returns></returns>
        public static double GetMonoIsotopicMass(double isotopeMz, int charge, int isotopeIndex)
        {
            var isotopeMass = (isotopeMz - Constants.Proton) * charge;
            var monoIsotopeMass = isotopeMass - isotopeIndex * Constants.C13MinusC12;
            return monoIsotopeMass;
        }
    }
}
