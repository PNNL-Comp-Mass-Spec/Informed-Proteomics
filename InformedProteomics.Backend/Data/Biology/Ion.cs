using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Biology
{
    public class Ion
    {
        public Composition Composition { get; private set; }
        public int Charge { get; private set; }

        public Ion(Composition composition, int charge)
        {
            Composition = composition;
            Charge = charge;
        }

        public double GetMz()
        {
            return (Composition.GetMass() + Charge * Constants.Proton) / Charge;
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
        public double GetBaseIsotopeMz()
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
        public IEnumerable<Tuple<int,float>> GetIsotopes(double relativeIntensityThreshold)
        {
            var isotopeIndex = -1;
            foreach (var isotopeRatio in Composition.GetIsotopomerEnvelop())
            {
                ++isotopeIndex;
                if (isotopeRatio > relativeIntensityThreshold)
                {
                    yield return new Tuple<int, float>(isotopeIndex, isotopeRatio);
                }
            }
        }

        /// <summary>
        /// Gets top n (numIsotopes) theoretical isotope peaks
        /// </summary>
        /// <param name="numIsotopes">number of isotopes</param>
        /// <returns>Enumerable of isotope peaks</returns>
        public IEnumerable<Tuple<int, float>> GetIsotopes(int numIsotopes)
        {
            var isotopes = Composition.GetIsotopomerEnvelop();
            var index = Enumerable.Range(0, isotopes.Length).ToArray();

            Array.Sort(index, (i,j) => isotopes[j].CompareTo(isotopes[i]));

            for (var i = 0; i < numIsotopes; i++)
            {
                yield return new Tuple<int, float>(index[i], isotopes[index[i]]);
            }
        }

    }
}
