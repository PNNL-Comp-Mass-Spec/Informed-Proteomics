using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Simple string processing functions
    /// </summary>
    public class SimpleStringProcessing
    {
        /// <summary>
        /// Random shuffle a string
        /// </summary>
        /// <param name="str"></param>
        /// <returns></returns>
        public static string Shuffle(string str)
        {
            var indices = Enumerable.Range(0, str.Length).OrderBy(r => _random.Next()).ToArray();
            var sflStr = new StringBuilder(str.Length);
            foreach (var index in indices) sflStr.Append(str[index]);
            return sflStr.ToString();
        }

        /// <summary>
        /// Perform a set number of random mutations on a string
        /// </summary>
        /// <param name="str"></param>
        /// <param name="numMutations"></param>
        /// <returns></returns>
        public static string Mutate(string str, int numMutations)
        {
            var length = str.Length;

            // Use a HashSet to assure that there are no duplicates
            var selectedIndexSet = new HashSet<int>();

            var maxIterations = numMutations * 5;
            var iteration = 0;
            while (selectedIndexSet.Count < numMutations && iteration < maxIterations)
            {
                // This .Add command proceeds gracefully if the HashSet already contains the newly generated random number
                selectedIndexSet.Add(_random.Next(length));
                iteration++;
            }

            var mutated = new StringBuilder(length);
            for (var i = 0; i < length; i++)
            {
                if (!selectedIndexSet.Contains(i))
                {
                    mutated.Append(str[i]);
                }
                else
                {
                    var mutatedResidue = str[i];
                    while (mutatedResidue == str[i])
                    {
                        mutatedResidue = AminoAcid.StandardAminoAcidCharacters[_random.Next(AminoAcid.StandardAminoAcidCharacters.Length)];
                    }
                    mutated.Append(mutatedResidue);
                }
            }
            return mutated.ToString();
        }

        /// <summary>
        /// Get the string between 2 periods, so A.PEPTIDE.J returns PEPTIDE
        /// </summary>
        /// <param name="str"></param>
        /// <returns></returns>
        public static string GetStringBetweenDots(string str)
        {
            //if (!Regex.IsMatch(str, @"^[A-Z" + FastaDatabase.Delimiter + @"]\.[A-Z]+\.[A-Z" + FastaDatabase.Delimiter + @"]$")) return null;
            var firstDotIndex = str.IndexOf('.');
            var lastDotIndex = str.LastIndexOf('.');
            if (firstDotIndex >= lastDotIndex) return null;
            return str.Substring(firstDotIndex + 1, lastDotIndex - firstDotIndex - 1);
        }

        /// <summary>
        /// Re-initialize the random number generator using the specified seed
        /// </summary>
        /// <param name="seed"></param>
        public static void DefineRandomNumberGeneratorSeed(int seed)
        {
            _random = new Random(seed);
        }

        private static Random _random = new Random();
    }
}
