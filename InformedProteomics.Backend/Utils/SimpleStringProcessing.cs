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
    public static class SimpleStringProcessing
    {
        /// <summary>
        /// Random shuffle a string
        /// </summary>
        /// <param name="str"></param>
        /// <returns>Shuffled string</returns>
        public static string Shuffle(string str)
        {
            var indices = Enumerable.Range(0, str.Length).OrderBy(_ => _random.Next()).ToArray();
            var sflStr = new StringBuilder(str.Length);
            foreach (var index in indices)
            {
                sflStr.Append(str[index]);
            }

            return sflStr.ToString();
        }

        /// <summary>
        /// Perform a set number of random mutations on a string
        /// </summary>
        /// <param name="textToMutate"></param>
        /// <param name="numMutations"></param>
        /// <remarks>The given number of residues will be changed to a different residue</remarks>
        /// <returns>String with randomly changed amino acids</returns>
        public static string Mutate(string textToMutate, int numMutations)
        {
            var textLength = textToMutate.Length;

            // Populate a HashSet with list of indices of the residues that we will mutate
            // Use a HashSet to assure that there are no duplicates
            var selectedIndexSet = new HashSet<int>();
            var maxValues = Math.Min(numMutations, textLength);

            while (selectedIndexSet.Count < maxValues)
            {
                // This .Add command proceeds gracefully if the HashSet already contains the newly generated random number
                selectedIndexSet.Add(_random.Next(textLength));
            }

            var mutated = new StringBuilder(textLength);
            for (var i = 0; i < textLength; i++)
            {
                if (!selectedIndexSet.Contains(i))
                {
                    mutated.Append(textToMutate[i]);
                }
                else
                {
                    var mutatedResidue = textToMutate[i];
                    while (mutatedResidue == textToMutate[i])
                    {
                        mutatedResidue = AminoAcid.StandardAminoAcidCharacters[_random.Next(AminoAcid.StandardAminoAcidCharacters.Length)];
                    }
                    mutated.Append(mutatedResidue);
                }
            }
            return mutated.ToString();
        }

        /// <summary>
        /// Get the string between two periods, so A.PEPTIDE.J returns PEPTIDE
        /// </summary>
        /// <param name="str"></param>
        /// <returns>Sequence without the prefix or suffix residue</returns>
        public static string GetStringBetweenDots(string str)
        {
            //if (!Regex.IsMatch(str, @"^[A-Z" + FastaDatabase.Delimiter + @"]\.[A-Z]+\.[A-Z" + FastaDatabase.Delimiter + @"]$")) return null;
            var firstDotIndex = str.IndexOf('.');
            var lastDotIndex = str.LastIndexOf('.');
            if (firstDotIndex >= lastDotIndex)
            {
                return null;
            }

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
