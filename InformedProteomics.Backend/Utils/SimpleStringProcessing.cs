using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Utils
{
    public class SimpleStringProcessing
    {
        public static string Shuffle(string str)
        {
            var indices = Enumerable.Range(0, str.Length).OrderBy(r => Random.Next()).ToArray();
            var sflStr = new StringBuilder(str.Length);
            foreach (var index in indices) sflStr.Append(str[index]);
            return sflStr.ToString();
        }

        public static string Mutate(string str, int numMutations)
        {
            var length = str.Length;
            var selectedIndexSet = new HashSet<int>();
            while(selectedIndexSet.Count < numMutations)
            {
                selectedIndexSet.Add(Random.Next(length));
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
                    char mutatedResidue = str[i];
                    while (mutatedResidue == str[i])
                    {
                        mutatedResidue = AminoAcid.StandardAminoAcidCharacters[Random.Next(AminoAcid.StandardAminoAcidCharacters.Length)];
                    }
                    mutated.Append(mutatedResidue);
                }
            }
            return mutated.ToString();
        }

        private static readonly Random Random = new Random();
    }
}
