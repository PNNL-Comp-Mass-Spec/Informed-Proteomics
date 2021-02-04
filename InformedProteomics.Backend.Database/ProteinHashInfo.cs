using System.Security.Cryptography;
using System.Text;

namespace InformedProteomics.Backend.Database
{
    internal class ProteinHashInfo
    {
        private static readonly SHA1 mSha1Hasher = new SHA1Managed();

        private static readonly StringBuilder mSha1StringBuilder = new StringBuilder();

        /// <summary>
        /// Number of times this protein name has been encountered in the source FASTA file
        /// </summary>
        public int ObservationCount { get; set; }

        /// <summary>
        /// Number of residues in the protein sequence
        /// </summary>
        public int SequenceLength { get; }

        /// <summary>
        /// SHA-1 Hash of the protein sequence
        /// </summary>
        public string SequenceHash { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        public ProteinHashInfo(string proteinSequence, int observationCount = 1)
        {
            ObservationCount = observationCount;
            SequenceLength = proteinSequence.Length;
            SequenceHash = Sha1Hash(proteinSequence);
        }

        /// <inheritdoc />
        public override string ToString()
        {
            return SequenceLength + " residues: " + SequenceHash;
        }

        /// <summary>
        /// Compute the SHA1 Hash of the given text
        /// </summary>
        /// <param name="text"></param>
        /// <returns>String representation of the SHA1 hash</returns>
        public static string Sha1Hash(string text)
        {
            var hash = mSha1Hasher.ComputeHash(Encoding.UTF8.GetBytes(text));

            mSha1StringBuilder.Clear();
            foreach (var b in hash)
            {
                mSha1StringBuilder.Append(b.ToString("X2"));
            }

            return mSha1StringBuilder.ToString();
        }
    }
}
