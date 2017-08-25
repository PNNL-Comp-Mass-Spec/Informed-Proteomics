using System;
using System.Collections.Generic;
using SuffixArray;

namespace InformedProteomics.Backend.Database
{
    // In memory algorithm
    public class SearchableDatabase
    {
        public SearchableDatabase(FastaDatabase fastaDatabase)
        {
            FastaDatabase = fastaDatabase;
            _sequence = fastaDatabase.GetSequence();
            _suffixArray = new int[_sequence.Length];
            SAIS.sufsort(_sequence, _suffixArray, _sequence.Length);

            var neighboringLcps = new byte[_suffixArray.Length];
            neighboringLcps[0] = 0;

            for (var i = 1; i < _suffixArray.Length; i++)
            {
                var lcp = IndexedDatabase.GetLcp(_sequence, _suffixArray[i - 1], _suffixArray[i]);
                neighboringLcps[i] = lcp;
            }

            _leftLcps = new byte[_suffixArray.Length];
            _rightLcps = new byte[_suffixArray.Length];

            InitializeLcps(neighboringLcps, _leftLcps, _rightLcps, 0, _suffixArray.Length-1);
        }

        public FastaDatabase FastaDatabase { get; }

        /// <summary>
        /// Find all occurrences of pattern in the fasta sequence. Return all matched indices.
        /// </summary>
        /// <param name="patternStr">pattern to search</param>
        /// <returns>the matched sequence indices.</returns>
        public IList<int> FindAllMatchedSequenceIndices(String patternStr)
        {
            return FindAllMatchedSequenceIndices(FastaDatabase.Encoding.GetBytes(patternStr));
        }

        /// <summary>
        /// Find all occurrences of pattern in the fasta sequence. Return all matched indices.
        /// </summary>
        /// <param name="pattern">pattern to search</param>
        /// <returns>the matched sequence indices.</returns>
        public IList<int> FindAllMatchedSequenceIndices(byte[] pattern)
        {
            var matchIndex = Search(pattern);
            if (matchIndex < 0) return new List<int>();

            var matchedIndices = new List<int>();

            for (var i = matchIndex; i < _suffixArray.Length; i++)
            {
                var index = _suffixArray[i];
                var lcp = GetLcp(pattern, index, 0);

                if (lcp == pattern.Length) matchedIndices.Add(index);
                else break;
            }
            return matchedIndices;
        }

        /// <summary>
        /// Suffix array based O(m) search, where m is the pattern length
        /// </summary>
        /// <param name="patternStr">pattern to search</param>
        /// <returns>the relative position in this suffix array.</returns>
        public int Search(string patternStr)
        {
            return Search(FastaDatabase.Encoding.GetBytes(patternStr));
        }

        /// <summary>
        /// Suffix array based O(m) search, where m is the pattern length
        /// </summary>
        /// <param name="pattern">pattern to search</param>
        /// <returns>the relative position in this suffix array.</returns>
        public int Search(byte[] pattern)
        {
            // check that the pattern is within the left boundary
            var leftResult = Compare(pattern, _suffixArray[0], 0);
            if (Math.Abs(leftResult) - 1 == pattern.Length)
                return 0;              // exact leftmost match of the first element
            if (leftResult < 0)
                return -1;             // insertion point is at position 0

            // check that the pattern is within the right boundary
            var rightResult = Compare(pattern, _suffixArray[_suffixArray.Length - 1], 0);
            if (rightResult > 0) return -_suffixArray.Length;     // insertion point is at the end of the array

            // initialize the longest common prefixes values
            var queryLeftLcp = GetLcp(pattern, _suffixArray[0], 0);
            var queryRightLcp = GetLcp(pattern, _suffixArray[_suffixArray.Length - 1], 0);

            // indices for the binary search
            var leftIndex = 0;
            var rightIndex = _sequence.Length - 1;

            // loop invariant: element at leftIndex < pattern <= element at rightIndex
            while (rightIndex - leftIndex > 1)
            {
                var middleIndex = (leftIndex + rightIndex) / 2;
                if (queryLeftLcp >= queryRightLcp)
                {
                    var leftMiddleLcp = _leftLcps[middleIndex];
                    if (leftMiddleLcp > queryLeftLcp)
                    {       // and queryMiddle == queryLeft
                        leftIndex = middleIndex;
                        // queryLeft = queryMiddle, already true
                    }
                    else if (queryLeftLcp > leftMiddleLcp)
                    {  // and queryMiddle == leftMiddle
                        // we can conclude that query < middle because queryMiddle < queryLeft
                        queryRightLcp = leftMiddleLcp;
                        rightIndex = middleIndex;
                    }
                    else
                    {
                        // queryLeft == leftMiddle == queryMiddle
                        var middleResult = Compare(pattern, _suffixArray[middleIndex], queryLeftLcp);
                        if (middleResult <= 0)
                        {      // pattern <= middle
                            queryRightLcp = middleResult == 0 ? (byte)Math.Min(pattern.Length, byte.MaxValue) : (byte)(-middleResult - 1);
                            rightIndex = middleIndex;
                        }
                        else
                        {                       // middle < pattern
                            queryLeftLcp = (byte)(middleResult - 1);
                            leftIndex = middleIndex;
                        }
                    }
                }
                else
                {
                    // queryRight > queryLeft
                    //int middleRightLcp = this.middleRightLcps.get(middleIndex);
                    var middleRightLcp = _rightLcps[middleIndex];
                    if (middleRightLcp > queryRightLcp)
                    {           // and queryMiddle == queryRight
                        rightIndex = middleIndex;
                        // queryRight = queryMiddle, already true
                    }
                    else if (queryRightLcp > middleRightLcp)
                    {      // and queryMiddle == middleRight
                        queryLeftLcp = middleRightLcp;
                        leftIndex = middleIndex;
                    }
                    else
                    {
                        // middleRight == queryRight == queryMiddle
                        //int middleResult = Math.min(pattern.compareTo(factory.makeSuffix(indices.get(middleIndex)), queryRightLcp), Byte.MAX_VALUE);
                        var middleResult = Compare(pattern, _suffixArray[middleIndex], queryRightLcp);
                        if (middleResult <= 0)
                        {      // pattern <= middle
                            queryRightLcp = middleResult == 0 ? (byte)Math.Min(pattern.Length, byte.MaxValue) : (byte)(-middleResult - 1);
                            rightIndex = middleIndex;
                        }
                        else
                        {                       // middle < pattern
                            queryLeftLcp = (byte)(middleResult - 1);
                            leftIndex = middleIndex;
                        }
                    }
                }
            }

            // evaluate the base cases, found!
            if (queryRightLcp == pattern.Length) return rightIndex;

            // not found
            return -rightIndex - 1;
        }

        private readonly byte[] _sequence;
        private readonly int[] _suffixArray;
        //private readonly byte[] _neighboringLcps;   // neighboring Lcps
        private readonly byte[] _leftLcps;   // left Lcps
        private readonly byte[] _rightLcps;   // right Lcps

        private static byte InitializeLcps(IList<byte> nLcps, IList<byte> lLcps, IList<byte> rLcps, int start, int end)
        {
            if (end - start == 1)
            {
                // _neighboringLcps[index] encodes the LCP(index-1, index)
//                Console.WriteLine("Returning {0}", nLcps[end]);
                return nLcps[end];
            }

            // recursion
            var middleIndex = (start + end) / 2;
            var lLcp = InitializeLcps(nLcps, lLcps, rLcps, start, middleIndex);
            lLcps[middleIndex] = lLcp;
            var rLcp = InitializeLcps(nLcps, lLcps, rLcps, middleIndex, end);
            rLcps[middleIndex] = rLcp;

            // return the smallest one
//            Console.WriteLine("Returning {0}", lLcp < rLcp ? lLcp : rLcp);
            return lLcp < rLcp ? lLcp : rLcp;
        }

        /// <summary>
        /// Compares two suffices (index1 and index2)
        /// </summary>
        /// <param name="index">suffix index</param>
        /// <param name="pattern">sequence to compare</param>
        /// <param name="startIndex">known common prefix</param>
        /// <returns>a positive number if 1 is larger,
        /// a negative if 1 is smaller and 0 if they are equal.
        /// The longest common prefix length can be retrieved by taking absolute value of the return value minus 1
        /// </returns>
        private int Compare(IList<byte> pattern, int index, byte startIndex)
        {
            for (var offset = startIndex; offset <= byte.MaxValue; offset++)
            {
                if (offset >= pattern.Count) return -offset - 1;
                if (index + offset >= _sequence.Length) return offset + 1;

                var byte1 = pattern[offset];
                var byte2 = _sequence[index + offset];
                if (byte1 > byte2) return offset + 1;
                if (byte1 < byte2) return -offset - 1;
            }

            return byte.MaxValue;
        }

        public byte GetLcp(IList<byte> pattern, int index, byte startIndex)
        {
            for (var offset = startIndex; offset < pattern.Count; offset++)
            {
                if (pattern[offset] != _sequence[index + offset]) return offset;
                if (offset == byte.MaxValue) return offset;
            }

            return (byte)(pattern.Count);
        }
    }
}
