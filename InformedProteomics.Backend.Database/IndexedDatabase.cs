using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;
using SuffixArray;

namespace InformedProteomics.Backend.Database
{
    /// <summary>
    /// Provide fast, indexed access to Fasta database information using Permuted Longest Common Prefix data
    /// </summary>
    public class IndexedDatabase
    {
        // Ignore Spelling: Lcp, cterm, foreach, const, ntt, sep, nmc

        /// <summary>
        /// File extension to use for Permuted Longest Common Prefix file
        /// </summary>
        // ReSharper disable once StringLiteralTypo
        public static readonly string PermutedLongestCommonPrefixFileExtension = ".icplcp";

        /// <summary>
        /// Encoding to use for writing and reading indexed database files
        /// </summary>
        public static readonly Encoding Encoding = FastaDatabase.Encoding;

        private readonly string _pLcpFilePath;

        /// <summary>
        /// The Fasta database that will be indexed
        /// </summary>
        protected readonly FastaDatabase FastaDatabase;

        /// <summary>
        /// Cached Permuted Longest Common Prefix data
        /// </summary>
        protected byte[] PLcp;

        /// <summary>
        /// Constructor - build the index
        /// </summary>
        /// <param name="fastaDatabase"></param>
        public IndexedDatabase(FastaDatabase fastaDatabase)
        {
            FastaDatabase = fastaDatabase;
            var databaseFilePath = FastaDatabase.GetFastaFilePath();
            //var databaseFilePathNoExt = Path.Combine(Path.GetDirectoryName(databaseFilePath), Path.GetFileNameWithoutExtension(databaseFilePath));
            //_pLcpFilePath = databaseFilePathNoExt + PermutedLongestCommonPrefixFileExtension;
            _pLcpFilePath = Path.ChangeExtension(databaseFilePath, PermutedLongestCommonPrefixFileExtension);
            var lastWriteTimeHash = FastaDatabase.GetLastWriteTimeHash();

            if (!File.Exists(_pLcpFilePath) || !FastaDatabase.CheckHashCodeBinaryFile(_pLcpFilePath, lastWriteTimeHash))
            {
                Console.Write("Generating " + _pLcpFilePath + " ... ");
                CreatePermutedLongestCommonPrefixFile();
                Console.WriteLine("Done");
            }
        }

        /// <summary>
        /// Read in the Permuted Longest Common Prefix file
        /// </summary>
        public void Read()
        {
            using (var fileStream = new FileStream(_pLcpFilePath, FileMode.Open, FileAccess.Read))
            {
                PLcp = new byte[fileStream.Length - sizeof(int)];
                fileStream.Read(PLcp, 0, PLcp.Length);
            }
        }

        /// <summary>
        /// Get the annotation and offset data from the database for all sequences that comply with the parameters
        /// </summary>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <param name="numTolerableTermini"></param>
        /// <param name="numMissedCleavages"></param>
        /// <param name="enzyme"></param>
        /// <returns></returns>
        public IEnumerable<AnnotationAndOffset> AnnotationsAndOffsets(int minLength, int maxLength, int numTolerableTermini,
                                                      int numMissedCleavages, Enzyme enzyme)
        {
            return AnnotationsAndOffsets(minLength, maxLength, numTolerableTermini, numMissedCleavages, enzyme.Residues,
                               enzyme.IsNTerm);
        }

        /// <summary>
        /// Get the annotation and offset data from the database for all sequences that comply with the parameters
        /// </summary>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <param name="numTolerableTermini"></param>
        /// <param name="numMissedCleavages"></param>
        /// <param name="enzyme"></param>
        /// <param name="threads"></param>
        /// <param name="cancellationToken"></param>
        /// <returns></returns>
        public IEnumerable<AnnotationAndOffset> AnnotationsAndOffsetsParallel(int minLength, int maxLength, int numTolerableTermini,
                                                      int numMissedCleavages, Enzyme enzyme, int threads, CancellationToken? cancellationToken = null)
        {
            return AnnotationsAndOffsetsParallel(
                minLength, maxLength, numTolerableTermini, numMissedCleavages, enzyme.Residues, enzyme.IsNTerm, threads, cancellationToken);
        }

        /// <summary>
        /// Get the annotation and offset data from the database for all sequences that comply with the parameters
        /// </summary>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <returns></returns>
        public IEnumerable<AnnotationAndOffset> AnnotationsAndOffsetsNoEnzyme(int minLength, int maxLength)
        {
            return AnnotationsAndOffsets(minLength, maxLength, 0, 0, null, false);
        }

        /// <summary>
        /// Get the annotation and offset data from the database for all sequences that comply with the parameters
        /// </summary>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <param name="threads"></param>
        /// <param name="cancellationToken"></param>
        /// <returns></returns>
        public IEnumerable<AnnotationAndOffset> AnnotationsAndOffsetsNoEnzymeParallel(int minLength, int maxLength, int threads = 0, CancellationToken? cancellationToken = null)
        {
            return AnnotationsAndOffsetsParallel(minLength, maxLength, 0, 0, null, false, threads, cancellationToken);
        }

        /// <summary>
        /// Get the annotation and offset data from the database for all sequences that comply with the parameters
        /// </summary>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <returns></returns>
        public IEnumerable<AnnotationAndOffset> IntactSequenceAnnotationsAndOffsets(int minLength, int maxLength)
        {
            return IntactSequenceAnnotationsAndOffsets(minLength, maxLength, 0);
        }

        /// <summary>
        /// Get the annotation and offset data from the database for all sequences that comply with the parameters
        /// </summary>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <param name="numCTermCleavages"></param>
        /// <returns></returns>
        public IEnumerable<AnnotationAndOffset> IntactSequenceAnnotationsAndOffsets(int minLength, int maxLength, int numCTermCleavages)
        {
            foreach (var seqWithOffset in SequencesWithOffsetNoCleavage())
            {
                var seqArr = Encoding.GetString(seqWithOffset.Sequence);
                //var seqArr = seqWithOffset.Sequence;
                var offset = seqWithOffset.Offset;

                for (var i = 0; i <= numCTermCleavages; i++)
                {
                    var length = seqArr.Length - i;
                    if (length >= minLength && length <= maxLength)
                    {
                        yield return new AnnotationAndOffset(
                            offset,
                            "_." + seqArr.Substring(0, length) + "." + (i == 0 ? '_' : seqArr[length])
                            );
                    }
                }
            }
        }

        /// <summary>
        /// Get the annotation and offset data from the database for all sequences that comply with the parameters
        /// </summary>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <param name="numCTermCleavages"></param>
        /// <returns></returns>
        public IEnumerable<AnnotationAndOffset> IntactSequenceAnnotationsAndOffsetsWithCTermCleavagesLargerThan(int minLength, int maxLength, int numCTermCleavages)
        {
            foreach (var seqWithOffset in SequencesWithOffsetNoCleavage())
            {
                var seqArr = Encoding.GetString(seqWithOffset.Sequence);
                //var seqArr = seqWithOffset.Sequence;
                var offset = seqWithOffset.Offset;

                for (var i = numCTermCleavages + 1; i <= seqArr.Length - minLength; i++)
                {
                    var length = seqArr.Length - i;
                    if (length <= maxLength)
                    {
                        yield return new AnnotationAndOffset(
                            offset,
                            "_." + seqArr.Substring(0, length) + "." + (i == 0 ? '_' : seqArr[length])
                            );
                    }
                }
            }
        }

        /// <summary>
        /// Get the annotation and offset data from the database for all sequences that comply with the parameters
        /// </summary>
        /// <param name="minSequenceLength"></param>
        /// <param name="maxSequenceLength"></param>
        /// <param name="maxNumNTermCleavages"></param>
        /// <param name="maxNumCTermCleavages"></param>
        /// <returns></returns>
        public IEnumerable<AnnotationAndOffset> SequenceAnnotationsAndOffsetsWithNTermOrCTermCleavageNoLargerThan(
            int minSequenceLength, int maxSequenceLength, int maxNumNTermCleavages, int maxNumCTermCleavages)
        {
            foreach (
                var annotationAndOffset in
                IntactSequenceAnnotationsAndOffsets(minSequenceLength, int.MaxValue, maxNumCTermCleavages))
            {
                // numCTermCleavages <= maxNumCTermCleavages
                var annotation = annotationAndOffset.Annotation;
                var offset = annotationAndOffset.Offset;
                var length = (annotation.Length - 4);
                var numNTermCleavage = 0;
                int cleavedLength;
                while ((cleavedLength = length - numNTermCleavage) >= minSequenceLength)
                {
                    if (cleavedLength <= maxSequenceLength)
                    {
                        // ReSharper disable once CommentTypo
                        //  if (numNTermCleavage <= maxNumNTermCleavages) "both" else "cterm"
                        var cleavedAnnotation = numNTermCleavage == 0
                            ? annotation
                            : annotation[1 + numNTermCleavage] + "." + annotation.Substring(2 + numNTermCleavage);
                        yield return new AnnotationAndOffset(offset + numNTermCleavage, cleavedAnnotation);
                    }
                    ++numNTermCleavage;
                }
            }

            foreach (
                var annotationAndOffset in
                IntactSequenceAnnotationsAndOffsetsWithCTermCleavagesLargerThan(minSequenceLength, int.MaxValue, maxNumCTermCleavages))
            {
                // numCTermCleavages > maxNumCTermCleavages
                var annotation = annotationAndOffset.Annotation;
                var offset = annotationAndOffset.Offset;
                var length = (annotation.Length - 4);
                for (var numNTermCleavage = 0; numNTermCleavage <= maxNumNTermCleavages; numNTermCleavage++)
                {
                    var cleavedLength = length - numNTermCleavage;
                    if (cleavedLength >= minSequenceLength && cleavedLength <= maxSequenceLength)
                    {
                        // N only
                        var cleavedAnnotation = numNTermCleavage == 0
                            ? annotation
                            : annotation[1 + numNTermCleavage] + "." + annotation.Substring(2 + numNTermCleavage);
                        yield return new AnnotationAndOffset(offset + numNTermCleavage, cleavedAnnotation);
                    }
                }
            }
        }

        /// <summary>
        /// Length of the longest sequence
        /// </summary>
        /// <returns></returns>
        public int GetLongestSequenceLength()
        {
            return SequencesWithOffsetNoCleavage().Select(seqWithOffset => seqWithOffset.Sequence.Length).Max();
        }

        private IEnumerable<byte> PLcps()
        {
            if (PLcp == null)
            {
                const int bufferSize = 1 << 16;
                var buffer = new byte[bufferSize];
                var count = bufferSize;
                var numBytesRead = 0;

                using (var fileStream = new FileStream(_pLcpFilePath, FileMode.Open, FileAccess.Read))
                {
                    var numBytesToRead = fileStream.Length - sizeof(int);
                    while (count > 0)
                    {
                        count = fileStream.Read(buffer, 0, bufferSize);
                        for (var i = 0; i < count && numBytesRead++ < numBytesToRead; i++)
                        {
                            yield return buffer[i];
                        }
                    }
                }
            }
            else
            {
                foreach (var lcp in PLcp) yield return lcp;
            }
        }

        private IEnumerable<SequenceAndOffset> SequencesWithOffsetNoCleavage()
        {
            List<byte> buf = null;
            var curOffset = 0L;

            var offset = -1L;
            foreach (var residue in FastaDatabase.Characters())
            {
                ++offset;
                if (residue == FastaDatabaseConstants.Delimiter)
                {
                    if (buf != null && buf.Count > 0)
                        yield return new SequenceAndOffset(buf.ToArray(), curOffset);

                    buf = new List<byte>();
                    curOffset = offset;
                }
                else
                {
                    buf?.Add(residue);
                }
            }
        }

        /// <summary>
        /// Estimate the total number of peptides that will be used in processing - essential for reasonably accurate progress reporting
        /// </summary>
        /// <param name="mode"></param>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <param name="numNTermCleavages"></param>
        /// <param name="numCTermCleavages"></param>
        /// <returns></returns>
        public long EstimateTotalPeptides(int mode, int minLength = 21, int maxLength = 300, int numNTermCleavages = 1,
            int numCTermCleavages = 0)
        {
            InternalCleavageType cleavageType;
            if (mode == 0)
            {
                cleavageType = InternalCleavageType.MultipleInternalCleavages;
            }
            else if (mode == 2)
            {
                cleavageType = InternalCleavageType.NoInternalCleavage;
            }
            else
            {
                cleavageType = InternalCleavageType.SingleInternalCleavage;
            }

            return EstimateTotalPeptides(cleavageType, minLength, maxLength, numNTermCleavages, numCTermCleavages);
        }

        /// <summary>
        /// Estimate the total number of peptides that will be used in processing - essential for reasonably accurate progress reporting
        /// </summary>
        /// <param name="mode"></param>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <param name="numNTermCleavages"></param>
        /// <param name="numCTermCleavages"></param>
        /// <returns></returns>
        public long EstimateTotalPeptides(InternalCleavageType mode, int minLength = 21, int maxLength = 300, int numNTermCleavages = 1, int numCTermCleavages = 0)
        {
            long count = 0;
            if (mode == InternalCleavageType.MultipleInternalCleavages)
            {
                var curSequence = new LinkedList<byte>();
                var lcpList = new LinkedList<byte>();
                var lcpEnum = PLcps().GetEnumerator();
                var fEnum = FastaDatabase.Characters().GetEnumerator();

                // Use "IntWrapper" to allow modifying the value inside of the for each loop
                var sequences = new Queue<IntWrapper>();
                bool read;
                while ((read = fEnum.MoveNext()) || curSequence.Count >= minLength)
                {
                    if (read)
                    {
                        lcpEnum.MoveNext();
                        curSequence.AddLast(fEnum.Current);
                        lcpList.AddLast(lcpEnum.Current);

                        if (fEnum.Current == FastaDatabaseConstants.Delimiter)
                        {
                            sequences.Enqueue(new IntWrapper(curSequence.Count - 1));
                        }

                        if (curSequence.Count < maxLength + 2) continue;
                    }

                    if (sequences.Count > 0 && sequences.Peek().Value == 0)
                    {
                        sequences.Dequeue();
                    }

                    var min = minLength > lcpList.First.Value ? minLength : lcpList.First.Value;
                    if (sequences.Count == 0 || sequences.Peek().Value >= maxLength + 2)
                    {
                        count += maxLength + 2 - min;
                    }
                    else if (sequences.Peek().Value >= min)
                    {
                        count += sequences.Peek().Value - min;
                    }

                    curSequence.RemoveFirst();
                    lcpList.RemoveFirst();
                    foreach (var sequence in sequences)
                    {
                        --sequence.Value;
                    }
                }

                lcpEnum.Dispose();
                fEnum.Dispose();
            }
            else
            {
                // mode 2
                foreach (var sequenceItem in SequencesWithOffsetNoCleavage())
                {
                    var seqLength = sequenceItem.Sequence.Length;
                    // mode 2
                    for (var i = 0; i <= numCTermCleavages; i++)
                    {
                        // mode 2
                        if (mode == InternalCleavageType.NoInternalCleavage && minLength <= seqLength - i && seqLength - i <= maxLength)
                        {
                            count++;
                        }
                        // mode 1 #1
                        if (mode == InternalCleavageType.SingleInternalCleavage)
                        {
                            for (var j = 0; minLength <= seqLength - i - j; j++)
                            {
                                if (seqLength - i - j <= maxLength)
                                {
                                    count++;
                                }
                            }
                        }
                    }
                    if (mode == InternalCleavageType.SingleInternalCleavage)
                    {
                        // mode 1 #2
                        for (var i = numCTermCleavages + 1; i <= seqLength - minLength; i++)
                        {
                            for (var j = 0; j <= numNTermCleavages; j++)
                            {
                                if (minLength <= seqLength - i - j &&
                                    seqLength - i - j <= maxLength)
                                {
                                    count++;
                                }
                            }
                        }
                    }
                }
            }
            return count;
        }

        private IEnumerable<AnnotationAndOffset> AnnotationsAndOffsets(int minLength, int maxLength, int numTolerableTermini,
                                                      int numMissedCleavages, IEnumerable<char> enzymaticResidues,
                                                      bool isNTermEnzyme
            )
        {
            var isCleavable = new bool[128];
            if (enzymaticResidues != null)
            {
                // Could be run in parallel, but probably not worth the cost.
                foreach (var residue in enzymaticResidues)
                {
                    isCleavable[residue] = true;
                    if (isCleavable.Length > FastaDatabaseConstants.Delimiter)
                        isCleavable[FastaDatabaseConstants.Delimiter] = true;
                }
            }

            var isStandardAminoAcid = new bool[256];
            // Could be run in parallel, but probably not worth the cost.
            foreach (var residue in AminoAcid.StandardAminoAcidCharacters)
            {
                isStandardAminoAcid[residue] = true;
            }

            var seqBuild = new StringBuilder(maxLength + 3);

            // pre, peptide sequence, next
            // Should probably be run in parallel
            foreach (var seqAndLcp in SequencesWithLcpAndOffset(minLength, maxLength + 2))
            {
                var seqArr = Encoding.GetString(seqAndLcp.Sequence);
                //var seqArr = seqAndLcp.Sequence;
                var lcp = seqAndLcp.Lcp;
                var offset = seqAndLcp.Offset;
                seqBuild.Clear();
                seqBuild.Append(seqArr[0] + ".");

                if (enzymaticResidues != null)
                {
                    var ntt = 0;
                    var nmc = 0;
                    if (!isNTermEnzyme) // C-term enzyme
                    {
                        if (isCleavable[seqArr[0]]) ++ntt;
                        if (ntt < numTolerableTermini - 1) continue;

                        // Could be run in parallel, but probably not worth the cost.
                        for (var i = 1; i < seqArr.Length - 1; i++)
                        {
                            var code = seqArr[i];
                            if (!isStandardAminoAcid[code]) break;
                            if (isCleavable[code]) ++nmc;

                            seqBuild.Append(code);
                            if (i >= minLength && i >= lcp)
                            {
                                if (ntt + (isCleavable[code] || seqArr[i + 1] == FastaDatabaseConstants.Delimiter ? 1 : 0) >= numTolerableTermini)
                                {
                                    yield return new AnnotationAndOffset(offset, seqBuild + "." + seqArr[i + 1]);
                                }
                            }
                            if (nmc > numMissedCleavages) break;
                        }
                    }
                    else // N-term enzyme
                    {
                        if (seqArr[0] == FastaDatabaseConstants.Delimiter || isCleavable[seqArr[1]]) ++ntt;
                        if (ntt < numTolerableTermini - 1) continue;

                        // Could be run in parallel, but probably not worth the cost.
                        for (var i = 2; i < seqArr.Length - 1; i++)
                        {
                            var code = seqArr[i];
                            if (!isStandardAminoAcid[code]) break;
                            if (isCleavable[code]) ++nmc;
                            if (nmc > numMissedCleavages) break;

                            seqBuild.Append(code);
                            if (i >= minLength && i >= lcp)
                            {
                                if (ntt + (isCleavable[seqArr[i + 1]] ? 1 : 0) >= numTolerableTermini)
                                {
                                    yield return new AnnotationAndOffset(offset, seqBuild + "." + seqArr[i + 1]);
                                }
                            }
                        }
                    }
                }
                else    // No enzyme
                {
                    // Could be run in parallel, but probably not worth the cost.
                    for (var i = 1; i < seqArr.Length - 1; i++)
                    {
                        var code = seqArr[i];
                        if (!isStandardAminoAcid[code]) break;

                        seqBuild.Append(code);
                        if (i >= minLength && i >= lcp)
                        {
                            yield return new AnnotationAndOffset(offset, seqBuild + "." + seqArr[i + 1]);
                        }
                    }
                }
            }
        }

        internal class IntWrapper
        {
            public int Value;

            public IntWrapper(int value)
            {
                Value = value;
            }
        }

        private IEnumerable<AnnotationAndOffset> AnnotationsAndOffsetsParallel(
            int minLength, int maxLength, int numTolerableTermini,
            int numMissedCleavages, IEnumerable<char> enzymaticResidues,
            bool isNTermEnzyme, int threads = 0, CancellationToken? cancellationToken = null
            )
        {
            var isCleavable = new bool[128];

            var residues = enzymaticResidues?.ToList() ?? new List<char>();

            // Could be run in parallel, but probably not worth the cost.
            foreach (var residue in residues)
            {
                isCleavable[residue] = true;
                if (isCleavable.Length > FastaDatabaseConstants.Delimiter)
                    isCleavable[FastaDatabaseConstants.Delimiter] = true;
            }

            var isStandardAminoAcid = new bool[256];
            // Could be run in parallel, but probably not worth the cost.
            foreach (var residue in AminoAcid.StandardAminoAcidCharacters)
            {
                isStandardAminoAcid[residue] = true;
            }

            // Try to get the number of physical cores in the system - requires System.Management.dll and a WMI query, but the performance penalty for
            // using the number of logical processors in a hyperthreaded system is significant, and worse than the penalty for using fewer than all physical cores.
            if (threads <= 0 || threads > ParallelizationUtils.NumPhysicalCores)
            {
                threads = ParallelizationUtils.NumPhysicalCores;
            }

            var token = cancellationToken ?? CancellationToken.None;

            // pre, peptide sequence, next
            //return SequencesWithLcpAndOffset(minLength, maxLength + 2).AsParallel().WithDegreeOfParallelism(48).WithExecutionMode(ParallelExecutionMode.ForceParallelism).SelectMany(seqAndLcp => AnnotationsAndOffsetsParallelInternal(minLength, numTolerableTermini, numMissedCleavages, enzymaticResidues, isNTermEnzyme, seqAndLcp, isCleavable, isStandardAminoAcid));
            //return SequencesWithLcpAndOffset(minLength, maxLength + 2).AsParallel().WithExecutionMode(ParallelExecutionMode.ForceParallelism).SelectMany(seqAndLcp => AnnotationsAndOffsetsParallelInternal(minLength, numTolerableTermini, numMissedCleavages, enzymaticResidues, isNTermEnzyme, seqAndLcp, isCleavable, isStandardAminoAcid));
            //return SequencesWithLcpAndOffset(minLength, maxLength + 2).AsParallel().WithDegreeOfParallelism(4).WithExecutionMode(ParallelExecutionMode.ForceParallelism).SelectMany(seqAndLcp => AnnotationsAndOffsetsParallelInternal(minLength, numTolerableTermini, numMissedCleavages, enzymaticResidues, isNTermEnzyme, seqAndLcp, isCleavable, isStandardAminoAcid));
            return SequencesWithLcpAndOffset(minLength, maxLength + 2).AsParallel()
                .WithDegreeOfParallelism(threads).WithCancellation(token)
                .SelectMany(seqAndLcp =>
                    AnnotationsAndOffsetsParallelInternal(
                        minLength, numTolerableTermini, numMissedCleavages, residues, isNTermEnzyme, seqAndLcp, isCleavable, isStandardAminoAcid));
        }

        private IEnumerable<AnnotationAndOffset> AnnotationsAndOffsetsParallelInternal(int minLength, int numTolerableTermini,
            int numMissedCleavages,
            IEnumerable<char> enzymaticResidues,
            bool isNTermEnzyme,
            SequenceLcpAndOffset seqAndLcp,
            IReadOnlyList<bool> isCleavable,
            IReadOnlyList<bool> isStandardAminoAcid
            )
        {
            var seqArr = Encoding.GetString(seqAndLcp.Sequence);
            //var seqArr = seqAndLcp.Sequence;
            var lcp = seqAndLcp.Lcp;
            var offset = seqAndLcp.Offset;
            var seqBuild = new StringBuilder(seqArr.Length + 3);
            seqBuild.Append(seqArr[0] + ".");

            if (enzymaticResidues != null)
            {
                var ntt = 0;
                var nmc = 0;
                if (!isNTermEnzyme) // C-term enzyme
                {
                    if (isCleavable[seqArr[0]]) ++ntt;
                    if (!(ntt < numTolerableTermini - 1))
                    {
                        // Could be run in parallel, but probably not worth the cost.
                        for (var i = 1; i < seqArr.Length - 1; i++)
                        {
                            var code = seqArr[i];
                            if (!isStandardAminoAcid[code]) break;
                            if (isCleavable[code]) ++nmc;

                            seqBuild.Append(code);
                            if (i >= minLength && i >= lcp)
                            {
                                if (ntt + (isCleavable[code] || seqArr[i + 1] == FastaDatabaseConstants.Delimiter ? 1 : 0) >= numTolerableTermini)
                                {
                                    yield return new AnnotationAndOffset(offset, seqBuild + "." + seqArr[i + 1]);
                                }
                            }
                            if (nmc > numMissedCleavages) break;
                        }
                    }
                }
                else // N-term enzyme
                {
                    if (seqArr[0] == FastaDatabaseConstants.Delimiter || isCleavable[seqArr[1]]) ++ntt;
                    if (!(ntt < numTolerableTermini - 1))
                    {
                        // Could be run in parallel, but probably not worth the cost.
                        for (var i = 2; i < seqArr.Length - 1; i++)
                        {
                            var code = seqArr[i];
                            if (!isStandardAminoAcid[code]) break;
                            if (isCleavable[code]) ++nmc;
                            if (nmc > numMissedCleavages) break;

                            seqBuild.Append(code);
                            if (i >= minLength && i >= lcp)
                            {
                                if (ntt + (isCleavable[seqArr[i + 1]] ? 1 : 0) >= numTolerableTermini)
                                {
                                    yield return new AnnotationAndOffset(offset, seqBuild + "." + seqArr[i + 1]);
                                }
                            }
                        }
                    }
                }
            }
            else    // No enzyme
            {
                // Could be run in parallel, but probably not worth the cost.
                for (var i = 1; i < seqArr.Length - 1; i++)
                {
                    var code = seqArr[i];
                    if (!isStandardAminoAcid[code]) break;

                    seqBuild.Append(code);
                    if (i >= minLength && i >= lcp)
                    {
                        yield return new AnnotationAndOffset(offset, seqBuild + "." + seqArr[i + 1]);
                    }
                }
            }
        }

        private IEnumerable<SequenceLcpAndOffset> SequencesWithLcpAndOffset(int minLength, int maxLength)
        {
            var curSequence = new LinkedList<byte>();
            var lcpList = new LinkedList<byte>();
            var lcpEnum = PLcps().GetEnumerator();

            var offset = -1L;
            // Data dependency: cannot run in parallel; Also, not a significantly costly operation (< 10 seconds)
            foreach (var residue in FastaDatabase.Characters())
            {
                lcpEnum.MoveNext();
                curSequence.AddLast(residue);
                lcpList.AddLast(lcpEnum.Current);

                if (curSequence.Count < maxLength) continue;

                /**/
                var seqArr = new byte[curSequence.Count];
                curSequence.CopyTo(seqArr, 0);
                yield return new SequenceLcpAndOffset(seqArr, lcpList.First.Value, ++offset);
                //yield return new SequenceLcpAndOffset(curSequence.ToArray(), lcpList.First.Value, ++offset);
                /*/
                var seqArr = Encoding.GetString(curSequence.ToArray());
                var sep = seqArr.IndexOf('_', 1); // Is valid if separator is the prefix
                ++offset;
                if (sep == -1)
                {
                    yield return new SequenceLcpAndOffset(seqArr, lcpList.First.Value, offset);
                }
                else if (sep > minLength)
                {
                    yield return new SequenceLcpAndOffset(seqArr.Substring(0, sep + 1), lcpList.First.Value, offset);
                }
                /**/
                curSequence.RemoveFirst();
                lcpList.RemoveFirst();
            }

            lcpEnum.Dispose();

            // Data dependency: cannot run in parallel; Also, not a significantly costly operation (< 10 seconds)
            while (curSequence.Count >= minLength)
            {
                /**/
                var seqArr = new byte[curSequence.Count];
                curSequence.CopyTo(seqArr, 0);
                yield return new SequenceLcpAndOffset(seqArr, lcpList.First.Value, ++offset);
                //yield return new SequenceLcpAndOffset(curSequence.ToArray(), lcpList.First.Value, ++offset);
                /*/
                var seqArr = Encoding.GetString(curSequence.ToArray());
                var sep = seqArr.IndexOf('_', 1); // Is valid if separator is the prefix
                ++offset;
                if (sep == -1)
                {
                    yield return new SequenceLcpAndOffset(seqArr, lcpList.First.Value, offset);
                }
                else if (sep > minLength)
                {
                    yield return new SequenceLcpAndOffset(seqArr.Substring(0, sep + 1), lcpList.First.Value, offset);
                }
                /**/
                curSequence.RemoveFirst();
                lcpList.RemoveFirst();
            }
        }

        private void CreatePermutedLongestCommonPrefixFile()
        {
            if (File.Exists(_pLcpFilePath))
                File.Delete(_pLcpFilePath);

            var sequence = FastaDatabase.GetSequence();
            //Console.WriteLine("Annotation: {0}", System.Text.Encoding.ASCII.GetString(sequence));
            var suffixArray = new int[sequence.Length - 1];
            SAIS.sufsort(sequence, suffixArray, sequence.Length - 1);

            var prevIndex = sequence.Length - 1;

            var pLcp = new byte[suffixArray.Length];

            // Data dependency: cannot run in parallel
            foreach (var index in suffixArray)
            {
                var lcp = GetLcp(sequence, prevIndex, index);
                pLcp[index] = lcp;
                //Console.WriteLine("{0}\t{1}", System.Text.Encoding.ASCII.GetString(sequence, offset, sequence.Length - offset - 1), lcp);
                prevIndex = index;
            }

            using (var fs = new FileStream(_pLcpFilePath, FileMode.OpenOrCreate, FileAccess.Write))
            {
                foreach (var lcp in pLcp)
                {
                    //Console.WriteLine("LCP: {0}", lcp);
                    fs.WriteByte(lcp);
                }
                fs.Write(BitConverter.GetBytes(FastaDatabase.FileFormatId), 0, sizeof(int));
                fs.Write(BitConverter.GetBytes(FastaDatabase.GetLastWriteTimeHash()), 0, sizeof(int));
            }
        }

        /// <summary>
        /// Get the longest common prefix for the supplied sequence
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="index1"></param>
        /// <param name="index2"></param>
        /// <returns></returns>
        public static byte GetLcp(IList<byte> sequence, int index1, int index2)
        {
            var lcp = (byte)0;

            while (sequence[index1 + lcp] == sequence[index2 + lcp])
            {
                ++lcp;
                if (lcp == byte.MaxValue) break;
            }

            return lcp;
        }

        //public long NumSequences(int minLength, int maxLength, int numTolerableTermini,
        //                         int numMissedCleavages, char[] enzymaticResidues,
        //                         bool isNTermEnzyme)
        //{
        //    var numSequences = 0L;

        //    {
        //        var isCleavable = new bool[128];
        //        if (enzymaticResidues != null)
        //        {
        //            foreach (var residue in enzymaticResidues)
        //            {
        //                isCleavable[residue] = true;
        //                isCleavable[FastaDatabase.Delimiter] = true;
        //            }
        //        }

        // ReSharper disable once CommentTypo
        //        const string standardAAStr = "ACDEFGHIKLMNPQRSTVWY";
        //        var isStandardAminoAcid = new bool[128];
        //        foreach (var residue in standardAAStr)
        //        {
        //            isStandardAminoAcid[residue] = true;
        //        }

        //        // pre, peptide sequence, next
        //        foreach (var seqAndLcp in SequencesWithLcpAndOffset(minLength, maxLength + 2))
        //        {
        //            var seqArr = seqAndLcp.Annotation;
        //            var lcp = seqAndLcp.Lcp;
        //            var ntt = 0;
        //            var nmc = 0;

        //            if (enzymaticResidues != null)
        //            {
        //                if (!isNTermEnzyme) // C-term enzyme
        //                {
        //                    if (isCleavable[seqArr[0]]) ++ntt;
        //                    if (ntt < numTolerableTermini - 1) continue;

        //                    for (var i = 1; i < seqArr.Length - 1; i++)
        //                    {
        //                        var code = seqArr[i];
        //                        if (!isStandardAminoAcid[code]) break;
        //                        if (isCleavable[code]) ++nmc;

        //                        if (i >= minLength && i >= lcp)
        //                        {
        //                            if (ntt + (isCleavable[code] || i < seqArr.Length - 1 && seqArr[i + 1] == FastaDatabase.Delimiter ? 1 : 0) >= numTolerableTermini)
        //                            {
        //                                numSequences++;
        //                            }
        //                        }
        //                        if (nmc > numMissedCleavages) break;
        //                    }
        //                }
        //                else // N-term enzyme
        //                {
        //                }
        //            }
        //            else // no enzyme
        //            {
        //            }
        //        }
        //    }
        //    return numSequences;
        //}

        //public long NumSequences(int minLength, int maxLength, int numTolerableTermini,
        //                         int numMissedCleavages, Enzyme theEnzyme)
        //{
        //    return NumSequences(minLength, maxLength, numTolerableTermini, numMissedCleavages,
        //                       theEnzyme.Residues, theEnzyme.IsNTerm);
        //}
    }

    internal class SequenceLcpAndOffset
    {
        public SequenceLcpAndOffset(byte[] sequence, short lcp, long offset)
        {
            Sequence = sequence;
            Lcp = lcp;
            Offset = offset;
        }

        public byte[] Sequence { set; get; }
        public short Lcp { set; get; }
        public long Offset { set; get; }
    }

    /// <summary>
    /// Container for holding a sequence and its offset
    /// </summary>
    public class SequenceAndOffset
    {
        /// <summary>
        /// Constructor - set the data
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="offset"></param>
        public SequenceAndOffset(byte[] sequence, long offset)
        {
            Sequence = sequence;
            Offset = offset;
        }

        /// <summary>
        /// Sequence
        /// </summary>
        public byte[] Sequence { set; get; }

        /// <summary>
        /// Offset
        /// </summary>
        public long Offset { set; get; }
    }
}
