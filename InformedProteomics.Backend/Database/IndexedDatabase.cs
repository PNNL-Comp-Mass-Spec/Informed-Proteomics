using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using SuffixArray;

namespace InformedProteomics.Backend.Database
{
    public class IndexedDatabase
    {
        public static readonly string PermutedLongestCommonPrefixFileExtension = ".icplcp";
        public static readonly Encoding Encoding = FastaDatabase.Encoding;

        private readonly string _pLcpFilePath;

        protected readonly FastaDatabase FastaDatabase;
        protected byte[] PLcp;

        public IndexedDatabase(FastaDatabase fastaDatabase)
        {
            FastaDatabase = fastaDatabase;
            var databaseFilePath = FastaDatabase.GetFastaFilePath();
            //var databaseFilePathNoExt = Path.GetDirectoryName(databaseFilePath) + Path.DirectorySeparatorChar + Path.GetFileNameWithoutExtension(databaseFilePath);
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

        public void Read()
        {
            using (var fileStream = new FileStream(_pLcpFilePath, FileMode.Open, FileAccess.Read))
            {
                PLcp = new byte[fileStream.Length - sizeof(int)];
                fileStream.Read(PLcp, 0, PLcp.Length);
            }
        }

        public IEnumerable<AnnotationAndOffset> AnnotationsAndOffsets(int minLength, int maxLength, int numTolerableTermini,
                                                      int numMissedCleavages, Enzyme enzyme)
        {
            return AnnotationsAndOffsets(minLength, maxLength, numTolerableTermini, numMissedCleavages, enzyme.Residues,
                               enzyme.IsNTerm);
        }

        public IEnumerable<AnnotationAndOffset> AnnotationsAndOffsetsParallel(int minLength, int maxLength, int numTolerableTermini,
                                                      int numMissedCleavages, Enzyme enzyme, int threads, CancellationToken? cancellationToken = null)
        {
            return AnnotationsAndOffsetsParallel(minLength, maxLength, numTolerableTermini, numMissedCleavages, enzyme.Residues,
                               enzyme.IsNTerm, threads, cancellationToken);
        }

        public IEnumerable<AnnotationAndOffset> AnnotationsAndOffsetsNoEnzyme(int minLength, int maxLength)
        {
            return AnnotationsAndOffsets(minLength, maxLength, 0, 0, null, false);
        }

        public IEnumerable<AnnotationAndOffset> AnnotationsAndOffsetsNoEnzymeParallel(int minLength, int maxLength, int threads = 0, CancellationToken? cancellationToken = null)
        {
            return AnnotationsAndOffsetsParallel(minLength, maxLength, 0, 0, null, false, threads, cancellationToken);
        }

        public IEnumerable<AnnotationAndOffset> IntactSequenceAnnotationsAndOffsets(int minLength, int maxLength)
        {
            return IntactSequenceAnnotationsAndOffsets(minLength, maxLength, 0);
        }

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

        public IEnumerable<AnnotationAndOffset> SequenceAnnotationsAndOffsetsWithNtermOrCtermCleavageNoLargerThan(
            int minSequenceLength, int maxSequenceLength
            , int maxNumNTermCleavages, int maxNumCTermCleavages)
        {
            foreach (
                var annotationAndOffset in
                    IntactSequenceAnnotationsAndOffsets(minSequenceLength, int.MaxValue,
                        maxNumCTermCleavages))
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
                    IntactSequenceAnnotationsAndOffsetsWithCTermCleavagesLargerThan(minSequenceLength, int.MaxValue,
                        maxNumCTermCleavages))
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
                if (residue == FastaDatabase.Delimiter)
                {
                    if (buf != null && buf.Count > 0) yield return new SequenceAndOffset(buf.ToArray(), curOffset);
                    //if (buf != null && buf.Count > 0) yield return new SequenceAndOffset(Encoding.GetString(buf.ToArray()), curOffset);
                    buf = new List<byte>();
                    curOffset = offset;
                }
                else
                {
                    if (buf != null) buf.Add(residue);
                }
            }
        }

        public long EstimateTotalPeptides(int mode, int minLength = 21, int maxLength = 300, int numNTermCleavages = 1, int numCTermCleavages = 0)
        {
            long count = 0;
            if (mode == 0)
            {
                var curSequence = new LinkedList<byte>();
                var lcpList = new LinkedList<byte>();
                var lcpEnum = PLcps().GetEnumerator();
                var fEnum = FastaDatabase.Characters().GetEnumerator();

                var seps = new Queue<IntWrapper>();
                bool read = false;
                while ((read = fEnum.MoveNext()) || curSequence.Count >= minLength)
                {
                    if (read)
                    {
                        lcpEnum.MoveNext();
                        curSequence.AddLast(fEnum.Current);
                        lcpList.AddLast(lcpEnum.Current);

                        if (fEnum.Current == FastaDatabase.Delimiter)
                        {
                            seps.Enqueue(new IntWrapper(curSequence.Count - 1));
                        }

                        if (curSequence.Count < maxLength + 2) continue;
                    }

                    if (seps.Count > 0 && seps.Peek().Value == 0)
                    {
                        seps.Dequeue();
                    }

                    var min = minLength > lcpList.First.Value ? minLength : lcpList.First.Value;
                    if (seps.Count == 0 || seps.Peek().Value >= maxLength + 2)
                    {
                        count += maxLength + 2 - min;
                    }
                    else if (seps.Peek().Value >= min)
                    {
                        count += seps.Peek().Value - min;
                    }

                    curSequence.RemoveFirst();
                    lcpList.RemoveFirst();
                    foreach (var sep in seps)
                    {
                        --sep.Value;
                    }
                }
            }
            else
            {
                // mode 2
                foreach (var sonc in SequencesWithOffsetNoCleavage())
                {
                    var seqLength = sonc.Sequence.Length;
                    // mode 2
                    for (int i = 0; i <= numCTermCleavages; i++)
                    {
                        // mode 2
                        if (mode == 2 && minLength <= seqLength - i && seqLength - i <= maxLength)
                        {
                            count++;
                        }
                        // mode 1 #1
                        if (mode == 1)
                        {
                            for (int j = 0; minLength <= seqLength - i - j; j++)
                            {
                                if (seqLength - i - j <= maxLength)
                                {
                                    count++;
                                }
                            }
                        }
                    }
                    if (mode == 1)
                    {
                        // mode 1 #2
                        for (int i = numCTermCleavages + 1; i <= seqLength - minLength; i++)
                        {
                            for (int j = 0; j <= numNTermCleavages; j++)
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
                    if (isCleavable.Length > FastaDatabase.Delimiter)
                        isCleavable[FastaDatabase.Delimiter] = true;
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
            foreach (var seqAndLcp in SequencesWithLcpAndOffset(minLength, maxLength+2))
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
                        if (ntt < numTolerableTermini-1) continue;

                        // Could be run in parallel, but probably not worth the cost.
                        for (var i = 1; i < seqArr.Length-1; i++)
                        {
                            var code = seqArr[i];
                            if (!isStandardAminoAcid[code]) break;
                            if (isCleavable[code]) ++nmc;

                            seqBuild.Append(code);
                            if (i >= minLength && i >= lcp)
                            {
                                if (ntt + (isCleavable[code] || seqArr[i+1] == FastaDatabase.Delimiter ? 1 : 0) >= numTolerableTermini)
                                {
                                    yield return new AnnotationAndOffset(offset, seqBuild + "." + seqArr[i + 1]);
                                }
                            }
                            if (nmc > numMissedCleavages) break;
                        }
                    }
                    else // N-term enzyme
                    {
                        if (seqArr[0] == FastaDatabase.Delimiter || isCleavable[seqArr[1]]) ++ntt;
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
                    for (var i = 1; i < seqArr.Length-1; i++)
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

        private IEnumerable<AnnotationAndOffset> AnnotationsAndOffsetsParallel(int minLength, int maxLength, int numTolerableTermini,
                                                      int numMissedCleavages, IEnumerable<char> enzymaticResidues,
                                                      bool isNTermEnzyme, int threads = 0, CancellationToken? cancellationToken = null
            )
        {
            var isCleavable = new bool[128];
            if (enzymaticResidues != null)
            {
                // Could be run in parallel, but probably not worth the cost.
                foreach (var residue in enzymaticResidues)
                {
                    isCleavable[residue] = true;
                    if (isCleavable.Length > FastaDatabase.Delimiter)
                        isCleavable[FastaDatabase.Delimiter] = true;
                }
            }

            var isStandardAminoAcid = new bool[256];
            // Could be run in parallel, but probably not worth the cost.
            foreach (var residue in AminoAcid.StandardAminoAcidCharacters)
            {
                isStandardAminoAcid[residue] = true;
            }

            // Try to get the number of physical cores in the system - requires System.Management.dll and a WMI query, but the performance penalty for 
            // using the number of logical processors in a hyperthreaded system is significant, and worse than the penalty for using fewer than all physical cores.
            int coreCount = 0;
            try
            {
                foreach (var item in new System.Management.ManagementObjectSearcher("Select NumberOfCores from Win32_Processor").Get())
                {
                    coreCount += int.Parse(item["NumberOfCores"].ToString());
                }
                //Console.WriteLine("Number Of Cores: {0}", coreCount);
            }
            catch (Exception)
            {
                // Use the logical processor count, divided by 2 to avoid the greater performance penalty of over-threading.
                coreCount = (int)(Math.Ceiling(System.Environment.ProcessorCount / 2.0));
            }

            if (threads <= 0 || threads > coreCount)
            {
                threads = coreCount;
            }
            //int prevThreads, prevPorts;
            //ThreadPool.GetMinThreads(out prevThreads, out prevPorts);
            //ThreadPool.SetMinThreads(8, prevPorts);
            CancellationToken token = cancellationToken != null ? cancellationToken.Value : CancellationToken.None;
            // pre, peptide sequence, next
            //return SequencesWithLcpAndOffset(minLength, maxLength + 2).AsParallel().WithDegreeOfParallelism(48).WithExecutionMode(ParallelExecutionMode.ForceParallelism).SelectMany(seqAndLcp => AnnotationsAndOffsetsParallelInternal(minLength, numTolerableTermini, numMissedCleavages, enzymaticResidues, isNTermEnzyme, seqAndLcp, isCleavable, isStandardAminoAcid));
            //return SequencesWithLcpAndOffset(minLength, maxLength + 2).AsParallel().WithExecutionMode(ParallelExecutionMode.ForceParallelism).SelectMany(seqAndLcp => AnnotationsAndOffsetsParallelInternal(minLength, numTolerableTermini, numMissedCleavages, enzymaticResidues, isNTermEnzyme, seqAndLcp, isCleavable, isStandardAminoAcid));
            //return SequencesWithLcpAndOffset(minLength, maxLength + 2).AsParallel().WithDegreeOfParallelism(4).WithExecutionMode(ParallelExecutionMode.ForceParallelism).SelectMany(seqAndLcp => AnnotationsAndOffsetsParallelInternal(minLength, numTolerableTermini, numMissedCleavages, enzymaticResidues, isNTermEnzyme, seqAndLcp, isCleavable, isStandardAminoAcid));
            return SequencesWithLcpAndOffset(minLength, maxLength + 2).AsParallel().WithDegreeOfParallelism(threads).WithCancellation(token).SelectMany(seqAndLcp => AnnotationsAndOffsetsParallelInternal(minLength, numTolerableTermini, numMissedCleavages, enzymaticResidues, isNTermEnzyme, seqAndLcp, isCleavable, isStandardAminoAcid));
        }

        private IEnumerable<AnnotationAndOffset> AnnotationsAndOffsetsParallelInternal(int minLength, int numTolerableTermini,
            int numMissedCleavages, IEnumerable<char> enzymaticResidues, bool isNTermEnzyme, SequenceLcpAndOffset seqAndLcp, bool[] isCleavable, bool[] isStandardAminoAcid
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
                    if (!(ntt < numTolerableTermini-1))
                    {
                        // Could be run in parallel, but probably not worth the cost.
                        for (var i = 1; i < seqArr.Length-1; i++)
                        {
                            var code = seqArr[i];
                            if (!isStandardAminoAcid[code]) break;
                            if (isCleavable[code]) ++nmc;

                            seqBuild.Append(code);
                            if (i >= minLength && i >= lcp)
                            {
                                if (ntt + (isCleavable[code] || seqArr[i+1] == FastaDatabase.Delimiter ? 1 : 0) >= numTolerableTermini)
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
                    if (seqArr[0] == FastaDatabase.Delimiter || isCleavable[seqArr[1]]) ++ntt;
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
                for (var i = 1; i < seqArr.Length-1; i++)
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
            var suffixArray = new int[sequence.Length-1];
            SAIS.sufsort(sequence, suffixArray, sequence.Length-1);
            
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

        public static byte GetLcp(IList<byte> sequence, int index1, int index2)
        {
            var lcp = (byte)0;

            while (sequence[index1+lcp] == sequence[index2+lcp])
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
        //                                              int numMissedCleavages, Enzyme enzyme)
        //{
        //    return NumSequences(minLength, maxLength, numTolerableTermini, numMissedCleavages, enzyme.Residues,
        //                       enzyme.IsNTerm);
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

    public class SequenceAndOffset
    {
        public SequenceAndOffset(byte[] sequence, long offset)
        {
            Sequence = sequence;
            Offset = offset;
        }

        public byte[] Sequence { set; get; }
        public long Offset { set; get; }
    }
}
