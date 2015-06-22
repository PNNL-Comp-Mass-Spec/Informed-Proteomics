using System;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Database;

namespace ParallelMSPathFinderT
{
    class Program
    {
        private const uint EnableExtendedFlags = 0x0080;
        [DllImport("kernel32.dll")]
        public static extern bool SetConsoleMode(IntPtr hConsoleHandle, uint dwMode);

        private static void Main(string[] args)
        {
            var handle = Process.GetCurrentProcess().MainWindowHandle;
            SetConsoleMode(handle, EnableExtendedFlags);
            TestCountingPeptides();
        }

        private static void TestCountingPeptides()
        {
            var aaSet = new AminoAcidSet();

            var sw = new Stopwatch();
            sw.Start();

            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_002166_F86E3B2F.fasta";
            const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_003456_9B916A8B.fasta";
            //            const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004208_295531A4.fasta";
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);
            indexedDb.Read();
            //var numPeptides = indexedDb.AnnotationsAndOffsetsNoEnzyme(7, 150).LongCount();
            var peptides =
                indexedDb.AnnotationsAndOffsets(7, 40, 2, 2, Enzyme.Trypsin);

            Parallel.ForEach(peptides, annotationAndOffset =>
            //foreach(var annotationAndOffset in peptides)
            {
                var annotation = annotationAndOffset.Annotation;
                var offset = annotationAndOffset.Offset;

                var graph = SequenceGraph.CreateGraph(aaSet, annotation);
            }
                )
            ;

//            Console.WriteLine("NumPeptides: {0}", numPeptides);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }
    }
}
