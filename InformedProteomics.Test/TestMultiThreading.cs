using System;
using System.Threading;
using System.Threading.Tasks;
using InformedProteomics.Backend.Database;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestMultiThreading
    {
        [Test]
        public void TestSequenceEnumeration()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            const string dbFile = @"C:\cygwin\home\kims336\Data\TopDownQCShew\database\ID_002216_235ACCEA.fasta";
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);
            var numSequences = 0L;
            var annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(30, 250);
            Parallel.ForEach(
                annotationsAndOffsets, 
//                new ParallelOptions { MaxDegreeOfParallelism = 2},
                annotationAndOffset =>
                {
                    Interlocked.Increment(ref numSequences);
                    //++numSequences;
                }
                );

            Console.WriteLine("NumPeptides: {0}", numSequences);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
            
        }
    }
}
