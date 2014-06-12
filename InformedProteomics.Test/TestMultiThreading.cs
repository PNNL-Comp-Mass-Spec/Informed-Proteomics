using System;
using System.Collections.Generic;
using System.Linq;
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
        public void TestSequenceEnumerationParallel2()
        {
            var sw = new System.Diagnostics.Stopwatch();

            const string dbFile = @"C:\cygwin\home\kims336\Data\TopDownQCShew\database\ID_002216_235ACCEA.fasta";
            var db = new FastaDatabase(dbFile);
            db.Read();
            var indexedDb = new IndexedDatabase(db);

            sw.Start();
            //var num = db.Characters().AsParallel().Count();
            var annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(7, 30);
            var num =
                annotationsAndOffsets.AsParallel().Count(annotationsAndOffset => annotationsAndOffset.Annotation.IndexOf('W') >= 0);

            //var n = annotationsAndOffsets.AsParallel().Select(a => a.Annotation.Length > 10);
            //var arr = n.ToArray();
            //var numSequences = arr.Length;
                
            //var numSequences = annotationsAndOffsets.AsParallel().LongCount();

            //from n in numbers
              //where Enumerable.Range(2, (int)Math.Sqrt(n)).All(i => n % i > 0)
              //select n;

            Console.WriteLine("NumPeptides: {0}", num);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }

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

        [Test]
        public void TestSequenceEnumerationSerial()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            const string dbFile = @"C:\cygwin\home\kims336\Data\TopDownQCShew\database\ID_002216_235ACCEA.fasta";
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);
            var numSequences = 0L;
            var annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(30, 250);
            foreach (var annotationsAndOffset in annotationsAndOffsets)
                {
                    //Interlocked.Increment(ref numSequences);
                    ++numSequences;
                }

            Console.WriteLine("NumPeptides: {0}", numSequences);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }

        [Test]
        public void TestPrimesParallel()
        {
            var numbers = Enumerable.Range(3, 10000000 - 3);

            //var parallelQuery =
            //  from n in numbers.AsParallel()
            //  where Enumerable.Range(2, (int)Math.Sqrt(n)).All(i => n % i > 0)
            //  select n;
            var parallelQuery =
                numbers.AsParallel().Count(n => Enumerable.Range(2, (int) Math.Sqrt(n)).All(i => n%i > 0));

//            var primes = parallelQuery.ToArray();
        }

        [Test]
        public void TestPrimesSerial()
        {
            var numbers = Enumerable.Range(3, 10000000 - 3);

            var parallelQuery =
                numbers.Count(n => Enumerable.Range(2, (int)Math.Sqrt(n)).All(i => n % i > 0));

            //var serialQuery =
            //  from n in numbers
            //  where Enumerable.Range(2, (int)Math.Sqrt(n)).All(i => n % i > 0)
            //  select n;

            //var primes = serialQuery.ToArray();
        }

    }
}
