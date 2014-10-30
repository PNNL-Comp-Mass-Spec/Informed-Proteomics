using System;
using System.Collections.Generic;
using System.Diagnostics;
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
            var arr = db.Characters().ToArray();

            sw.Start();
            //var annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(7, 30);
            //            var num = annotationsAndOffsets.AsParallel().LongCount(annotationsAndOffset => annotationsAndOffset.Annotation.IndexOf('W') >= 0);
            //var num = annotationsAndOffsets.LongCount(annotationsAndOffset => annotationsAndOffset.Annotation.IndexOf('W') >= 0);
            //var num = arr.AsParallel().Where(c => c == 'W').LongCount();
            var num = 0;
            var sum = 0L;
            //foreach (var c in arr)
            for (var a = 0; a < arr.Length; a++)
            {
                var c = arr[a];
                for (var i = 0; i < c * 10000; i++) sum += i;
                //                Interlocked.Increment(ref num);
                if (++num == 1000) break;
            }

            Console.WriteLine("NumPeptides: {0}", sum);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }

        [Test]
        public void TestSumParallel()
        {
            //var array = Enumerable.Range(0, short.MaxValue).ToArray();
            const string dbFile = @"C:\cygwin\home\kims336\Data\TopDownQCShew\database\ID_002216_235ACCEA.fasta";
            var db = new FastaDatabase(dbFile);
            db.Read();
            //var indexedDb = new IndexedDatabase(db);
            //indexedDb.Read();
            //var peptides = indexedDb.AnnotationsAndOffsetsNoEnzyme(7, 30);
            var charArray = db.Characters().Select(c => (int)c);

            // Test methods.
            Console.WriteLine(SumAsParallel(charArray));
            Console.WriteLine(SumDefault(charArray));

            const int m = 100;
            var s1 = Stopwatch.StartNew();
            for (var i = 0; i < m; i++)
            {
                SumDefault(charArray);
            }
            s1.Stop();
            var s2 = Stopwatch.StartNew();
            for (var i = 0; i < m; i++)
            {
                SumAsParallel(charArray);
            }
            s2.Stop();
            Console.WriteLine((s1.Elapsed.TotalMilliseconds * 1000000 /
                m).ToString("0.00 ns"));
            Console.WriteLine((s2.Elapsed.TotalMilliseconds * 1000000 /
                m).ToString("0.00 ns"));
            Console.Read();
        }

        static int SumDefault(IEnumerable<int> array)
        {
            return array.Sum();
        }

        static int SumAsParallel(IEnumerable<int> array)
        {
            return array.AsParallel().AsUnordered().Sum();
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
                numbers.AsParallel().Count(n => Enumerable.Range(2, (int)Math.Sqrt(n)).All(i => n % i > 0));

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
