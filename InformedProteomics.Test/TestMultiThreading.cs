using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
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

            const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_002216_235ACCEA.fasta";
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
            const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_002216_235ACCEA.fasta";
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
        [TestCase(1.5, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_002216_235ACCEA.fasta", 188961836)]  // 1.5MB
        [TestCase(3, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_005133_8491EFA2.fasta", 323719193)]  // 3MB
        [TestCase(6, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004530_B63BD900.fasta", 595227563)]  // 6MB
        [TestCase(15, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004208_295531A4.fasta", 1882434687)]  // 15MB
        public void TestSequenceEnumeration(double size, string dbFile, int expected)
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_002216_235ACCEA.fasta";  // 1.5MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_005133_8491EFA2.fasta";  // 3MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004530_B63BD900.fasta";  // 6MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004208_295531A4.fasta";  // 15MB
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);
            var numSequences = 0L;
            var timeDB = sw.Elapsed;
            Console.WriteLine("Read DB in " + timeDB.TotalSeconds + " Seconds");
            var estimatedAnnOff = indexedDb.EstimateTotalPeptides(0, 30, 250);
            var timeEstimate = sw.Elapsed;
            Console.WriteLine("Read Estimate in " + (timeEstimate - timeDB).TotalSeconds + " Seconds");
            //int coreCount = 0;
            //foreach (var item in new System.Management.ManagementObjectSearcher("Select NumberOfCores from Win32_Processor").Get())
            //{
            //    coreCount += int.Parse(item["NumberOfCores"].ToString());
            //}
            //Console.WriteLine("Number Of Cores: {0}", coreCount);
            //Console.WriteLine("Processors: " + System.Environment.ProcessorCount);
            Console.WriteLine("Estimated results: " + estimatedAnnOff);
            var annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzymeParallel(30, 250);
            var timeGetAnn = sw.Elapsed;
            Console.WriteLine("Read Annotations in " + (timeGetAnn - timeEstimate).TotalSeconds + " Seconds");
            /*/Parallel.ForEach(
                annotationsAndOffsets,
                //                new ParallelOptions { MaxDegreeOfParallelism = 2},
                annotationAndOffset =>
                {
                    Interlocked.Increment(ref numSequences);
                    //++numSequences;
                }
                );/**/
            //annotationsAndOffsets.Select(annotationsAndOffset => annotationsAndOffset.)
            // Below, original: 110, 109(total) seconds
            // Parallelizing AnnotationsAndOffsetsNoEnzyme: 86 seconds
            // Parallelizing AnnotationsAndOffsetsNoEnzyme, yield returns: 79.6, 94, 60, 60 seconds
            //
            // 3MB
            // serial: 
            // Parallel2: 107, 
            //
            // 6MB
            // serial: 
            // Parallel2: 
            //
            // 15MB
            // serial: 
            // Parallel2: 
            //using (var ofstream = new FileStream(Path.Combine(@"F:\InformedProteomicsTestFiles", Path.GetFileNameWithoutExtension(dbFile) + "_par.txt"), FileMode.Create))
            //using (var fout = new StreamWriter(ofstream))
            //{
            //    foreach (var annOff in annotationsAndOffsets)
            //    {
            //        numSequences++;
            //        fout.WriteLine(annOff.Annotation);
            //    }
            //}
            numSequences = annotationsAndOffsets.Count();
            var timeParForEach = sw.Elapsed;
            Console.WriteLine("Parallel ForEach in " + (timeParForEach - timeGetAnn).TotalSeconds + " Seconds");

            Console.WriteLine("NumPeptides: {0}", numSequences);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
            //Assert.AreEqual(188961836, numSequences);
            Assert.AreEqual(expected, numSequences);
        }

        [Test]
        [TestCase(1.5, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_002216_235ACCEA.fasta", 2399)]  // 1.5MB
        [TestCase(3, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_005133_8491EFA2.fasta", 3711)]  // 3MB
        [TestCase(6, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004530_B63BD900.fasta", 8898)]  // 6MB
        [TestCase(15, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004208_295531A4.fasta", 6334)]  // 15MB
        public void TestSequenceEnumerationIntact(double size, string dbFile, int expected)
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_002216_235ACCEA.fasta";  // 1.5MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_005133_8491EFA2.fasta";  // 3MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004530_B63BD900.fasta";  // 6MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004208_295531A4.fasta";  // 15MB
            const int numCTermCleavages = 0;
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);
            var numSequences = 0L;
            var timeDB = sw.Elapsed;
            Console.WriteLine("Read DB in " + timeDB.TotalSeconds + " Seconds");
            var estimatedAnnOff = indexedDb.EstimateTotalPeptides(2, 21, 300, 1, numCTermCleavages);
            var timeEstimate = sw.Elapsed;
            Console.WriteLine("Read Estimate in " + (timeEstimate - timeDB).TotalSeconds + " Seconds");
            Console.WriteLine("Estimated results: " + estimatedAnnOff);
            var annotationsAndOffsets = indexedDb.IntactSequenceAnnotationsAndOffsets(21, 300, numCTermCleavages);
            var timeGetAnn = sw.Elapsed;
            Console.WriteLine("Read Annotations in " + (timeGetAnn - timeEstimate).TotalSeconds + " Seconds");
            /*/Parallel.ForEach(
                annotationsAndOffsets,
                //                new ParallelOptions { MaxDegreeOfParallelism = 2},
                annotationAndOffset =>
                {
                    Interlocked.Increment(ref numSequences);
                    //++numSequences;
                }
                );/**/
            //using (var ofstream = new FileStream(Path.Combine(@"F:\InformedProteomicsTestFiles", Path.GetFileNameWithoutExtension(dbFile) + "_par.txt"), FileMode.Create))
            //using (var fout = new StreamWriter(ofstream))
            //{
            //    foreach (var annOff in annotationsAndOffsets)
            //    {
            //        numSequences++;
            //        fout.WriteLine(annOff.Annotation);
            //    }
            //}
            numSequences = annotationsAndOffsets.Count();
            var timeParForEach = sw.Elapsed;
            Console.WriteLine("Parallel ForEach in " + (timeParForEach - timeGetAnn).TotalSeconds + " Seconds");

            Console.WriteLine("NumPeptides: {0}", numSequences);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
            //Assert.AreEqual(188961836, numSequences);
            Assert.AreEqual(expected, numSequences);
        }

        [Test]
        [TestCase(1.5, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_002216_235ACCEA.fasta", 2700388)]  // 1.5MB
        [TestCase(3, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_005133_8491EFA2.fasta", 4165765)]  // 3MB
        [TestCase(6, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004530_B63BD900.fasta", 9146396)]  // 6MB
        [TestCase(15, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004208_295531A4.fasta", 14862126)]  // 15MB
        public void TestSequenceEnumerationNCTerm(double size, string dbFile, int expected)
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_002216_235ACCEA.fasta";  // 1.5MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_005133_8491EFA2.fasta";  // 3MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004530_B63BD900.fasta";  // 6MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004208_295531A4.fasta";  // 15MB
            const int numNTermCleavages = 1;
            const int numCTermCleavages = 0;
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);
            var numSequences = 0L;
            var timeDB = sw.Elapsed;
            Console.WriteLine("Read DB in " + timeDB.TotalSeconds + " Seconds");
            var estimatedAnnOff = indexedDb.EstimateTotalPeptides(1, 21, 300, numNTermCleavages, numCTermCleavages);
            var timeEstimate = sw.Elapsed;
            Console.WriteLine("Read Estimate in " + (timeEstimate - timeDB).TotalSeconds + " Seconds");
            Console.WriteLine("Estimated results: " + estimatedAnnOff);
            var annotationsAndOffsets = indexedDb.SequenceAnnotationsAndOffsetsWithNtermOrCtermCleavageNoLargerThan(21, 300, numNTermCleavages, numCTermCleavages);
            var timeGetAnn = sw.Elapsed;
            Console.WriteLine("Read Annotations in " + (timeGetAnn - timeEstimate).TotalSeconds + " Seconds");
            /*/Parallel.ForEach(
                annotationsAndOffsets,
                //                new ParallelOptions { MaxDegreeOfParallelism = 2},
                annotationAndOffset =>
                {
                    Interlocked.Increment(ref numSequences);
                    //++numSequences;
                }
                );/**/
            //using (var ofstream = new FileStream(Path.Combine(@"F:\InformedProteomicsTestFiles", Path.GetFileNameWithoutExtension(dbFile) + "_par.txt"), FileMode.Create))
            //using (var fout = new StreamWriter(ofstream))
            //{
            //    foreach (var annOff in annotationsAndOffsets)
            //    {
            //        numSequences++;
            //        fout.WriteLine(annOff.Annotation);
            //    }
            //}
            //foreach (var sao in annotationsAndOffsets)
            //{
            //    numSequences++;
            //}
            numSequences = annotationsAndOffsets.Count();
            var timeParForEach = sw.Elapsed;
            Console.WriteLine("Parallel ForEach in " + (timeParForEach - timeGetAnn).TotalSeconds + " Seconds");

            Console.WriteLine("NumPeptides: {0}", numSequences);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
            //Assert.AreEqual(188961836, numSequences);
            Assert.AreEqual(expected, numSequences);
        }

        [Test]
        [TestCase(1.5, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_002216_235ACCEA.fasta", 188961836)]  // 1.5MB
        [TestCase(3, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_005133_8491EFA2.fasta", 323719193)]  // 3MB
        [TestCase(6, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004530_B63BD900.fasta", 595227563)]  // 6MB
        [TestCase(15, @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004208_295531A4.fasta", 1882434687)]  // 15MB
        public void TestSequenceEnumerationSerial(double size, string dbFile, int expected)
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_002216_235ACCEA.fasta";  // 1.5MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_005133_8491EFA2.fasta";  // 3MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004530_B63BD900.fasta";  // 6MB
            //const string dbFile = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSPathFinderT\ID_004208_295531A4.fasta";  // 15MB
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);
            indexedDb.Read();
            var numSequences = 0L;
            var timeDB = sw.Elapsed;
            Console.WriteLine("Read DB in " + timeDB.TotalSeconds + " Seconds");
            var estimatedAnnOff = indexedDb.EstimateTotalPeptides(0, 30, 250);
            var timeEstimate = sw.Elapsed;
            Console.WriteLine("Read Estimate in " + (timeEstimate - timeDB).TotalSeconds + " Seconds");
            Console.WriteLine("Estimated results: " + estimatedAnnOff);
            var annotationsAndOffsets = indexedDb.AnnotationsAndOffsetsNoEnzyme(30, 250);
            var timeGetAnn = sw.Elapsed;
            Console.WriteLine("Read Annotations in " + (timeGetAnn - timeEstimate).TotalSeconds + " Seconds");
            //foreach (var annotationsAndOffset in annotationsAndOffsets)
            //{
            //    //Interlocked.Increment(ref numSequences);
            //    ++numSequences;
            //}
            using (
                var ofstream =
                    new FileStream(
                        Path.Combine(@"F:\InformedProteomicsTestFiles",
                            Path.GetFileNameWithoutExtension(dbFile) + "_old.txt"), FileMode.Create))
            using (var fout = new StreamWriter(ofstream))
            {
                foreach (var annOff in annotationsAndOffsets)
                {
                    numSequences++;
                    fout.WriteLine(annOff.Annotation);
                }
            }
            //numSequences = annotationsAndOffsets.Count();
            var timeParForEach = sw.Elapsed;
            Console.WriteLine("Parallel ForEach in " + (timeParForEach - timeGetAnn).TotalSeconds + " Seconds");

            Console.WriteLine("NumPeptides: {0}", numSequences);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
            //Assert.AreEqual(188961836, numSequences);
            Assert.AreEqual(expected, numSequences);
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
