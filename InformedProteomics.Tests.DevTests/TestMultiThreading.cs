using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Database;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.DevTests
{
    [TestFixture]
    internal class TestMultiThreading
    {
        [Test]
        public void TestSequenceEnumerationParallel2()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var sw = new Stopwatch();

            var fastaFile = Utils.GetTestFile(methodName, Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"MSPathFinderT\ID_002216_235ACCEA.fasta"));

            var db = new FastaDatabase(fastaFile.FullName);
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

            Console.WriteLine(@"{0:f4} sec", sw.Elapsed.TotalSeconds);
        }

        [Test]
        public void TestSumParallel()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            //var array = Enumerable.Range(0, short.MaxValue).ToArray();
            var fastaFile = Utils.GetTestFile(methodName, Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"MSPathFinderT\ID_002216_235ACCEA.fasta"));

            var db = new FastaDatabase(fastaFile.FullName);
            db.Read();
            //var indexedDb = new IndexedDatabase(db);
            //indexedDb.Read();
            //var peptides = indexedDb.AnnotationsAndOffsetsNoEnzyme(7, 30);
            var charArray = db.Characters().Select(c => (int)c).ToList();

            // Test methods.
            var defaultSum = SumAsParallel(charArray);
            var parallelSum = SumAsParallel(charArray);

            Console.WriteLine("Default sum {0}", defaultSum);
            Console.WriteLine("Parallel sum {0}", parallelSum);

            Assert.AreEqual(parallelSum, defaultSum);

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

            Console.WriteLine("{0:F2} msec/sum, on average for default", s1.Elapsed.TotalMilliseconds / m);
            Console.WriteLine("{0:F2} msec/sum, on average for parallel", s2.Elapsed.TotalMilliseconds / m);
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
        [TestCase(@"TEST_FOLDER\MSPathFinderT\ID_002216_235ACCEA.fasta", 188961836)] // 1.5MB
        public void TestSequenceEnumeration(string dbFile, int expectedPeptideCount)
        {
            TestSequenceEnumerationWork(dbFile, expectedPeptideCount);
        }


        [Test]
        [TestCase(@"TEST_FOLDER\MSPathFinderT\ID_005133_8491EFA2.fasta", 323719193)] // 3MB
        //[TestCase(@"TEST_FOLDER\MSPathFinderT\ID_004530_B63BD900.fasta", 595227563)]  // 6MB
        //[TestCase(@"TEST_FOLDER\MSPathFinderT\ID_004208_295531A4.fasta", 1882434687)]  // 15MB
        [Category("Long_Running")]
        public void TestSequenceEnumerationAddnl(string dbFile, int expectedPeptideCount)
        {
            TestSequenceEnumerationWork(dbFile, expectedPeptideCount);
        }

        private void TestSequenceEnumerationWork(string dbFile, int expectedPeptideCount)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName, dbFile);

            var fastaFile = Utils.GetTestFile(methodName, dbFile.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));

            var sw = new Stopwatch();
            sw.Start();

            var db = new FastaDatabase(fastaFile.FullName);
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
            //using (var ofstream = new FileStream(Path.Combine(@"F:\InformedProteomicsTestFiles", Path.GetFileNameWithoutExtension(fastaFile) + "_par.txt"), FileMode.Create))
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

            Console.WriteLine(@"{0:f4} sec", sw.Elapsed.TotalSeconds);

            Assert.AreEqual(expectedPeptideCount, numSequences);
        }

        [Test]
        [TestCase(1.5, @"TEST_FOLDER\MSPathFinderT\ID_002216_235ACCEA.fasta", 2399)]  // 1.5MB
        [TestCase(3, @"TEST_FOLDER\MSPathFinderT\ID_005133_8491EFA2.fasta", 3711)]  // 3MB
        [TestCase(6, @"TEST_FOLDER\MSPathFinderT\ID_004530_B63BD900.fasta", 8898)]  // 6MB
        //[TestCase(15, @"TEST_FOLDER\MSPathFinderT\ID_004208_295531A4.fasta", 6334)]  // 15MB
        public void TestSequenceEnumerationIntact(double size, string dbFile, int expected)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName, dbFile);

            var fastaFile = Utils.GetTestFile(methodName, dbFile.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));

            var sw = new Stopwatch();
            sw.Start();

            const int numCTermCleavages = 0;
            var db = new FastaDatabase(fastaFile.FullName);
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
            //using (var ofstream = new FileStream(Path.Combine(@"F:\InformedProteomicsTestFiles", Path.GetFileNameWithoutExtension(fastaFile) + "_par.txt"), FileMode.Create))
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

            Console.WriteLine(@"{0:f4} sec", sw.Elapsed.TotalSeconds);
            //Assert.AreEqual(188961836, numSequences);
            Assert.AreEqual(expected, numSequences);
        }

        [Test]
        [TestCase(1.5, @"TEST_FOLDER\MSPathFinderT\ID_002216_235ACCEA.fasta", 2700388)]  // 1.5MB
        [TestCase(3, @"TEST_FOLDER\MSPathFinderT\ID_005133_8491EFA2.fasta", 4165765)]  // 3MB
        [TestCase(6, @"TEST_FOLDER\MSPathFinderT\ID_004530_B63BD900.fasta", 9146396)]  // 6MB
        [TestCase(15, @"TEST_FOLDER\MSPathFinderT\ID_004208_295531A4.fasta", 14862126)]  // 15MB
        public void TestSequenceEnumerationNCTerm(double size, string dbFile, int expected)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName, dbFile);

            var fastaFile = Utils.GetTestFile(methodName, dbFile.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));

            var sw = new Stopwatch();
            sw.Start();

            const int numNTermCleavages = 1;
            const int numCTermCleavages = 0;
            var db = new FastaDatabase(fastaFile.FullName);
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
            //using (var ofstream = new FileStream(Path.Combine(@"F:\InformedProteomicsTestFiles", Path.GetFileNameWithoutExtension(fastaFile) + "_par.txt"), FileMode.Create))
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

            Console.WriteLine(@"{0:f4} sec", sw.Elapsed.TotalSeconds);
            //Assert.AreEqual(188961836, numSequences);
            Assert.AreEqual(expected, numSequences);
        }

        [Test]
        [TestCase(1.5, @"TEST_FOLDER\MSPathFinderT\ID_002216_235ACCEA.fasta", 188961836)]  // 1.5MB
        //[TestCase(3, @"TEST_FOLDER\MSPathFinderT\ID_005133_8491EFA2.fasta", 323719193)]  // 3MB
        //[TestCase(6, @"TEST_FOLDER\MSPathFinderT\ID_004530_B63BD900.fasta", 595227563)]  // 6MB
        //[TestCase(15, @"TEST_FOLDER\MSPathFinderT\ID_004208_295531A4.fasta", 1882434687)]  // 15MB
        public void TestSequenceEnumerationSerial(double size, string dbFile, int expected)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName, dbFile);

            var fastaFile = Utils.GetTestFile(methodName, dbFile.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));

            var sw = new Stopwatch();
            sw.Start();

            var db = new FastaDatabase(fastaFile.FullName);
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
            //using (
            //    var ofstream =
            //        new FileStream(
            //            Path.Combine(@"F:\InformedProteomicsTestFiles",
            //                Path.GetFileNameWithoutExtension(fastaFile) + "_old.txt"), FileMode.Create))
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

            Console.WriteLine(@"{0:f4} sec", sw.Elapsed.TotalSeconds);
            //Assert.AreEqual(188961836, numSequences);
            Assert.AreEqual(expected, numSequences);
        }

        [Test]
        [TestCase(10000, 1228, 20)]
        [TestCase(1000000, 78497, 20)]
        [TestCase(10000000, 664578, 0)]
        public void TestPrimesParallel(int ceiling, int primeCountExpected, int resultsToPreview)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var sw = new Stopwatch();
            sw.Start();

            var maxValue = ceiling - 3;
            var numbers = Enumerable.Range(3, maxValue);

            if (resultsToPreview <= 0)
            {

                var primeCountParallel =
                    numbers.AsParallel().Count(n => Enumerable.Range(2, (int)Math.Sqrt(n)).All(i => n % i > 0));

                sw.Stop();
                Console.WriteLine("Took {0:F2} seconds to find {1:N0} primes between 3 and {2}", sw.Elapsed.TotalSeconds, primeCountParallel,
                                  maxValue);

                Assert.AreEqual(primeCountExpected, primeCountParallel);

            }
            else
            {
                var index = 0;
                foreach (var item in
                    from n in numbers.AsParallel()
                    where Enumerable.Range(2, (int)Math.Sqrt(n)).All(i => n % i > 0)
                    select n)
                {
                    Console.WriteLine("{0}: {1}", index, item);
                    index++;

                    if (index >= resultsToPreview)
                        break;
                }
            }

        }

        [Test]
        [TestCase(10000, 1228, 20)]
        [TestCase(1000000, 78497, 20)]
        [TestCase(10000000, 664578, 0)]
        public void TestPrimesSerial(int ceiling, int primeCountExpected, int resultsToPreview)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var sw = new Stopwatch();
            sw.Start();

            var maxValue = ceiling - 3;
            var numbers = Enumerable.Range(3, maxValue);


            if (resultsToPreview <= 0)
            {
                var primeCountSerial =
                    numbers.Count(n => Enumerable.Range(2, (int)Math.Sqrt(n)).All(i => n % i > 0));


                sw.Stop();
                Console.WriteLine("Took {0:F2} seconds to find {1:N0} primes between 3 and {2}", sw.Elapsed.TotalSeconds, primeCountSerial, maxValue);

                Assert.AreEqual(primeCountExpected, primeCountSerial);

            }
            else
            {
                var index = 0;
                foreach (var item in
                    from n in numbers
                    where Enumerable.Range(2, (int)Math.Sqrt(n)).All(i => n % i > 0)
                    select n)
                {
                    Console.WriteLine("{0}: {1}", index, item);
                    index++;

                    if (index >= resultsToPreview)
                        break;
                }

            }

        }

    }
}
