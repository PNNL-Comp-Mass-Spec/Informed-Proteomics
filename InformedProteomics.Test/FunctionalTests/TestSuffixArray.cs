using System;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Database;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestSuffixArray
    {
        [Test]
        public void TestSearching()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            const string dbFile = @"\\protoapps\UserData\Sangtae\TestData\Databases\Short.fasta";
            var db = new FastaDatabase(dbFile);
            var searchableDb = new SearchableDatabase(db);
            const string pattern = "NSGSHFCGGSLINSQWVVSAAH";
            var position = searchableDb.Search(pattern);
            Assert.True(position >= 0);
            Console.WriteLine("Position: {0}", position);
            Console.WriteLine("Matched indices: {0}", string.Join(",", searchableDb.FindAllMatchedSequenceIndices(pattern)));
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }

        [Test]
        public void TestEnumeratingPeptides()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            const string dbFile = @"\\protoapps\UserData\Sangtae\TestData\Databases\Short.fasta";
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);
            foreach (var annotationAndOffset in indexedDb.AnnotationsAndOffsetsNoEnzyme(10, 20))
            {
                Console.WriteLine(annotationAndOffset.Annotation);
            }
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }

        [Test]
        public void TestEnumeratingProteins()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            const string dbFile = @"\\protoapps\UserData\Sangtae\TestData\Databases\Short.fasta";
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);
            foreach (var annotationAndOffset in indexedDb.IntactSequenceAnnotationsAndOffsetsWithCTermCleavagesLargerThan(10, 300, 3))
            {
                Console.WriteLine(annotationAndOffset.Annotation);
            }
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }

        [Test]
        public void TestCountingPeptides()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            const string dbFile = @"D:\Research\Data\CommonContaminants\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";

//            const string dbFile = @"C:\cygwin\home\kims336\Data\QCShew\ID_003456_9B916A8B.fasta";
//            const string dbFile = @"H:\Research\DDAPlus\database\Yeast_SGD_withContam.fasta";
//            const string dbFile = @"H:\Research\CPTAC_Phospho\database\ID_004208_295531A4.fasta";
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);
            //var numPeptides = indexedDb.IntactSequenceAnnotationsAndOffsets(21, 300, 0).LongCount()*31;
            var numPeptides = indexedDb
                    .SequenceAnnotationsAndOffsetsWithNtermOrCtermCleavageNoLargerThan(
                        21, 300, 1, 0).LongCount();
            //var numPeptides = indexedDb.AnnotationsAndOffsetsNoEnzyme(7, 150).LongCount();
            //var numPeptides =
            //    indexedDb.AnnotationsAndOffsets(7, 40, 2, 2, Enzyme.Trypsin).LongCount();


            //var numPeptides = indexedDb.AnnotationsAndOffsets(6, 40, 2, 2, Enzyme.Trypsin).LongCount();
            //var numPeptides = indexedDb.IntactSequenceAnnotationsAndOffsets(30, 250, 0).LongCount();
            //    .Select(annotationAndSequence => annotationAndSequence.Annotation.Length - 4)
            //    .Aggregate(0L, (current, length) => current + Math.Min(length - 29, 30));

            Console.WriteLine("NumPeptides: {0}", numPeptides);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }

        [Test]
        public void TestCountingProteoformsCloseToNTermOrCTerm()
        {
            const int minSequenceLength = 21;   // 21
            const int maxSequenceLength = 300;  // 300
            const int maxNumNTermCleavages = 1;
            const int maxNumCTermCleavages = 0;

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            //const string dbFile = @"C:\cygwin\home\kims336\Data\TopDownQCShew\database\ID_002216_235ACCEA.fasta";
            const string dbFile = @"\\protoapps\UserData\Sangtae\TestData\Databases\Short.fasta";
            var db = new FastaDatabase(dbFile);
            var indexedDb = new IndexedDatabase(db);

            var both = 0L;
            var nTermOnly = 0L;
            var cTermOnly = 0L;
            
            foreach (
                var annotationAndOffset in
                    indexedDb.IntactSequenceAnnotationsAndOffsets(minSequenceLength, int.MaxValue,
                        maxNumCTermCleavages))
            {
                // numCTermCleavages <= maxNumCTermCleavages
                var annotation = annotationAndOffset.Annotation;
                var length = (annotation.Length - 4);
                var numNTermCleavage = 0;
                int cleavedLength;
                while ((cleavedLength = length - numNTermCleavage) >= minSequenceLength)
                {
                    if (cleavedLength <= maxSequenceLength)
                    {
                        if (numNTermCleavage <= maxNumNTermCleavages)
                        {
                            ++both;
                        }
                        else
                        {
                            ++cTermOnly;
                        }
                        var anno = numNTermCleavage == 0 
                            ? annotation 
                            : string.Format("{0}.{1}", annotation[1 + numNTermCleavage], annotation.Substring(2 + numNTermCleavage));
                        Console.WriteLine(anno);
                    }
                    ++numNTermCleavage;
                }
            }

            foreach (
                var annotationAndOffset in
                    indexedDb.IntactSequenceAnnotationsAndOffsetsWithCTermCleavagesLargerThan(minSequenceLength, int.MaxValue,
                        maxNumCTermCleavages))
            {
                // numCTermCleavages > maxNumCTermCleavages
                var annotation = annotationAndOffset.Annotation;
                var length = (annotation.Length - 4);
                for (var numNTermCleavage = 0; numNTermCleavage <= maxNumNTermCleavages; numNTermCleavage++)
                {
                    var cleavedLength = length - numNTermCleavage;
                    if (cleavedLength >= minSequenceLength && cleavedLength <= maxSequenceLength)
                    {
                        ++nTermOnly;
                        var anno = numNTermCleavage == 0
                            ? annotation
                            : string.Format("{0}.{1}", annotation[1 + numNTermCleavage], annotation.Substring(2 + numNTermCleavage));
                        Console.WriteLine(anno);
                    }
                }
            }

            Console.WriteLine("Both: {0}", both);
            Console.WriteLine("N-term only: {0}", nTermOnly);
            Console.WriteLine("C-term only: {0}", cTermOnly);
            Console.WriteLine("All: {0}", both + nTermOnly + cTermOnly);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }

        [Test]
        public void TestGettingProteinLengthAndPosition()
        {
            const string dbFile = @"\\protoapps\UserData\Sangtae\TestData\Short.fasta";
            var db = new FastaDatabase(dbFile);
            db.Read();
            var indexedDb = new IndexedDatabase(db);
            foreach (var peptideAnnotationAndOffset in indexedDb.AnnotationsAndOffsets(6, 20, 2, 0, Enzyme.Trypsin))
            {
                var annotation = peptideAnnotationAndOffset.Annotation;
                var offset = peptideAnnotationAndOffset.Offset;
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}",
                    annotation, 
                    offset, 
                    db.GetProteinName(offset), 
                    db.GetProteinLength(db.GetProteinName(offset)), 
                    db.GetZeroBasedPositionInProtein(offset)+1);
            }
        }

        //[Test]
        //public void GetTrypticPeptideMassMzDistribution()
        //{
        //    const int minPeptideLength = 6;
        //    const int maxPeptideLength = 30;
        //    const int numTolerableTermini = 2;
        //    const int numMissedCleavages = 1;
        //    var enzyme = Enzyme.Trypsin;

        //    const string dbFilePath = @"\\protoapps\UserData\Sangtae\TestData\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
        //    var targetDb = new FastaDatabase(dbFilePath);

        //    var indexedDbTarget = new IndexedDatabase(targetDb);
        //    var aminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);

        //    var hist = new int[31];
        //    var numPeptides = 0;
        //    foreach (
        //        var annotationAndOffset in
        //            indexedDbTarget.FullSequenceAnnotationsAndOffsets(minPeptideLength, maxPeptideLength, numTolerableTermini,
        //                                               numMissedCleavages, enzyme))
        //    {
        //        var annotation = annotationAndOffset.Annotation;
        //        var pepSequence = annotation.Substring(2, annotation.Length - 4);
        //        var pepComposition = aminoAcidSet.GetComposition(pepSequence) + Composition.H2O;
        //        numPeptides++;
        //        for (var charge = 2; charge < 3; charge++)
        //        {
        //            var ion = new Ion(pepComposition, charge);
        //            var precursorMz = ion.GetMonoIsotopicMz();

        //            var massIndex = (int)Math.Round(precursorMz / 100.0);
        //            if (massIndex > 30) massIndex = 30;
        //            hist[massIndex]++;
        //        }
        //    }

        //    Console.WriteLine("AveragineMass\t#Pep\tRatio");
        //    for (var i = 1; i < hist.Length; i++)
        //    {
        //        Console.WriteLine("{0}\t{1}\t{2}", i*100, hist[i], hist[i]/(float)numPeptides);
        //    }
        //}

        //public void TestReadingBigFile()
        //{
        //    var sw = new System.Diagnostics.Stopwatch();
        //    sw.Start();
        //    const string bigDbFile = @"C:\cygwin\home\kims336\Research\SuffixArray\uniprot2012_7_ArchaeaBacteriaFungiSprotTrembl_2012-07-11.fasta";
        //    var lastLine = File.ReadLines(bigDbFile).Last();
        //    sw.Stop();

        //    var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
        //    Console.WriteLine(@"{0:f4} sec", sec);
        //    Console.WriteLine(lastLine);
        //}

        //[Test]
        //public void TestGeneratingDecoyDatabase()
        //{
        //    var sw = new System.Diagnostics.Stopwatch();
        //    sw.Start();

        //    const string dbFile = @"\\protoapps\UserData\Sangtae\TestData\BSA.fasta";
        //    var db = new FastaDatabase(dbFile);
        //    db.Decoy(Enzyme.Trypsin);
        //    sw.Stop();
        //    var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
        //    Console.WriteLine(@"{0:f4} sec", sec);
        //}
    }
}