using System;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Database;
using NUnit.Framework;
using SuffixArray;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestSuffixArray
    {
        [Test]
        public void GetTrypticPeptideMassMzDistribution()
        {
            const int minPeptideLength = 6;
            const int maxPeptideLength = 30;
            const int numTolerableTermini = 2;
            const int numMissedCleavages = 1;
            var enzyme = Enzyme.Trypsin;

            const string dbFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            var targetDb = new FastaDatabase(dbFilePath);

            var indexedDbTarget = new IndexedDatabase(targetDb);
            var aminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);

            var hist = new int[31];
            var numPeptides = 0;
            foreach (
                var annotationAndOffset in
                    indexedDbTarget.AnnotationsAndOffsets(minPeptideLength, maxPeptideLength, numTolerableTermini,
                                                       numMissedCleavages, enzyme))
            {
                var annotation = annotationAndOffset.Annotation;
                var pepSequence = annotation.Substring(2, annotation.Length - 4);
                var pepComposition = aminoAcidSet.GetComposition(pepSequence) + Composition.H2O;
                var mass = pepComposition.GetMass();
                numPeptides++;
                for (var charge = 2; charge < 3; charge++)
                {
                    var ion = new Ion(pepComposition, charge);
                    var precursorMz = ion.GetMonoIsotopicMz();

                    var massIndex = (int)Math.Round(precursorMz / 100.0);
                    if (massIndex > 30) massIndex = 30;
                    hist[massIndex]++;
                }
            }

            Console.WriteLine("Mass\t#Pep\tRatio");
            for (var i = 1; i < hist.Length; i++)
            {
                Console.WriteLine("{0}\t{1}\t{2}", i*100, hist[i], hist[i]/(float)numPeptides);
            }
        }

        [Test]
        public void TestReadingIndexedDatabase()
        {
            //const string dbFile = @"C:\cygwin\home\kims336\Research\SuffixArray\uniprot_sprot_bacterial_ALLEntries_fungal_decoy_2009-05-28.fasta";
            //const string dbFile = @"D:\Research\Data\CompRef\H_sapiens_M_musculus_Trypsin_NCBI_Build37_2011-12-02.fasta";
            const string dbFile = @"..\..\..\TestFiles\BSA.fasta";

            var database = new FastaDatabase(dbFile);
            database.Read();
            database.PrintSequence();
            Console.WriteLine("Done");
        }

        [Test]
        public void TestReadingBigFile()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            const string bigDbFile = @"C:\cygwin\home\kims336\Research\SuffixArray\uniprot2012_7_ArchaeaBacteriaFungiSprotTrembl_2012-07-11.fasta";
            var lastLine = File.ReadLines(bigDbFile).Last();
            sw.Stop();

            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            System.Console.WriteLine(@"{0:f4} sec", sec);
            System.Console.WriteLine(lastLine);
        }

        [Test]
        public void TestSuffixArrayGeneration()
        {
            const string dbFile = @"C:\cygwin\home\kims336\Data\SuffixArray\test.fasta";
            //const string dbFile = @"D:\Research\Data\CompRef\H_sapiens_M_musculus_Trypsin_NCBI_Build37_2011-12-02.fasta";
            //const string dbFile = @"..\..\..\TestFiles\BSA.fasta";

            var fileStream = new FileStream(dbFile, FileMode.Open, FileAccess.Read);
            var text = new byte[fileStream.Length];
            var suffixArray = new int[fileStream.Length];

            fileStream.Read(text, 0, text.Length);

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            SAIS.sufsort(text, suffixArray, text.Length);
            sw.Stop();

            Console.WriteLine("Text: " + Encoding.UTF8.GetString(text));
            foreach (var index in suffixArray)
            {
                Console.WriteLine(Encoding.ASCII.GetString(text, index, text.Length-index));
            }

            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            System.Console.WriteLine(@"{0:f4} sec", sec);
        }

        [Test]
        public void TestGeneratingDecoyDatabase()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            const string dbFile = @"C:\cygwin\home\kims336\Data\SuffixArray\BSA.fasta";
            var db = new FastaDatabase(dbFile);
            var decoyDb = db.Decoy(Enzyme.Trypsin);
            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            System.Console.WriteLine(@"{0:f4} sec", sec);
        }
    }
}