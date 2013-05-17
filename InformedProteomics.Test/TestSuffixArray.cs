using System;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Database;
using NUnit.Framework;
using SuffixArray;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestSuffixArray
    {

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
            //const string dbFile = @"C:\cygwin\home\kims336\Research\SuffixArray\uniprot_sprot_bacterial_ALLEntries_fungal_decoy_2009-05-28.fasta";
            //const string dbFile = @"D:\Research\Data\CompRef\H_sapiens_M_musculus_Trypsin_NCBI_Build37_2011-12-02.fasta";
            const string dbFile = @"..\..\..\TestFiles\BSA.fasta";
            var fileStream = new FileStream(dbFile, FileMode.Open, FileAccess.Read);
            var text = new byte[fileStream.Length];
            var suffixArray = new int[fileStream.Length];

            fileStream.Read(text, 0, text.Length);

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            SAIS.sufsort(text, suffixArray, text.Length);
            sw.Stop();

            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            System.Console.WriteLine(@"{0:f4} sec", sec);
        }
    }
}