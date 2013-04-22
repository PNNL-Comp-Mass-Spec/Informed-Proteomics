using System.IO;
using NUnit.Framework;
using SuffixArray;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestSuffixArray
    {
        [Test]
        public void TestSuffixArrayGeneration()
        {
            const string dbFile = @"..\..\..\TestFiles\Shewanella.fasta";
            //const string dbFile = @"..\..\..\TestFiles\BSA.fasta";
            //const string dbFile = @"C:\cygwin\home\kims336\Research\SuffixArray\uniprot_sprot_bacterial_ALLEntries_fungal_decoy_2009-05-28.fasta";
            //const string dbFile = @"D:\Research\Data\CompRef\H_sapiens_M_musculus_Trypsin_NCBI_Build37_2011-12-02.fasta";
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