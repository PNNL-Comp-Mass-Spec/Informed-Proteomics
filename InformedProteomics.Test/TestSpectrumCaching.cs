using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestSpectrumCaching
    {
        [Test]
        public void TestAveragine()
        {
            for (var nominalMass = 1000; nominalMass <= 1000; nominalMass++)
            {
                Console.WriteLine("{0}\t{1}", nominalMass, 
                    string.Join(",", Averagine.GetIsotopomerEnvelopeFromNominalMass(nominalMass).Envolope.Select( v => string.Format("{0:f3}",v))));
            }
        }

        [Test]
        public void TestMs2Caching()
        {
            const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDown\raw\DataFiles\SBEP_STM_001_02272012_Aragon.raw";
            var run = LcMsRun.GetLcMsRun(rawFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

            const int minLength = 30;    // 7
            const int maxLength = 250; // 1000
            const int minPrecursorIonCharge = 3; // 3
            const int maxPrecursorIonCharge = 30;// 67
            const int minProductIonCharge = 1; 
            const int maxProductIonCharge = 10;

            var runCache = new CachedLcMsRun(run, 
                minPrecursorIonCharge, maxPrecursorIonCharge, minProductIonCharge, maxProductIonCharge,
                600.0, 1800.0, new Tolerance(10), new Tolerance(10));

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            runCache.CacheProductSpectra();
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void TestMs1Signature()
        {
            const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDown\raw\DataFiles";

            const string resultPath = @"C:\cygwin\home\kims336\Data\TopDown\raw\CorrMatches_N30";
            foreach (var resultFilePath in Directory.GetFiles(resultPath, "*.tsv"))
            {
                Console.WriteLine(resultFilePath);
            }

        }

        public void TestNominalMassErrors()
        {
            const int minLength = 300;
            const int maxLength = 400;

            var sw = new System.Diagnostics.Stopwatch();

//            const string dbFile = @"\\protoapps\UserData\Sangtae\TestData\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            const string dbFile = @"C:\cygwin\home\kims336\Data\TopDownJia\database\ID_003962_71E1A1D4.fasta";
            //const string dbFile = @"C:\cygwin\home\kims336\Data\TopDownJia\database\TargetProteins.fasta";
            var db = new FastaDatabase(dbFile);
            db.Read();
            var indexedDb = new IndexedDatabase(db);
            var numSequences = 0L;
            sw.Start();

            var hist = new long[11];
            var aaSet = new AminoAcidSet();
            foreach (var peptideAnnotationAndOffset in indexedDb.AnnotationsAndOffsetsNoEnzyme(minLength, maxLength))
            {
                ++numSequences;
                var annotation = peptideAnnotationAndOffset.Annotation;
                var sequenceStr = annotation.Substring(2, annotation.Length - 4);
                var sequenceComp = aaSet.GetComposition(sequenceStr);
                var mass = sequenceComp.Mass;
                var nominalMass = sequenceComp.NominalMass;
                var error = (int) Math.Round(mass*Constants.RescalingConstant) - nominalMass;
                var errorBin = error + hist.Length/2;
                if (errorBin < 0) errorBin = 0;
                if (errorBin >= hist.Length) errorBin = hist.Length - 1;
                hist[errorBin]++;
            }

            Console.WriteLine("NumSequences: {0}", numSequences);
            for (var i = 0; i < hist.Length; i++)
            {
                Console.WriteLine("{0}\t{1}\t{2}", i - hist.Length/2, hist[i], hist[i]/(double)numSequences);
            }

            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }
    }
}
