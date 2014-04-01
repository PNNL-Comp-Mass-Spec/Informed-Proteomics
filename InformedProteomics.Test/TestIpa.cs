using System;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestIpa
    {
        [Test]
        public void TestEdrnTargets()
        {
            const string dbFilePath = @"D:\Research\Data\EDRN\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            var targetDb = new FastaDatabase(dbFilePath);
            var indexedDbTarget = new IndexedDatabase(targetDb);
            var annotations = indexedDbTarget.FullSequenceAnnotationsAndOffsets(6, 30, 2, 1, Enzyme.Trypsin);

            //const string spikedInPeptideFile = @"D:\Research\Data\EDRN\SpikedPeptides.txt";
            //var spikedInPeptides = File.ReadAllLines(spikedInPeptideFile);
            var spikedInPeptides = new[] {"FNTANDDNVTQVR"};

            var aminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);
            const int minCharge = 2;
            const int maxCharge = 4;
            var tolerance = new Tolerance(10, ToleranceUnit.Ppm);

            const string rawFilePath = @"D:\Research\Data\EDRN\RawFiles\DIA\342935_EDRN_Serum_07_DIA_2_12Nov13_Samwise_13-07-28.raw";
            var run = LcMsRun.GetLcMsRun(rawFilePath, MassSpecDataType.XCaliburRun);

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            run.ComputeMs1Features(minCharge, maxCharge, tolerance);
            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Finding MS1 feature: {0:f4} sec", sec);

            sw.Reset();
            sw.Start();

            Console.WriteLine("Peptide\tMonoMz\tCharge\tNumMs2Matches\tMs2Scans");
            var numPeptides = 0;
            var numValidScans = 0;
            foreach(var peptide in annotations.Select(a => a.Annotation.Substring(2, a.Annotation.Length-4)))
//            foreach (var peptide in spikedInPeptides)
            {
                if (++numPeptides%1000 == 0)
                {
                    Console.WriteLine("Processed {0} peptides.", numPeptides);
                    if (numPeptides == 100000) break;
                }
                var seqGraph = SequenceGraph.CreatePeptideGraph(aminoAcidSet, peptide);
                foreach (var sequenceComposition in seqGraph.GetSequenceCompositions())
                {
                    var peptideComposition = sequenceComposition + Composition.H2O;
                    peptideComposition.GetIsotopomerEnvelop();
                    for (var precursorCharge = minCharge; precursorCharge <= maxCharge; precursorCharge++)
                    {
                        var precursorIon = new Ion(peptideComposition, precursorCharge);
                        var validMs2ScanList = run.GetFragmentationSpectraScanNums(precursorIon)
                            .Where(scanNum => run.CheckMs1Signature(precursorIon, scanNum, tolerance)).ToList();
                        numValidScans += validMs2ScanList.Count;
                        //Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}",
                        //    peptide,
                        //    precursorIon.GetMonoIsotopicMz(),
                        //    precursorCharge,
                        //    validMs2ScanList.Count,
                        //    string.Join(",", validMs2ScanList)
                        //    );
                    }
                }
            }

            Console.WriteLine("NumValidScans: {0}", numValidScans);

            sw.Stop();
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }
    }
}
