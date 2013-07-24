using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Database;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestIcTopDown
    {

        [Test]
        public void TestTopDownEnumeration()
        {
            // Search parameters
            const int numMissedCleavages = 25;
            const int minLength = 7;
            const int numMaxModsPerProtein = 6;
            //const string dbFilePath = @"C:\cygwin\home\kims336\Data\SuffixArray\BSA.fasta";
            const string dbFilePath =
                @"C:\cygwin\home\kims336\Data\SuffixArray\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            const int maxProteinLength = 10000;

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            // Configure amino acids
            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.ProteinNTerm, false);
            var dehydro = new SearchModification(Modification.PyroGluQ, 'C', SequenceLocation.Everywhere, false);
            var cysteinylC = new SearchModification(Modification.CysteinylC, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.GlutathioneC, 'C', SequenceLocation.Everywhere, false);
            //var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
                {
                    //pyroGluQ,
                    dehydro,
                    cysteinylC,
                    glutathioneC
                };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var targetDb = new FastaDatabase(dbFilePath);
            var indexedDb = new IndexedDatabase(targetDb);

            var numProteins = 0;
            long totalProtCompositions = 0;
            long totalProductCompositions = 0;
            foreach (var protAnnotation in indexedDb.SequencesAsStrings(numMissedCleavages, minLength))
            {
                ++numProteins;

                if (numProteins % 1000 == 0)
                {
                    Console.WriteLine("Processed {0} proteins", numProteins);
                }
                if (numProteins > 10000)
                    break;

                //if (numProteins < 28367) continue;
                //Console.WriteLine("{0} {1}", protAnnotation.Length, protAnnotation);

                var protSeq = protAnnotation.Substring(2, protAnnotation.Length - 4);
                var nTerm = protAnnotation[0] == FastaDatabase.Delimiter
                                      ? AminoAcid.ProteinNTerm
                                      : AminoAcid.PeptideNTerm;
                var cTerm = protAnnotation[protAnnotation.Length-1] == FastaDatabase.Delimiter
                                      ? AminoAcid.ProteinCTerm
                                      : AminoAcid.PeptideCTerm;

                var seqGraph = new SequenceGraph(aaSet, maxProteinLength);
                seqGraph.AddAminoAcid(nTerm.Residue);
                var isValidSequence = protSeq.All(seqGraph.AddAminoAcid);

                if (!isValidSequence) continue;

                seqGraph.AddAminoAcid(cTerm.Residue);

                //var numProtCompositions = seqGraph.GetSequenceCompositions().Length;
                var numProtCompositions = seqGraph.GetNumCompositions();
                totalProtCompositions += numProtCompositions;

                //long numProductCompositions = 0;
                //foreach (var scoringGraph in seqGraph.GetScoringGraphs())
                //{
                //    numProductCompositions += scoringGraph.GetCompositions().Count();
                //}
                //totalProductCompositions += numProductCompositions;
                //Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}",
                //    numProteins, numProtCompositions, numProductCompositions,
                //    totalProtCompositions, totalProductCompositions);
            }

            sw.Stop();
            Console.WriteLine("NumProteins: {0}", numProteins);
            Console.WriteLine("NumProteinCompositions: {0}", totalProtCompositions);
            Console.WriteLine("NumScoringGraphs: {0}", totalProductCompositions);
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void TestListProteinSequences()
        {
            // Search parameters
            const int numMissedCleavages = 25;
            const int minLength = 7;
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\SuffixArray\BSA.fasta";

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            var targetDb = new FastaDatabase(dbFilePath);
            var indexedDb = new IndexedDatabase(targetDb);

            var protIndex = 0;
            long numProtCompositions = 0;
            foreach (var protAnnotation in indexedDb.SequencesAsStrings(numMissedCleavages, minLength))
            {
                Console.WriteLine("{0} {1}", ++protIndex, protAnnotation);
            }

            sw.Stop();
            Console.WriteLine("NumProteins: {0}", protIndex);
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }
    }
}
