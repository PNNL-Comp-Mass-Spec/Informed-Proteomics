using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestTopDownMs1Scoring
    {
        [Test]
        public void TestTopDownScoringForAllXics()
        {
            // Search parameters
            const int numNTermCleavages = 1;  // 30
            const int minLength = 7;
            const int maxLength = 1000;
            const int minCharge = 5; // 3
            const int maxCharge = 15; // 67
            const int numMaxModsPerProtein = 0; // 6
            var precursorTolerance = new Tolerance(10);
            const string dbFilePath = @"..\..\..\TestFiles\sprot.Ecoli.2012_07.fasta";
            //const string dbFilePath = @"..\..\..\TestFiles\sprot.Ecoli.2012_07.icdecoy.KR.fasta";
            
            //const string dbFilePath = @"..\..\..\TestFiles\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            // const string dbFilePath =
            //    @"C:\cygwin\home\kims336\Data\TopDown\ID_003558_56D73071.fasta";

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            Console.Write("Reading raw file...");
            const string specFilePath = @"C:\workspace\TopDown\E_coli_iscU_60_mock.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

            // Configure amino acid set
            //            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.ProteinNTerm, false);
            var dehydro = new SearchModification(Modification.PyroGluQ, 'C', SequenceLocation.Everywhere, false);
            var cysteinylC = new SearchModification(Modification.CysteinylC, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.GlutathioneC, 'C', SequenceLocation.Everywhere, false);
            //            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
                {
                    //pyroGluQ,
                    dehydro,
                    cysteinylC,
                    glutathioneC,
                    //oxM
                };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var targetDb = new FastaDatabase(dbFilePath);
         //   targetDb.CreateDecoyDatabase(Enzyme.Trypsin);
         //   System.Environment.Exit(1);
            var indexedDb = new IndexedDatabase(targetDb);

            var numProteins = 0;
            long totalProtCompositions = 0;
            long numXics = 0;
            TopDownScorer.MaxCharge = 25;
            TopDownScorer.MinCharge = 8;

            sw.Reset();
            sw.Start();
            Console.WriteLine("Generating XICs...");

            foreach (var protAnnotationAndOffset in indexedDb.AnnotationsAndOffsets(minLength, maxLength))
            {

                ++numProteins;
                //if (numProteins > 2000) break;

                if (numProteins % 1000 == 0)
                {
                    Console.WriteLine("Processed {0} proteins", numProteins);
                }

                //Console.WriteLine(protAnnotation);


                var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotationAndOffset.Annotation);
                
                //Console.WriteLine(seqGraph.GetSequenceCompositions()[0]);
                
                if (seqGraph == null) continue;

                for (var nTermCleavages = 0; nTermCleavages <= numNTermCleavages; nTermCleavages++)
                {
                    var protCompositions = seqGraph.GetSequenceCompositionsWithNTermCleavage(nTermCleavages);
                    foreach (var protComposition in protCompositions)
                    {
                        totalProtCompositions++;
                       // Console.WriteLine(protComposition);
                        var scorer = new TopDownScorer(protComposition, run, precursorTolerance, null);
                        var score = scorer.GetScore();

                        Console.WriteLine(score);


                    }
                }
            }

            sw.Stop();
            Console.WriteLine("NumProteins: {0}", numProteins);
            Console.WriteLine("NumProteinCompositions: {0}", totalProtCompositions);
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void TestMsAlignPlusResults()
        {
            TopDownScorer.MaxCharge = 25;
            TopDownScorer.MinCharge = 8;

            const string specFilePath = @"C:\workspace\TopDown\E_coli_iscU_60_mock.raw";
            const string msAlignPlusResultPath = @"C:\workspace\TopDown\E_coli_iscU_60_mock_MSAlign_ResultTable_sam.txt";

            var dehydro = new SearchModification(Modification.PyroGluQ, 'C', SequenceLocation.Everywhere, false);
            var cysteinylC = new SearchModification(Modification.CysteinylC, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.GlutathioneC, 'C', SequenceLocation.Everywhere, false);


            var searchModifications = new List<SearchModification>
                {
                    //pyroGluQ,
                    dehydro,
                    cysteinylC,
                    glutathioneC,
                    //oxM
                };
            var aaSet = new AminoAcidSet(searchModifications, 0);
            var precursorTolerance = new Tolerance(10);
            var run = new LcMsRun(new XCaliburReader(specFilePath));
            var writer = new StreamWriter(msAlignPlusResultPath+ ".txt");
            var reader = new StreamReader(msAlignPlusResultPath);

            string s;

            while ((s=reader.ReadLine())!=null)
            {
                if (s.StartsWith("Data_file_name	")) 
                { 
                    writer.WriteLine(s+"\tScore");
                    continue; 
                }
                var token = s.Split('\t');
                var annotation = token[13];
              //  Console.WriteLine("***\t" + annotation);
                var seqGraph = SequenceGraph.CreateGraph(aaSet, annotation);
                if (seqGraph == null)
                {
                    writer.WriteLine(s+"\tN/A");
                    continue;
                }

                var protCompositions = seqGraph.GetSequenceCompositions();
                
                var scorer = new TopDownScorer(protCompositions[0], run, precursorTolerance, null);
                var score = scorer.GetScore();

                writer.WriteLine(s+"\t"+score);
                Console.WriteLine(score);

            }
            

            writer.Close();
            reader.Close();
        }
        
        [Test]
        public void TestTopDownScoring()
        {
            TopDownScorer.MaxCharge = 25;
            TopDownScorer.MinCharge = 8;

            const string specFilePath = @"C:\workspace\TopDown\E_coli_iscU_60_mock.raw";
            const string protAnnotation = "A.AHAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK.";
            var dehydro = new SearchModification(Modification.PyroGluQ, 'C', SequenceLocation.Everywhere, false);
            var cysteinylC = new SearchModification(Modification.CysteinylC, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.GlutathioneC, 'C', SequenceLocation.Everywhere, false);

            
            var searchModifications = new List<SearchModification>
                {
                    //pyroGluQ,
                    dehydro,
                    cysteinylC,
                    glutathioneC,
                    //oxM
                };
            //var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var aaSet = new AminoAcidSet(searchModifications, 0);

            
            var precursorTolerance = new Tolerance(10);
            //Console.WriteLine(aaSet.GetAminoAcid('C').GetComposition());
            // Create a sequence graph
            //var protSeq = protAnnotation.Substring(2, protAnnotation.Length - 4);

            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);

          //  TopDownScorer.MaxCharge = 60;
          //  TopDownScorer.MinCharge = 3;
            var run = new LcMsRun(new XCaliburReader(specFilePath));
            
            foreach (var protComposition in seqGraph.GetSequenceCompositions())
            {
                var mostAbundantIsotopeIndex = protComposition.GetMostAbundantIsotopeZeroBasedIndex();
                Console.WriteLine("Composition\t{0}", protComposition);
                Console.WriteLine("MostAbundantIsotopeIndex\t{0}", mostAbundantIsotopeIndex);

                Console.WriteLine(new Ion(protComposition + Composition.H2O, 11).GetIsotopeMz(mostAbundantIsotopeIndex));

                Console.WriteLine();

                //for (var charge = TopDownScorer.MinCharge; charge <= TopDownScorer.MaxCharge; charge++)
                //{
                var scorer = new TopDownScorer(protComposition, run, precursorTolerance, null);
                var score = scorer.GetScore();

                Console.WriteLine(score);
                //var precursorIon = new Ion(protComposition + Composition.H2O, charge);
                //var xic = run.GetExtractedIonChromatogram(precursorIon.GetIsotopeMz(mostAbundantIsotopeIndex), precursorTolerance);
                //Console.WriteLine(xic[0].ScanNum + " " + xic[1].ScanNum);

                //Console.WriteLine("ScanNum\t{0}", string.Join("\t", xic.Select(p => p.ScanNum.ToString())));
                //Console.WriteLine("Charge " + charge + "\t" + string.Join("\t", xic.Select(p => p.Intensity.ToString())));
                // }

                Console.WriteLine("\nCharge\tm/z");
                for (var charge = 9; charge <= 18; charge++)
                {
                    var precursorIon = new Ion(protComposition + Composition.H2O, charge);
                    Console.WriteLine("{0}\t{1}", charge, precursorIon.GetIsotopeMz(mostAbundantIsotopeIndex));
                }
            }

            // sw.Stop();
            // var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            // Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }
    }
}
