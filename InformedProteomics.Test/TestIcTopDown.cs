using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
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
    internal class TestIcTopDown
    {
        [Test]
        public void TestTopDownSearchOneProtein()
        {
            // Parameters
            const int minCharge = 14;
            const int maxCharge = 14;

            // Configure amino acids
            var aaSet = new AminoAcidSet();

            const string protAnnotation = "A.HAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK.-";
            const int scanNum = 810;

            // Create a sequence graph
            //var protSeq = protAnnotation.Substring(2, protAnnotation.Length - 4);
            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            if (seqGraph == null)
            {
                Console.WriteLine("Invalid sequence: {0}", protAnnotation);
                return;
            }

            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);

            var precursorFilter = new PrecursorFilter(run, new Tolerance(10));

            foreach (var seqComposition in seqGraph.GetSequenceCompositions())
            {
                Console.WriteLine("Composition: {0}, Mass: {1}", seqComposition, seqComposition.GetMass());

                for (var charge = minCharge; charge <= maxCharge; charge++)
                {
                    var precursorIon = new Ion(seqComposition + Composition.H2O, charge);
                    Console.WriteLine("Charge: {0}, m/z: {1}", charge, precursorIon.GetMz());

                    foreach (var isotope in precursorIon.GetIsotopes(relativeIntensityThreshold: 0.1))
                    {
                        Console.WriteLine("{0}\t{1}\t{2}", isotope.Item1, precursorIon.GetIsotopeMz(isotope.Item1), isotope.Item2);
                    }

                    var isValid = precursorFilter.IsValid(precursorIon, scanNum);
                    Console.WriteLine("IsValid: {0}", isValid);

                    foreach (var ms2ScanNum in run.GetFragmentationSpectraScanNums(precursorIon))
                    {
                        Console.WriteLine("ScanNum: {0}", ms2ScanNum);
                    }
                }
            }
        }

        [Test]
        public void TestMs1BasedFilteringOneProtein()
        {
            // Search parameters
            const int minCharge = 9;
            const int maxCharge = 19; // target charge: 14
            var precursorTolerance = new Tolerance(20);

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            Console.Write("Reading raw file...");
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);

            var aaSet = new AminoAcidSet();

            var sequences = new[] { "A.HAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK.-" };

            long totalProtCompositions = 0;
            long numPrecursorIons = 0;

            sw.Reset();
            sw.Start();

            var precursorFilter = new PrecursorFilter(run, new Tolerance(10));

            var outputFilePath = Path.ChangeExtension(specFilePath, ".icresult");
            using (var writer = new StreamWriter(outputFilePath))
            {
                foreach (var protAnnotation in sequences)
                {
                    var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
                    if (seqGraph == null) continue;

                    var protCompositions = seqGraph.GetSequenceCompositionsWithNTermCleavage(0);
                    foreach (var protComposition in protCompositions)
                    {
                        totalProtCompositions++;
                        for (var charge = minCharge; charge <= maxCharge; charge++)
                        {
                            numPrecursorIons++;
                            var precursorIon = new Ion(protComposition + Composition.H2O, charge);
                            foreach (var ms2ScanNum in run.GetFragmentationSpectraScanNums(precursorIon))
                            {
                                if (!precursorFilter.IsValid(precursorIon, ms2ScanNum)) continue;
                                writer.WriteLine("{0}\t{1}\t{2}\t{3}", protAnnotation, charge, precursorIon.GetMz(), ms2ScanNum);
                            }
                        }
                    }
                }
            }

            sw.Stop();
            Console.WriteLine("NumProteinCompositions: {0}", totalProtCompositions);
            Console.WriteLine("NumPrecursorIons: {0}", numPrecursorIons);
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void CountProteins()
        {
            const string filePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.icresult";

            var protSet = new HashSet<string>();
            foreach (var line in File.ReadLines(filePath))
            {
                var token = line.Split('\t');
                if (token.Length != 4) continue;
                protSet.Add(token[0]);
            }
            Console.WriteLine("NumProteins: {0}", protSet.Count);
        }

        [Test]
        public void TestMs1BasedFiltering()
        {
            // Search parameters
            const int numNTermCleavages = 0;  // 30
            const int minLength = 7;
            const int maxLength = 1000;
            const int minCharge = 3; // 3
            const int maxCharge = 67; // 67
            const int numMaxModsPerProtein = 0; // 6
            var precursorTolerance = new Tolerance(20);
            //const string dbFilePath = @"..\..\..\TestFiles\BSA.fasta";
            const string dbFilePath =
                @"C:\cygwin\home\kims336\Data\TopDown\ID_003558_56D73071.fasta";

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            Console.Write("Reading raw file...");
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
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
                    //dehydro,
                    //cysteinylC,
                    //glutathioneC,
                    //oxM
                };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            //var targetDb = new FastaDatabase(dbFilePath);
            //var indexedDb = new IndexedDatabase(targetDb);
            //var sequences = indexedDb.SequencesAsStrings(0, minLength, maxLength);

            var sequences = new[] { "A.HAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK.-" };

            var numProteins = 0;
            long totalProtCompositions = 0;
            long numPrecursorIons = 0;

            sw.Reset();
            sw.Start();

            var precursorFilter = new PrecursorFilter(run, new Tolerance(10));

            var outputFilePath = Path.ChangeExtension(specFilePath, ".icresult");
            using (var writer = new StreamWriter(outputFilePath))
            {
                foreach (var protAnnotation in sequences)
                {
                    ++numProteins;

                    if (numProteins % 1000 == 0)
                    {
                        Console.WriteLine("Processed {0} proteins", numProteins);
                    }

                    //Console.WriteLine(protAnnotation);

                    var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
                    if (seqGraph == null) continue;

                    for (var nTermCleavages = 0; nTermCleavages <= numNTermCleavages; nTermCleavages++)
                    {
                        var protCompositions = seqGraph.GetSequenceCompositionsWithNTermCleavage(nTermCleavages);
                        foreach (var protComposition in protCompositions)
                        {
                            totalProtCompositions++;
                            for (var charge = minCharge; charge <= maxCharge; charge++)
                            {
                                numPrecursorIons++;
                                var precursorIon = new Ion(protComposition + Composition.H2O, charge);
                                foreach (var ms2ScanNum in run.GetFragmentationSpectraScanNums(precursorIon))
                                {
                                    if (!precursorFilter.IsValid(precursorIon, ms2ScanNum)) continue;
                                    writer.WriteLine("{0}\t{1}\t{2}\t{3}", protAnnotation, charge, precursorIon.GetMz(), ms2ScanNum);
                                }
                            }
                        }
                    }
                }                
            }

            sw.Stop();
            Console.WriteLine("NumProteins: {0}", numProteins);
            Console.WriteLine("NumProteinCompositions: {0}", totalProtCompositions);
            Console.WriteLine("NumPrecursorIons: {0}", numPrecursorIons);
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void TestGeneratingAllXics()
        {
            // Search parameters
            const int numNTermCleavages = 1;  // 30
            const int minLength = 7;
            const int maxLength = 1000;
            const int minCharge = 3; // 3
            const int maxCharge = 67; // 67
            const int numMaxModsPerProtein = 0; // 6
            var precursorTolerance = new Tolerance(20);
            //const string dbFilePath = @"..\..\..\TestFiles\BSA.fasta";
            const string dbFilePath =
                @"C:\cygwin\home\kims336\Data\TopDown\ID_003558_56D73071.fasta";

            var sw = new System.Diagnostics.Stopwatch();

            sw.Start();
            Console.Write("Reading raw file...");
            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
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
            var indexedDb = new IndexedDatabase(targetDb);

            var numProteins = 0;
            long totalProtCompositions = 0;
            long numXics = 0;

            sw.Reset();
            sw.Start();
            Console.WriteLine("Generating XICs...");

            foreach (var protAnnotation in indexedDb.SequencesAsStrings(0, minLength, maxLength))
            {
                ++numProteins;

                if (numProteins % 1000 == 0)
                {
                    Console.WriteLine("Processed {0} proteins", numProteins);
                }

                //Console.WriteLine(protAnnotation);

                var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
                if (seqGraph == null) continue;

                for (var nTermCleavages = 0; nTermCleavages <= numNTermCleavages; nTermCleavages++)
                {
                    var protCompositions = seqGraph.GetSequenceCompositionsWithNTermCleavage(nTermCleavages);
                    foreach (var protComposition in protCompositions)
                    {
                        totalProtCompositions++;
                        var mostAbundantIsotopeIndex = protComposition.GetMostAbundantIsotopeZeroBasedIndex();
                        for (var charge = minCharge; charge <= maxCharge; charge++)
                        {
                            numXics++;
                            var precursorIon = new Ion(protComposition + Composition.H2O, charge);
                            run.GetExtractedIonChromatogram(precursorIon.GetIsotopeMz(mostAbundantIsotopeIndex), precursorTolerance);
                        }
                    }
                }
            }

            sw.Stop();
            Console.WriteLine("NumProteins: {0}", numProteins);
            Console.WriteLine("NumProteinCompositions: {0}", totalProtCompositions);
            Console.WriteLine("NumXics: {0}", numXics);
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }
        
        [Test]
        public void TestGeneratingXics()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            const string specFilePath = @"C:\workspace\TopDown\E_coli_iscU_60_mock.raw";
            var run = new LcMsRun(new XCaliburReader(specFilePath));
            const string protAnnotation = "A.HAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK.-";
            var aaSet = new AminoAcidSet();

            var precursorTolerance = new Tolerance(10);

            // Create a sequence graph
            var protSeq = protAnnotation.Substring(2, protAnnotation.Length - 4);
            var seqGraph = SequenceGraph.CreateGraph(aaSet, protSeq);
            foreach (var protComposition in seqGraph.GetSequenceCompositions())
            {
                var mostAbundantIsotopeIndex = protComposition.GetMostAbundantIsotopeZeroBasedIndex();
                Console.WriteLine("Composition\t{0}", protComposition);
                Console.WriteLine("MostAbundantIsotopeIndex\t{0}", mostAbundantIsotopeIndex);
                Console.WriteLine();

                for (var charge = 10; charge <= 14; charge++)
                {
                    var precursorIon = new Ion(protComposition+Composition.H2O, charge);
                    var xic = run.GetExtractedIonChromatogram(precursorIon.GetIsotopeMz(mostAbundantIsotopeIndex), precursorTolerance);
                    //Console.WriteLine(xic[0].ScanNum + " " + xic[1].ScanNum);
                    
                    Console.WriteLine("ScanNum\t{0}", string.Join("\t", xic.Select(p => p.ScanNum.ToString())));
                    Console.WriteLine("Charge " + charge + "\t" + string.Join("\t", xic.Select(p => p.Intensity.ToString())));
                }

                Console.WriteLine("\nCharge\tm/z");
                for (var charge = 9; charge <= 18; charge++)
                {
                    var precursorIon = new Ion(protComposition + Composition.H2O, charge);
                    Console.WriteLine("{0}\t{1}", charge, precursorIon.GetIsotopeMz(mostAbundantIsotopeIndex));
                }
            }

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }



        [Test]
        public void TestHistonEnumeration()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            const int numNTermCleavages = 0;
            const int numMaxModsPerProtein = 11;
            //const string protAnnotation = "-.MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG.-";  // Histone H4
            const string protAnnotation =
                "-.MARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA.-"; // Histone H3.2

            var acetylR = new SearchModification(Modification.Acetylation, 'R', SequenceLocation.Everywhere, false);
            var acetylK = new SearchModification(Modification.Acetylation, 'K', SequenceLocation.Everywhere, false);
            var methylR = new SearchModification(Modification.Methylation, 'R', SequenceLocation.Everywhere, false);
            var methylK = new SearchModification(Modification.Methylation, 'K', SequenceLocation.Everywhere, false);
            var diMethylR = new SearchModification(Modification.DiMethylation, 'R', SequenceLocation.Everywhere, false);
            var diMethylK = new SearchModification(Modification.DiMethylation, 'K', SequenceLocation.Everywhere, false);
            var triMethylR = new SearchModification(Modification.TriMethylation, 'R', SequenceLocation.Everywhere, false);
            var phosphoS = new SearchModification(Modification.Phosphorylation, 'S', SequenceLocation.Everywhere, false);
            var phosphoT = new SearchModification(Modification.Phosphorylation, 'T', SequenceLocation.Everywhere, false);
            var phosphoY = new SearchModification(Modification.Phosphorylation, 'Y', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
                {
                    acetylR, acetylK,
                    methylR, methylK,
                    diMethylR, diMethylK,
                    triMethylR,
                    phosphoS, phosphoT, phosphoY
                };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var protSeq = protAnnotation.Substring(2, protAnnotation.Length - 4);
            var seqGraph = new SequenceGraph(aaSet, protSeq.Length);
            seqGraph.AddAminoAcid(AminoAcid.ProteinNTerm.Residue);
            var isValidSequence = protSeq.All(seqGraph.AddAminoAcid);
            if (!isValidSequence)
            {
                Console.WriteLine("Invalid sequence: {0}", protSeq);
                return;
            }
            seqGraph.AddAminoAcid(AminoAcid.ProteinCTerm.Residue);

            Console.WriteLine("Num sequence compositions: {0}, {1}", seqGraph.GetNumCompositions(), seqGraph.GetNumDistinctCompositions()
                );

            Console.WriteLine("Num product compositions: {0}", seqGraph.GetNumProductCompositions()
                );

            //var numProductCompositions = 0L;

            //foreach (var seqComposition in seqGraph.GetSequenceCompositions())
            //{
            //    Console.WriteLine("Composition: {0}, Mass: {1}", seqComposition, seqComposition.GetMass());
            //    foreach (var scoringGraph in seqGraph.GetScoringGraphs())
            //    {
            //        //numProductCompositions += scoringGraph.GetCompositions().Count();
            //    }
            //}
            //Console.WriteLine("Num product compositions: {0}", numProductCompositions
            //    );

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }

        [Test]
        public void TestTopDownEnumeration()
        {
            // Search parameters
            const int numNTermCleavages = 0;
            const int minLength = 7;
            const int maxLength = 1000;
            const int numMaxModsPerProtein = 6;
            //const string dbFilePath = @"..\..\..\TestFiles\BSA.fasta";
            const string dbFilePath =
                @"..\..\..\TestFiles\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

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
            var indexedDb = new IndexedDatabase(targetDb);
            var maxProteinLength = indexedDb.GetLongestSequenceLength();

            var numProteins = 0;
            long totalProtCompositions = 0;
            long totalProductCompositions = 0;
            foreach (var protAnnotation in indexedDb.SequencesAsStrings(numNTermCleavages, minLength, maxLength))
            {
                ++numProteins;

                if (numProteins % 1000 == 0)
                {
                    Console.WriteLine("Processed {0} proteins", numProteins);
                }
                //if (numProteins > 100)
                //    break;

                //if (numProteins < 28367) continue;
                //Console.WriteLine("{0} {1}", protAnnotation.Length, protAnnotation);

                var protSeq = protAnnotation.Substring(2, protAnnotation.Length - 4);
                var nTerm = protAnnotation[0] == FastaDatabase.Delimiter
                                      ? AminoAcid.ProteinNTerm
                                      : AminoAcid.PeptideNTerm;
                var cTerm = protAnnotation[protAnnotation.Length-1] == FastaDatabase.Delimiter
                                      ? AminoAcid.ProteinCTerm
                                      : AminoAcid.PeptideCTerm;

                var seqGraph = new SequenceGraph(aaSet, protSeq.Length);
                seqGraph.AddAminoAcid(nTerm.Residue);
                var isValidSequence = protSeq.All(seqGraph.AddAminoAcid);

                if (!isValidSequence) continue;

                seqGraph.AddAminoAcid(cTerm.Residue);

                //var numProtCompositions = seqGraph.GetSequenceCompositions().Length;
                var numProtCompositions = seqGraph.GetNumCompositions();
                totalProtCompositions += numProtCompositions;

                long numProductCompositions = 0;
                //foreach (var scoringGraph in seqGraph.GetScoringGraphs())
                //{
                //    //numProductCompositions += scoringGraph.GetCompositions().Count();
                //}
                //totalProductCompositions += numProductCompositions;
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}",
                    numProteins, numProtCompositions, numProductCompositions,
                    totalProtCompositions, totalProductCompositions);
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
            const int maxLength = 1000;
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\SuffixArray\BSA.fasta";

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            var targetDb = new FastaDatabase(dbFilePath);
            var indexedDb = new IndexedDatabase(targetDb);

            var protIndex = 0;
            long numProtCompositions = 0;
            foreach (var protAnnotation in indexedDb.SequencesAsStrings(numMissedCleavages, minLength, maxLength))
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
