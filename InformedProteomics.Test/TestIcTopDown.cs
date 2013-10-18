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
            const int minCharge = 3;
            const int maxCharge = 67;
            var precursorIonTolerance = new Tolerance(10);
            var productIonTolerance = new Tolerance(10);

            // Configure amino acids
            var aaSet = new AminoAcidSet();

            const string protAnnotation = "A.HAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK.-";
            const int scanNum = 810;

            // Create a sequence graph
            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            if (seqGraph == null)
            {
                Console.WriteLine("Invalid sequence: {0}", protAnnotation);
                return;
            }

            const string specFilePath = @"C:\cygwin\home\kims336\Data\TopDown\E_coli_iscU_60_mock.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);

            var precursorFilter = new PrecursorFilter(run, precursorIonTolerance);

            var seqCompositionArr = seqGraph.GetSequenceCompositions();
            for (var modIndex=0; modIndex<seqCompositionArr.Length; modIndex++)
            {
                var seqComposition = seqCompositionArr[modIndex];
                Console.WriteLine("Composition: {0}, Mass: {1}", seqComposition, seqComposition.GetMass());

                for (var charge = minCharge; charge <= maxCharge; charge++)
                {
                    var precursorIon = new Ion(seqComposition + Composition.H2O, charge);
                    //Console.WriteLine("Charge: {0}, m/z: {1}", charge, precursorIon.GetMz());

                    //foreach (var isotope in precursorIon.GetIsotopes(relativeIntensityThreshold: 0.1))
                    //{
                    //    Console.WriteLine("{0}\t{1}\t{2}", isotope.Item1, precursorIon.GetIsotopeMz(isotope.Item1), isotope.Item2);
                    //}

                    var bestScore = Double.NegativeInfinity;
                    var bestScanNum = -1;
                    foreach (var ms2ScanNum in run.GetFragmentationSpectraScanNums(precursorIon))
                    {
                        if (!precursorFilter.IsValid(precursorIon, ms2ScanNum)) continue;


                        var spec = run.GetSpectrum(ms2ScanNum);
                        spec.FilterNoise(signalToNoiseRatio: 1.5);
                        var scorer = new MatchedPeakCounter(spec, productIonTolerance, 1, 6);
                        var score = seqGraph.GetScore(modIndex, 0, charge, scorer);

                        if (score > bestScore)
                        {
                            bestScore = score;
                            bestScanNum = ms2ScanNum;
                        }
                    }
                    if (bestScore > 10)
                    {
                        Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", protAnnotation, charge, precursorIon.GetMz(), bestScanNum, bestScore);
                    }
                }
            }
        }

        [Test]
        public void TestTopDownSearch()
        {
            // Search parameters
            const int maxNumNTermCleavages = 1;  // 30
            const int minLength = 7;
            const int maxLength = 1000;
            const int minCharge = 3; // 3
            const int maxCharge = 67; // 67
            const int numMaxModsPerProtein = 0; // 6

            var precursorTolerance = new Tolerance(10);
            var productIonTolerance = new Tolerance(10);

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

            var targetDb = new FastaDatabase(dbFilePath);
            var indexedDb = new IndexedDatabase(targetDb);
            var sequences = indexedDb.SequencesAsStrings(minLength, maxLength);

            //var sequences = new[] { "A.HAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK.-" };

            var numProteins = 0;
            long totalProtCompositions = 0;
            long numPrecursorIons = 0;

            sw.Reset();
            sw.Start();

            var precursorFilter = new PrecursorFilter(run, precursorTolerance);

            var outputFilePath = Path.ChangeExtension(specFilePath, ".icresult");
            using (var writer = new StreamWriter(outputFilePath))
            {
                foreach (var protAnnotation in sequences)
                {
                    ++numProteins;

                    //if (numProteins % 100 == 0)
                    {
                        Console.WriteLine("Processed {0} proteins", numProteins);
                    }

                    //Console.WriteLine(protAnnotation);

                    var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
                    if (seqGraph == null) continue;

                    for (var numNTermCleavage = 0; numNTermCleavage <= maxNumNTermCleavages; numNTermCleavage++)
                    {
                        var protCompositions = seqGraph.GetSequenceCompositionsWithNTermCleavage(numNTermCleavage);
                        for(var modIndex=0; modIndex<protCompositions.Length; modIndex++)
                        //foreach (var protComposition in protCompositions)
                        {
                            var protComposition = protCompositions[modIndex];
                            totalProtCompositions++;
                            for (var charge = minCharge; charge <= maxCharge; charge++)
                            {
                                numPrecursorIons++;
                                var precursorIon = new Ion(protComposition + Composition.H2O, charge);

                                var bestScore = Double.NegativeInfinity;
                                var bestScanNum = -1;
                                foreach (var ms2ScanNum in run.GetFragmentationSpectraScanNums(precursorIon))
                                {
                                    if (!precursorFilter.IsValid(precursorIon, ms2ScanNum)) continue;
                                    

                                    var spec = run.GetSpectrum(ms2ScanNum);
                                    spec.FilterNoise(signalToNoiseRatio: 1.5);
                                    var scorer = new MatchedPeakCounter(spec, productIonTolerance, 1, 6);
                                    var score = seqGraph.GetScore(modIndex, 0, charge, scorer);

                                    if (score > bestScore)
                                    {
                                        bestScore = score;
                                        bestScanNum = ms2ScanNum;
                                    }
                                }
                                if (bestScore > 10)
                                {
                                    writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", protAnnotation.Substring(numNTermCleavage + 2), charge, precursorIon.GetMz(), bestScanNum, bestScore);
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

            foreach (var protAnnotation in indexedDb.SequencesAsStrings(minLength, maxLength))
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

            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            if (seqGraph == null)
            {
                Console.WriteLine("Invalid sequence: {0}", protAnnotation);
                return;
            }

            Console.WriteLine("Num sequence compositions: {0}, {1}", seqGraph.GetNumSequenceCompositions(), seqGraph.GetNumDistinctSequenceCompositions()
                );

            Console.WriteLine("Num product compositions: {0}", seqGraph.GetNumFragmentCompositions()
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
            const int numNTermCleavages = 30;
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

            var numProteins = 0;
            long totalProtCompositions = 0;
            foreach (var protAnnotation in indexedDb.SequencesAsStrings(minLength, maxLength))
            {
                ++numProteins;

                if (numProteins % 1000 == 0)
                {
                    Console.WriteLine("Processed {0} proteins", numProteins);
                }

                var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
                if (seqGraph == null) continue;

                for (var nTermCleavage = 0; nTermCleavage <= numNTermCleavages; nTermCleavage++)
                {
                    totalProtCompositions += seqGraph.GetNumSequenceCompositionsWithNTermCleavage(nTermCleavage);
                }
            }

            sw.Stop();
            Console.WriteLine("NumProteins: {0}", numProteins);
            Console.WriteLine("NumProteinCompositions: {0}", totalProtCompositions);
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Elapsed Time: {0:f4} sec", sec);
        }
    }
}
