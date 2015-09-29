using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Scoring.BottomUp;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestSequenceGraph
    {
        [Test]
        public void TestBuildingReverseGraph()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string annotation = "_.MARTKQTARK._";

            // Configure amino acid set
            var methylK = new SearchModification(Modification.Methylation, 'K', SequenceLocation.Everywhere, false);
            //var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                methylK,
                //pyroGluQ,
                oxM
            };

            const int numMaxModsPerProtein = 2;

            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var seqGraph = SequenceGraph.CreateGraph(aaSet, annotation);
            foreach (var composition in seqGraph.GetSequenceCompositions())
            {
                Console.WriteLine("{0}\t{1}", composition, composition.Mass);
            }
        }

        [Test]
        public void TestBuildingSequenceGraphLongProtein()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            // Configure amino acid set
            const int numMaxModsPerProtein = 6;
            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.ProteinNTerm, false);
            var dehydro = new SearchModification(Modification.PyroGluQ, 'C', SequenceLocation.Everywhere, false);
            var cysteinylC = new SearchModification(Modification.Cysteinyl, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                pyroGluQ,
                //dehydro,
                //cysteinylC,
                //glutathioneC,
                //oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            //const string protAnnotation = "A.HAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK.-";
            //const string protAnnotation =
            //    "_.QQ._";

            const string protAnnotation =
                "_.MKLYNLKDHNEQVSFAQAVTQGLGKNQGLFFPHDLPEFSLTEIDEMLKLDFVTRSAKILSAFIGDEIPQEILEERVRAAFAFPAPVANVESDVGCLELFHGPTLAFKDFGGRFMAQMLTHIAGDKPVTILTATSGDTGAAVAHAFYGLPNVKVVILYPRGKISPLQEKLFCTLGGNIETVAIDGDFDACQALVKQAFDDEELKVALGLNSANSINISRLLAQICYYFEAVAQLPQETRNQLVVSVPSGNFGDLTAGLLAKSLGLPVKRFIAATNVNDTVPRFLHDGQWSPKATQATLSNAMDVSQPNNWPRVEELFRRKIWQLKELGYAAVDDETTQQTMRELKELGYTSEPHAAVAYRALRDQLNPGEYGLFLGTAHPAKFKESVEAILGETLDLPKELAERADLPLLSHNLPADFAALRKLMMNHQ._";

            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            var seqCompositions = seqGraph.GetSequenceCompositions();

            for (var modIndex = 0; modIndex < seqCompositions.Length; modIndex++)
            {
                var seqComposition = seqCompositions[modIndex];
                Console.WriteLine("SequenceComposition: {0}", seqComposition);

                foreach (var composition in seqGraph.GetFragmentCompositions(modIndex, 0))
                {
                    //if (composition.GetMass() > seqComposition.GetMass())
                    {
                        Console.WriteLine("***Seq: {0}, Frag: {1}", seqComposition, composition);
                    }
                }
            }
        }

        [Test]
        public void TestGraphWithModifications()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string annotation = "_.MIALNKTPQTIVFYKPYGVLCQFTDNSAHPRPTLKDYINLPDLYPVGRLDQDSEGLLLLTSNGKLQHRLAHREFAHQRTYFAQVEGSPTDEDLEPLRRGITFADYPTRPAIAKIITEPDFPPRNPPIRYRASIPTSWLSITLTEGRNRQVRRMTAAVGFPTLRLVRVQIQVTGRSPQQGKGKSAATWCLTLEGLSPGQWRPLTPWEENFCQQLLTGNPNGPWQKKFGDRR._";

            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var nitrosylC = new SearchModification(Modification.Nitrosyl, 'C', SequenceLocation.Everywhere, false);
            var nethylmaleimideC = new SearchModification(Modification.Nethylmaleimide, 'C', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                dehydroC,
                glutathioneC,
                nitrosylC,
                nethylmaleimideC,
                oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var seqGraph = SequenceGraph.CreateGraph(aaSet, annotation);
            var seqCompositions = seqGraph.GetSequenceCompositions();
            var modCombs = seqGraph.GetModificationCombinations();

            Console.WriteLine("*** Before cleavage: {0}", seqCompositions.Length);
            for (var modIndex = 0; modIndex < seqCompositions.Length; modIndex++)
            {
                var seqComposition = seqCompositions[modIndex];
                Console.WriteLine("SequenceComposition: {0}, ModComb: {1}", seqComposition, modCombs[modIndex]);
            }

            seqGraph.CleaveNTerm();
            seqCompositions = seqGraph.GetSequenceCompositions();
            modCombs = seqGraph.GetModificationCombinations();
            Console.WriteLine("*** After cleavage: {0}", seqCompositions.Length);
            for (var modIndex = 0; modIndex < seqCompositions.Length; modIndex++)
            {
                var seqComposition = seqCompositions[modIndex];
                Console.WriteLine("SequenceComposition: {0}, ModComb: {1}", seqComposition, modCombs[modIndex]);
            }
        }

        [Test]
        public void TestCreatingAminoAcidSet()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            // Configure amino acid set
            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.PeptideNTerm, false);
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                acetylN,
                pyroGluQ,
                oxM
            };

            const int numMaxModsPerProtein = 2;

            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            //var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            aaSet.Display();
        }

        [Test]
        public void TestNTermMods()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string annotation = "_.QARTKQTARK._";

            // Configure amino acid set
            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.ProteinNTerm, false);
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                acetylN,
                pyroGluQ,
                //oxM
            };

            const int numMaxModsPerProtein = 2;

            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            //aaSet.Display();
            var seqGraph = SequenceGraph.CreateGraph(aaSet, annotation);
            foreach (var composition in seqGraph.GetSequenceCompositions())
            {
                Console.WriteLine("{0}\t{1}", composition, composition.Mass);
            }

            Console.WriteLine("*** Cleave N-term");
            seqGraph.CleaveNTerm();
            foreach (var composition in seqGraph.GetSequenceCompositions())
            {
                Console.WriteLine("{0}\t{1}", composition, composition.Mass);
            }
        }

        [Test]
        public void TestReadingModFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string modFilePath = @"..\..\..\TestFiles\Mods.txt";

            if (!File.Exists(modFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, modFilePath);
            }
            
            var modFileParser = new ModFileParser(modFilePath);
            Console.WriteLine("MaxNumDynModsPerSequence: {0}", modFileParser.MaxNumDynModsPerSequence);
            foreach (var searhMod in modFileParser.SearchModifications)
            {
                Console.WriteLine("{0}\t{1}\t{2}\t{3}", searhMod.TargetResidue, searhMod.Location, searhMod.IsFixedModification, searhMod.Modification);
            }
            var aaSet = new AminoAcidSet(modFilePath);
            aaSet.Display();
        }

        [Test]
        public void TestParsingCompositions()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string compStr = "CH-1N-3O23S";
            var composition = Composition.ParseFromPlainString(compStr);
            Console.WriteLine(composition);
        }

        [Test]
        public void TestNumberOfProteoforms()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            //const string annotation = "_.AMCMC._";
            const string annotation = "_.MARTKQTARK._";

            // Configure amino acid set
            var methylK = new SearchModification(Modification.Methylation, 'K', SequenceLocation.Everywhere, false);
            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var carbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C', SequenceLocation.Everywhere, true);
            var methylC = new SearchModification(Modification.Methylation, 'C', SequenceLocation.Everywhere, true);
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.PeptideNTerm, false);

            var searchModifications = new List<SearchModification>
            {
                //carbamidomethylC,
                //methylC,
                methylK,
                //pyroGluQ,
                oxM,
                //acetylN
            };

            const int numMaxModsPerProtein = 4;

            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var seqGraph = SequenceGraph.CreateGraph(aaSet, annotation);
            var protCompositions = seqGraph.GetSequenceCompositions();
            var modCombs = seqGraph.GetModificationCombinations();

            Console.WriteLine("\n#Protoeoforms by mod combinations: ");
            for (var modIndex = 0; modIndex < protCompositions.Length; modIndex++)
            {
                Console.Write((modIndex == 0) ? "No modifications" : modCombs[modIndex].ToString());
                Console.Write("\t");
                Console.WriteLine("{0}", seqGraph.GetNumProteoformSequences(modIndex));
            }

            Console.WriteLine("\n#Protoeoforms by number of modificaionts: ");
            for (var nMod = 0; nMod <= numMaxModsPerProtein; nMod++)
            {
                Console.Write("#modificaitons = {0}", nMod);
                Console.Write("\t");
                Console.WriteLine("{0}", seqGraph.GetNumProteoformSequencesByNumMods(nMod));
            }
        }

        [Test]
        public void TestGettingSequence()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            //const string annotation = "_.AMCMC._";
            const string annotation = "_.MARTKQTARK._";

            // Configure amino acid set
            var methylK = new SearchModification(Modification.Methylation, 'K', SequenceLocation.Everywhere, false);
            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var carbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C', SequenceLocation.Everywhere, true);
            var methylC = new SearchModification(Modification.Methylation, 'C', SequenceLocation.Everywhere, true);
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.PeptideNTerm, false);

            var searchModifications = new List<SearchModification>
            {
                //carbamidomethylC,
                //methylC,
                methylK,
                //pyroGluQ,
                oxM,
                //acetylN
            };

            const int numMaxModsPerProtein = 2;

            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var seqGraph = SequenceGraph.CreateGraph(aaSet, annotation);
            var protCompositions = seqGraph.GetSequenceCompositions();
            
            for (var modIndex = 0; modIndex < protCompositions.Length; modIndex++)
            {
                seqGraph.SetSink(modIndex);

                var composition = protCompositions[modIndex];
                Console.WriteLine("{0}\t{1}", composition, composition.Mass);
                var curScoreAndModifications = seqGraph.GetFragmentScoreAndModifications(new DummyScorer());
                if (curScoreAndModifications != null) Console.WriteLine("Score: {0}, Modifications: {1}", curScoreAndModifications.Item1, curScoreAndModifications.Item2);
            }
        }

        [Test]
        public void TestCreatingHistoneGraph()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const int numMaxModsPerProtein = 11;

            // Histone H4
            const string annotation =
                "_.MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG._";

            // Histone H3.1
//            const string annotation =
//                "_.MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA._";

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
                acetylR,
                acetylK,
                methylR,
                methylK,
                diMethylR,
                diMethylK,
                triMethylR,
                phosphoS,
                phosphoT,
                phosphoY
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);
            var graph = SequenceGraph.CreateGraph(aaSet, annotation);

            var numFragCompositions = graph.GetNumFragmentCompositions();
            var numProteoforms = graph.GetNumProteoformCompositions();
            var numSeqCompositions = graph.GetNumProteoformCompositions();

            Console.WriteLine("NumFragmentCompositions: " + numFragCompositions);
            Console.WriteLine("NumProteoforms: " + numProteoforms);
            Console.WriteLine("NumSequenceCompositions: " + numSeqCompositions);
        }
    }


}
