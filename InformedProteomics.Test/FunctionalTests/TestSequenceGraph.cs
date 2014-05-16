using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestSequenceGraph
    {
        [Test]
        public void TestBuildingReverseGraph()
        {
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

            //for (var modIndex = 0; modIndex < seqCompositions.Length; modIndex++)
            const int modIndex = 4;
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
            const string modFilePath = @"\\protoapps\UserData\Sangtae\TestData\MiscFiles\Mods.txt";
            var modFileParser = new ModFileParser(modFilePath);
            Console.WriteLine("MaxNumDynModsPerSequence: {0}", modFileParser.MaxNumDynModsPerSequence);
            foreach (var searhMod in modFileParser.SearchModifications)
            {
                Console.WriteLine("{0}\t{1}\t{2}\t{3}", searhMod.TargetResidue, searhMod.Location, searhMod.IsFixedModification, searhMod.Modification);
            }
        }
    }
}
