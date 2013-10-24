using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using NUnit.Framework;

namespace InformedProteomics.Test
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
            foreach (var composition in seqGraph.GetSequenceCompositionsWithNTermCleavage(1))
            {
                Console.WriteLine("{0}\t{1}", composition, composition.GetMass());
            }
        }

        [Test]
        public void TestBuildingSequenceGraphLongProtein()
        {
            // Configure amino acid set
            const int numMaxModsPerProtein = 6;
            var pyroGluQ = new SearchModification(Modification.PyroGluQ, 'Q', SequenceLocation.ProteinNTerm, false);
            var dehydro = new SearchModification(Modification.PyroGluQ, 'C', SequenceLocation.Everywhere, false);
            var cysteinylC = new SearchModification(Modification.CysteinylC, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.GlutathioneC, 'C', SequenceLocation.Everywhere, false);
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
            const string protAnnotation =
                "_.QQ._";

            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            var seqCompositions = seqGraph.GetSequenceCompositionsWithNTermCleavage(0);
            Console.WriteLine("SequenceCompositions: {0}", seqCompositions.Length);
            foreach (var composition in seqCompositions)
            {
                Console.WriteLine("{0}\t{1}", composition, composition.GetMass());
            }

            Console.WriteLine("FragmentCompositions");
            var numNodes = 0;
            foreach (var composition in seqGraph.GetAllFragmentNodeCompositions())
            {
                numNodes++;

                //if (composition.GetMass() < 0)
                {
                    Console.WriteLine(composition);
                }
            }
            Console.WriteLine("NumNodes: {0}", numNodes);
        }
    }
}
