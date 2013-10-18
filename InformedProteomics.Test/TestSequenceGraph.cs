using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
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
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                methylK,
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
        public void SimpleTest()
        {
            const string protAnnotation = "A.HAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK.-";

            var aaSet = new AminoAcidSet();
            // Create a sequence graph
            var seqGraph = SequenceGraph.CreateGraph(aaSet, protAnnotation);
            foreach (var composition in seqGraph.GetSequenceCompositions())
            {
                var ion = new Ion(composition + Composition.H2O, 11);
                Console.WriteLine("{0}\t{1}\t{2}", ion.GetMz(), ion.Composition, ion.Composition.GetMass());
            }
        }
    }
}
