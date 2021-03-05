using System;
using System.Collections.Generic;
using System.Reflection;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Tests.Base;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Tests.FunctionalTests
{
    public class TestInformedTopDownScoring
    {
        [OneTimeSetUp]
        public void Setup()
        {
            // Verify that the test .pbf file exists
            // If it does not exist, yet the .mzML file exists, create the .pbf file
            Utils.GetPbfTestFilePath(true);
        }

        [Test]
        [TestCase(178, 6, "MLILTRRVGETLMIGDEVTVTVLGVKGNQVRIGVNAPKEVSVHREEIYQRIQSEKS", 114.36886559)]
        [TestCase(266, 15, "ALRLKDLVKKTERQLSDYQRQLSMVKTTESVQKATATITDSFASSNSKLLNAKDSLERIKARQQQFDDRLKAAETLAEEGSDKSLQAKLAEAGIGEQKSNANAVLERIKARKS", 63.5591347)]
        public void TestRescoring(int scanNum, int charge, string sequence, double expectedScore)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var pbfFilePath = Utils.GetPbfTestFilePath(false);
            var pbfFile = Utils.GetTestFile(methodName, pbfFilePath);

            // Configure amino acid set
            var acetylN = new SearchModification(Modification.Acetylation, '*', SequenceLocation.ProteinNTerm, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                dehydroC,
                glutathioneC,
                oxM,
                acetylN,
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);

            var composition = aaSet.GetComposition(sequence) + Composition.H2O;

            var run = PbfLcMsRun.GetLcMsRun(pbfFile.FullName, 0, 0);
            var informedScorer = new InformedTopDownScorer(run, aaSet, 1, 15, new Tolerance(10));
            var scores = informedScorer.GetScores(AminoAcid.ProteinNTerm, sequence, AminoAcid.ProteinCTerm, composition, charge, scanNum);
            Console.WriteLine("Total Score = " + scores.Score);
            Console.WriteLine("#Fragments = " + scores.NumMatchedFrags);

            Assert.AreEqual(expectedScore, scores.Score, 0.001);
        }
    }
}
