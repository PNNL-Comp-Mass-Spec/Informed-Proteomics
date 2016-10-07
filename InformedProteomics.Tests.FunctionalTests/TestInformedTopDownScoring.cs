using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Tests.FunctionalTests
{
    public class TestInformedTopDownScoring
    {
        [Test]
        public void TestRescoring()
        {
            //const string specFilePath = @"H:\Research\QCShew_TopDown\Production\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
            const string specFilePath = @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf";
            //const string sequence = "SGWYELSKSSNDQFKFVLKAGNGEVILTSELYTGKSGAMNGIESVQTNSPIEARYAKEVAKNDKPYFNLKAANHQIIGTSQMYSSTA";
            //const int scanNum = 4084;

            const string sequence = "SKTKHPLPEQWQKNQEAAKATQVAFDLDEKFQYSIRKAALDAGVSPSDQIRTILGLSVSRRPTRPRLTVSLNADDYVQLAEKYDLNADAQLEIKRRVLEDLVRFVAED";
            const int scanNum = 5448;
            const int charge = 11;

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

            var run = PbfLcMsRun.GetLcMsRun(specFilePath, 0, 0);
            var informedScorer = new InformedTopDownScorer(run, aaSet, 1, 15, new Tolerance(10));
            var scores = informedScorer.GetScores(AminoAcid.ProteinNTerm, sequence, AminoAcid.ProteinCTerm, composition, charge, scanNum);
            Console.WriteLine("Total Score = " + scores.Score);
            Console.WriteLine("#Fragments = " + scores.NumMatchedFrags);
        }
    }
}
