using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    class TestReadingMsgfPlusResults
    {
        [Test]
        public void TestReadingTmtResultFile()
        {
            const string filePath = @"\\protoapps\UserData\Sangtae\TestData\MSGFPlusResultTMT10.tsv";
            var parser = new TsvFileParser(filePath);
            var pepStrs = parser.GetData("Peptide");
            var formulaStrs = parser.GetData("Formula");

            Assert.True(pepStrs.Count == formulaStrs.Count);

            var peptides = pepStrs.Select(Sequence.GetSequenceFromMsGfPlusPeptideStr).ToList();
            var formulae = formulaStrs.Select(Composition.Parse).ToList();

            Assert.True(peptides.Count == formulae.Count);

            for (var i = 0; i < peptides.Count; i++)
            {
                Assert.True((peptides[i].Composition + Composition.H2O).Equals(formulae[i]));
            }
        }
    }
}
