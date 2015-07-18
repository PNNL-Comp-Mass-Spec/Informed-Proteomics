using System;
using System.IO;
using System.Linq;
using System.Reflection;
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
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);
           
            const string filePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MSGFPlusResultTMT10.tsv";
            if (!File.Exists(filePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, filePath);
                return;
            }

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
