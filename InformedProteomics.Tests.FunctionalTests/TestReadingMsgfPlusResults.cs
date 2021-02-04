using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.FunctionalTests
{
    [TestFixture]
    internal class TestReadingMsgfPlusResults
    {
        [Test]
        [Category("Local_Testing")]
        public void TestReadingTmtResultFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var filePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, "MSGFPlusResultTMT10.tsv");
            if (!File.Exists(filePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, filePath);
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
