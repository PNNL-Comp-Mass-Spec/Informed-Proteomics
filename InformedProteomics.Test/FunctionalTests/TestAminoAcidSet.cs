using InformedProteomics.Backend.Data.Sequence;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestAminoAcidSet
    {
        [Ignore]
        public void TestParsingManyMods()
        {
            const string modFilePath = @"H:\Research\Lewy\MSPathFinder_Mods.txt";
            var aaSet = new AminoAcidSet(modFilePath);
            aaSet.Display();
        }

        [Ignore]
        public void TestParsingGlycoMods()
        {
            const string modFilePath = @"C:\cygwin\home\kims336\Data\Debug\MSPathFinder_Mods.txt";
            var aaSet = new AminoAcidSet(modFilePath);
            aaSet.Display();
        }
    }
}
