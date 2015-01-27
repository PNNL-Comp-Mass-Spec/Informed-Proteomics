using InformedProteomics.Backend.Data.Sequence;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestAminoAcidSet
    {
        public void TestParsingManyMods()
        {
            const string modFilePath = @"H:\Research\Lewy\MSPathFinder_Mods.txt";
            var aaSet = new AminoAcidSet(modFilePath);
            aaSet.Display();
        }
    }
}
