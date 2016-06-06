using System.Reflection;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.FunctionalTests.FunctionalTests
{
    [TestFixture]
    public class TestTopDownConsole
    {
        [Test]
        public void TestManyModSearch()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string rawFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\Lewy_ManyMods\Lewy_intact_01.raw";
            const string featureFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\Lewy_ManyMods\Lewy_intact_01.ms1ft";
            const string databaseFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\Lewy_ManyMods\ID_004858_0EE8CF61.fasta";
            const string modFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\Lewy_ManyMods\Mods.txt";
            const string testOutputPath =
                @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\TopDown\Lewy_ManyMods\TestOutput";

            var args = new[]
            {
                "-s", rawFilePath,
                "-d", databaseFilePath,
                "-o", testOutputPath,
                "-m", "1",
                "-mod", modFilePath,
                "-tda", "1",
                "-feature", featureFilePath
            };

            MSPathFinderT.Program.Main(args);
        }
    }
}
