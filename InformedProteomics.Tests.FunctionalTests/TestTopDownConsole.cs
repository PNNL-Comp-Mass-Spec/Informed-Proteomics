using System.IO;
using System.Reflection;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.FunctionalTests
{
    [TestFixture]
    public class TestTopDownConsole
    {
        [Test]
        [Category("PNL_Domain")]
        [Ignore("Long running")]
        public void TestManyModSearch()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var rawFilePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\Lewy_ManyMods\Lewy_intact_01.raw");
            var featureFilePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\Lewy_ManyMods\Lewy_intact_01.ms1ft");
            var databaseFilePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\Lewy_ManyMods\ID_004858_0EE8CF61.fasta");
            var modFilePath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\Lewy_ManyMods\Lewy_DB_Mods.txt");
            var testOutputPath = Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, @"TopDown\Lewy_ManyMods\TestOutput");

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
