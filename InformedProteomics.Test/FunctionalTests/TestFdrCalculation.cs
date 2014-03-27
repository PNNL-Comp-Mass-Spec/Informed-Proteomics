using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestFdrCalculation
    {
        public void TestIcTopDown()
        {
            const string targetResultPath = @"\\protoapps\UserData\Sangtae\TestData\SBEP_STM_001_02272012_Aragon.icresult";
            const string decoyResultPath = @"\\protoapps\UserData\Sangtae\TestData\SBEP_STM_001_02272012_Aragon.decoy.icresult";
            var fdrCalculator = new FdrCalculator(targetResultPath, decoyResultPath);
            fdrCalculator.WriteTo(Path.ChangeExtension(targetResultPath, ".tsv"));
            Console.WriteLine("Done");
        }

        public void MergeTargetDecoyFiles()
        {
            const string dir = @"C:\cygwin\home\kims336\Data\TopDown\raw\CorrMatches";
            var rawFileNames = new HashSet<string>();
            foreach (var f in Directory.GetFiles(dir, "*.icresult"))
            {
                rawFileNames.Add(f.Substring(0, f.IndexOf('.')));
            }

            foreach (var rawFileName in rawFileNames)
            {
                var targetResultFilePath = rawFileName + ".icresult";
                var decoyResultFilePath = rawFileName + ".decoy.icresult";
                var mergedResultFilePath = rawFileName + ".tsv";

                Console.Write("Creating {0}...", mergedResultFilePath);
                var fdrCalculator = new FdrCalculator(targetResultFilePath, decoyResultFilePath);
                fdrCalculator.WriteTo(mergedResultFilePath);
                Console.WriteLine("Done");
            }
        }
    }
}
