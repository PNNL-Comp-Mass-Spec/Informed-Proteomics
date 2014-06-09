using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestFdrCalculation
    {
        public void TestIcTopDown()
        {
            const string targetResultPath = @"D:\Research\Data\Vlad\Ic\NewMod_M1\Alz_RA_C1_CID_11012013_SW_03Nov2013_IcTarget.tsv";
            const string decoyResultPath = @"D:\Research\Data\Vlad\Ic\NewMod_M1\Alz_RA_C1_CID_11012013_SW_03Nov2013_IcDecoy.tsv";
            const string tdaResultPath = @"D:\Research\Data\Vlad\Ic\NewMod_M1\Alz_RA_C1_CID_11012013_SW_03Nov2013_IcTda.tsv";
            //const string targetResultPath = @"C:\cygwin\home\kims336\Data\TopDown\raw\SBEP_STM_001_02272012_Aragon.icresult";
            //const string decoyResultPath = @"C:\cygwin\home\kims336\Data\TopDown\raw\SBEP_STM_001_02272012_Aragon.decoy.icresult";
            var fdrCalculator = new FdrCalculator(targetResultPath, decoyResultPath);
            fdrCalculator.WriteTo(tdaResultPath);
            Console.WriteLine("Done");
        }

        public void MergeTargetDecoyFiles()
        {
            const string dir = @"C:\cygwin\home\kims336\Data\TopDown\raw\Cache";
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
