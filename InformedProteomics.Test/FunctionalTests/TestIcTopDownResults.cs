using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestIcTopDownResults
    {
        [Test]
        public void ValidateIcResultsWithModifications()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string resultFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1_Rescored.tsv";

            if (!File.Exists(resultFilePath))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, resultFilePath);
                return;
            }

            var parser = new TsvFileParser(resultFilePath);
            var sequences = parser.GetData("Sequence");
            var modifications = parser.GetData("Modifications");
            var compositions = parser.GetData("Composition").Select(Composition.Parse).ToArray();
            var scanNums = parser.GetData("ScanNum").Select(s => Convert.ToInt32(s)).ToArray();
            var aaSet = new AminoAcidSet();
            for (var i = 0; i < parser.NumData; i++)
            {
                var sequenceComp = aaSet.GetComposition(sequences[i]) + Composition.H2O;

                var modComposition = Composition.Zero;
                var modsStr = modifications[i].Substring(1, modifications[i].Length - 2);
                var mods = modsStr.Split(',');
                foreach(var modStr in mods)
                {
                    if (modStr.Length == 0) continue;
                    var modName = modStr.Split()[0];
                    var mod = Modification.Get(modName);
                    modComposition += mod.Composition;
                }

                var compFromSeqAndMods = sequenceComp + modComposition;
                Assert.True(compFromSeqAndMods.Equals(compositions[i]));
            }
        }
    }
}
