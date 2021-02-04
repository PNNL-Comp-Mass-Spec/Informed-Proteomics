using System;
using System.Reflection;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.FunctionalTests
{
    [TestFixture]
    public class TestAminoAcidSet
    {
        [Test]
        [TestCase(@"TEST_FOLDER\TopDown\Lewy_ManyMods\Lewy_DB_Mods.txt")]
        [TestCase(@"TEST_FOLDER\TopDown\Lewy_ManyMods\Lewy_DB_Mods2.txt")]
        public void TestParsingManyMods(string modDefsFile)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;

            var modFile = Utils.GetTestFile(methodName, modDefsFile.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));

            if (!modFile.Exists)
            {
                Assert.Ignore("Ignoring test TestParsingManyMods since file not found: " + modFile.FullName);
            }

            var aaSet = new AminoAcidSet(modFile.FullName);
            //aaSet.Display();

            //SequenceLocation.ProteinNTerm
            var residue = AminoAcid.ProteinNTerm.Residue;
            var location = SequenceLocation.ProteinNTerm;
            var aa = aaSet.GetAminoAcid(residue, location);
            Console.Write("{0}\t{1}\t{2}", residue, aa.Mass, aa.Composition);
            foreach (var modIndex in aaSet.GetModificationIndices(residue, location))
            {
                var modification = aaSet.GetModificationParams().GetModification(modIndex);
                Console.WriteLine(modification.Mass);
                //Console.Write("\t" + _modificationParams.GetModification(modIndex));
            }
            Console.WriteLine();
            residue = AminoAcid.ProteinCTerm.Residue;
            location = SequenceLocation.ProteinCTerm;
            aa = aaSet.GetAminoAcid(residue, location);
            Console.Write("{0}\t{1}\t{2}", residue, aa.Mass, aa.Composition);
            foreach (var modIndex in aaSet.GetModificationIndices(residue, location))
            {
                var modification = aaSet.GetModificationParams().GetModification(modIndex);
                Console.WriteLine(modification.Mass);
                //Console.Write("\t" + _modificationParams.GetModification(modIndex));
            }

            //foreach (var aa in AminoAcid.StandardAminoAcidArr)
            //{
            /*
            var keys = _locationSpecificResidueMap[location].Keys.ToArray();
            Array.Sort(keys);
            foreach (var residue in keys)
            {
                var aa = GetAminoAcid(residue, location);
                Console.Write("{0}\t{1}\t{2}", residue, aa.Mass, aa.Composition);
                foreach (var modIndex in GetModificationIndices(residue, location))
                {
                    Console.Write("\t" + _modificationParams.GetModification(modIndex));
                }
                Console.WriteLine();
            }*/
            //}
        }

        [Ignore("File Missing, test obsolete, or long test")]
        [Category("Local_Testing")]
        public void TestParsingGlycoMods()
        {
            const string modFilePath = @"C:\cygwin\home\kims336\Data\Debug\MSPathFinder_Mods.txt";
            var aaSet = new AminoAcidSet(modFilePath);
            aaSet.Display();
        }
    }
}
