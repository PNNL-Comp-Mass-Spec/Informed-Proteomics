using System;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestUtils
    {
        [Test]
        public void TestCSharpSyntax()
        {
            const string path = @"C:\cygwin\home\kims336\Developments\InformedProteomics\InformedProteomics.Test\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            Console.WriteLine(Path.GetFileName(path));
            Console.WriteLine(Path.GetFileNameWithoutExtension(path));
            Console.WriteLine(Path.GetExtension(path));
            Console.WriteLine(Path.GetDirectoryName(path));
        }

        [Test]
        public void TestCompositions()
        {
            const string sequence = "PEPTIDE";
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var compositions = aaSet.GetCompositions(sequence);

            Composition[] compositionArr = compositions.ToArray();
            Assert.AreEqual(compositionArr.Count(), 1);
            Console.WriteLine(compositionArr[0]);
            Console.WriteLine(compositionArr[0].GetMass());
            Console.WriteLine(compositionArr[0].GetNominalMass());
            Assert.AreEqual(compositionArr[0].ToString(), "C34H51N7O14S0");
        }

        [Test]
        public void TestIonTypeGeneration()
        {
            var ionTypeFactory = new IonTypeFactory();
            int index = 0;
            foreach (var ionType in ionTypeFactory.GetAllKnownIonTypes())
            {
                Console.WriteLine(++index + ": " + ionType);
            }
            var yIon = ionTypeFactory.GetIonType("y2-H2O");
            Console.WriteLine(yIon.GetMz(0));
        }
	}
}
