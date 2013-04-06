using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestUtils
    {
        [Test]
        public void TestCompositions()
        {
            var comp1 = Modification.Oxidation.Composition;
            var comp2 = new Composition(1, 2, 3, 4, 5,
                                        new[]
                                            {
                                                new Tuple<Atom, short>(Atom.Get("Au"), 2),
                                                new Tuple<Atom, short>(Atom.Get("P"), 2)
                                            });
            var comp3 = comp1 + comp2;

            Console.WriteLine(comp1);
            Console.WriteLine(comp2);
            Console.WriteLine(comp3);
        }

        [Test]
        public void TestModificationParams()
        {
            var modifications = new[] { Modification.Acetylation, Modification.Phosphorylation, Modification.Oxidation}; //, Modification.PyroGluQ };
            var modParams = new ModificationParams(modifications, 3);
            int numCombinations = modParams.GetNumModificationCombinations();
            for (int modCombIndex = 0; modCombIndex < numCombinations; modCombIndex++)
            {
                var modCombination = modParams.GetModificationCombination(modCombIndex);
                Console.WriteLine("{0}: {1} {2} {3}", modCombIndex, modCombination.ToString(), modCombination.Composition, modCombination.Composition.GetMass());
                
            }
        }

        [Test]
        public void TestGeneratingCombinations()
        {
            const int n = 5;
            const int k = 3;
            var combinations = SimpleMath.GetCombinationsWithRepetition(n, k);
            int count = 0;
            foreach (var combination in combinations)
            {
                ++count;
                Console.WriteLine((count-1)+": "+string.Join(",", combination));
            }
            Console.WriteLine("Count: " + count);
            Console.WriteLine("Count2: " + SimpleMath.NChooseK(n + k - 1, k));
        }

        [Test]
        public void TestCSharpSyntax()
        {
            const string path = @"C:\cygwin\home\kims336\Developments\InformedProteomics\InformedProteomics.Test\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            Console.WriteLine(Path.GetFileName(path));
            Console.WriteLine(Path.GetFileNameWithoutExtension(path));
            Console.WriteLine(Path.GetExtension(path));
            Console.WriteLine(Path.GetDirectoryName(path));

            var test = new Dictionary<int, int>();
            test[1] = 2;
            Console.WriteLine(test.Equals(null));
        }

        [Test]
        public void TestPeptide()
        {
            const string sequence = "PEPTIDE";
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var compositions = aaSet.GetCompositions(sequence);

            Composition[] compositionArr = compositions.ToArray();
            Assert.AreEqual(compositionArr.Count(), 1);
            Console.WriteLine(compositionArr[0]);
            Console.WriteLine(compositionArr[0].GetMass());
            Console.WriteLine(compositionArr[0].GetNominalMass());
            // 2nd isotope
            Console.WriteLine(compositionArr[0].GetIsotopeMass(0));
            Console.WriteLine(compositionArr[0].GetIsotopeMass(1));
            Console.WriteLine(compositionArr[0].GetIsotopeMass(2));
            Assert.AreEqual(compositionArr[0].ToString(), "C34H51N7O14S0");

            foreach (var e in compositionArr[0].GetApproximatedIsotopomerEnvelop())
                Console.WriteLine(e);
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
