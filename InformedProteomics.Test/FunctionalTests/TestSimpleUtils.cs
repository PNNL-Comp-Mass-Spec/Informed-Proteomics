using System;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestSimpleUtils
    {
        [Test]
        public void TestStringShuffling()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string str = "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG";  // Histone H4
            var shuffled = SimpleStringProcessing.Shuffle(str);
            
            var strSorted = String.Concat(str.OrderBy(c => c));
            var shuffledSorted = String.Concat(shuffled.OrderBy(c => c));
            Assert.IsTrue(strSorted.Equals(shuffledSorted));
        }

        [Test]
        public void TestStringMutation()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string str = "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG";  // Histone H4
            const int numMutations = 3;
            var mutated = SimpleStringProcessing.Mutate(str, numMutations);

            Console.WriteLine(mutated);
            Assert.IsTrue(str.Length == mutated.Length);

            var numDiff = str.Where((t, i) => t != mutated[i]).Count();
            Console.WriteLine("Mutations: {0}", numDiff);

            //var strSorted = String.Concat(str.OrderBy(c => c));
            //var shuffledSorted = String.Concat(mutated.OrderBy(c => c));
            //Assert.IsTrue(strSorted.Equals(shuffledSorted));
        }

        [Test]
        public void TestCompositionOperations()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var comp = new Composition(1, 2, 3, 4, 5);
            var massComp = new CompositionWithDeltaMass(15.995);
            var sum = comp + massComp + comp + massComp;
            Console.WriteLine("{0}\t{1}\t{2}", sum, sum.Mass, sum.NominalMass);
            Console.WriteLine(sum.GetIsotopomerEnvelope().MostAbundantIsotopeIndex);
        }

        [Test]
        public void TestDblToString()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            TestValue(0, 0, "0");
            TestValue(0, 1, "0");
            TestValue(0, 2, "0");
            TestValue(0, 3, "0");

            TestValue(0.1, 0, "0");
            TestValue(0.1, 1, "0.1");
            TestValue(0.1, 2, "0.1");

            TestValue(0.1234, 0, "0");
            TestValue(0.1234, 1, "0.1");
            TestValue(0.1234, 2, "0.12");
            TestValue(0.1234, 4, "0.1234");
            TestValue(0.1234, 8, "0.1234");

            TestValue(0.987654321, 0, "1");
            TestValue(0.987654321, 1, "1.0");
            TestValue(0.987654321, 2, "0.99");
            TestValue(0.987654321, 4, "0.9877");
            TestValue(0.987654321, 8, "0.98765432");
            TestValue(0.987654321, 9, "0.987654321");
            TestValue(0.987654321, 12, "0.987654321");

            TestValue(-0.987654321, 0, "-1");
            TestValue(-0.987654321, 1, "-1.0");
            TestValue(-0.987654321, 2, "-0.99");
            TestValue(-0.987654321, 4, "-0.9877");
            TestValue(-0.987654321, 8, "-0.98765432");
            TestValue(-0.987654321, 9, "-0.987654321");
            TestValue(-0.987654321, 12, "-0.987654321");

            TestValue(0.00009876, 0, "0");
            TestValue(0.00009876, 1, "0.0");
            TestValue(0.00009876, 2, "0.0");
            TestValue(0.00009876, 3, "0.0");
            TestValue(0.00009876, 4, "0.0001");
            TestValue(0.00009876, 5, "0.0001");
            TestValue(0.00009876, 6, "0.000099");
            TestValue(0.00009876, 7, "0.0000988");
            TestValue(0.00009876, 8, "0.00009876");
            TestValue(0.00009876, 9, "0.00009876");

            TestValue(0.00004002, 0, "0");
            TestValue(0.00004002, 1, "0.0");
            TestValue(0.00004002, 2, "0.0");
            TestValue(0.00004002, 3, "0.0");
            TestValue(0.00004002, 4, "0.0");
            TestValue(0.00004002, 5, "0.00004");
            TestValue(0.00004002, 6, "0.00004");
            TestValue(0.00004002, 7, "0.00004");
            TestValue(0.00004002, 8, "0.00004002");

            TestValue(-0.00004002, 0, "0");
            TestValue(-0.00004002, 1, "0.0");
            TestValue(-0.00004002, 2, "0.0");
            TestValue(-0.00004002, 3, "0.0");
            TestValue(-0.00004002, 4, "0.0");
            TestValue(-0.00004002, 5, "-0.00004");
            TestValue(-0.00004002, 6, "-0.00004");
            TestValue(-0.00004002, 7, "-0.00004");
            TestValue(-0.00004002, 8, "-0.00004002");

        }

        private void TestValue(double value, byte digitsOfPrecision, string resultExpected)
        {
            var result = Misc.DblToString(value, digitsOfPrecision);
            Console.WriteLine(@"{0,12}, digits={1,2}: {2}", value, digitsOfPrecision, result);

            Assert.IsTrue(string.CompareOrdinal(result, resultExpected) == 0, "Result " + result + " did not match expected result (" + resultExpected + ")");
        }
    }
}
