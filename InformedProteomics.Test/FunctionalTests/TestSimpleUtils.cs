using System;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;
using PNNLOmics.Utilities;

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

            TestValue(1, 1, "1.0");
            TestValue(1, 3, "1.0");
            TestValue(5, 3, "5.0");

            TestValue(10, 0, "10");
            TestValue(10, 1, "10");
            TestValue(10, 2, "10");
            TestValue(10, 3, "10");

            TestValue(10.123, 0, "10");
            TestValue(10.123, 1, "10.1");
            TestValue(10.123, 2, "10.12");
            TestValue(10.123, 3, "10.123");

            TestValue(50, 0, "50");
            TestValue(50, 2, "50");
            TestValue(50, 4, "50");
            
            TestValue(50.653, 0, "51");
            TestValue(50.653, 1, "50.7");
            TestValue(50.653, 2, "50.65");
            TestValue(50.653, 3, "50.653");
            TestValue(50.653, 4, "50.653");

            TestValue(50.753, 0, "51");
            TestValue(50.753, 1, "50.8");
            TestValue(50.753, 2, "50.75");
            TestValue(50.753, 3, "50.753");
            TestValue(50.753, 4, "50.753");

            TestValue(110, 0, "110");
            TestValue(110, 1, "110");
            TestValue(110, 2, "110");

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
            var result = StringUtilities.DblToString(value, digitsOfPrecision);
            Console.WriteLine(@"{0,12}, digits={1,2}: {2}", value, digitsOfPrecision, result);

            Assert.IsTrue(string.CompareOrdinal(result, resultExpected) == 0, "Result " + result + " did not match expected result (" + resultExpected + ")");
        }
    }
}
