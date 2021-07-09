using System;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Tests.Base;
using NUnit.Framework;
using PRISM;

namespace InformedProteomics.Tests.UnitTests
{
    [TestFixture]
    public class TestSimpleUtils
    {
        [Test]
        public void TestStringShuffling()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string str = "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG";  // Histone H4
            var shuffled = SimpleStringProcessing.Shuffle(str);

            var strSorted = string.Concat(str.OrderBy(c => c));
            var shuffledSorted = string.Concat(shuffled.OrderBy(c => c));
            Assert.IsTrue(strSorted.Equals(shuffledSorted));
        }

        [Test]
        public void TestStringMutation()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string str = "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG";  // Histone H4
            const int numMutations = 3;
            var mutated = SimpleStringProcessing.Mutate(str, numMutations);

            Console.WriteLine(mutated);
            Assert.IsTrue(str.Length == mutated.Length);

            var numDiff = str.Where((t, i) => t != mutated[i]).Count();
            Console.WriteLine("Mutations: {0}", numDiff);

            //var strSorted = string.Concat(str.OrderBy(c => c));
            //var shuffledSorted = string.Concat(mutated.OrderBy(c => c));
            //Assert.IsTrue(strSorted.Equals(shuffledSorted));
        }

        [Test]
        public void TestCompositionOperations()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

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
            Utils.ShowStarting(methodName);

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

            TestValue(54.753, 0, "55");
            TestValue(54.753, 1, "54.8");
            TestValue(54.753, 2, "54.75");
            TestValue(54.753, 3, "54.753");
            TestValue(54.753, 4, "54.753");

            TestValue(110, 0, "110");
            TestValue(110, 1, "110");
            TestValue(110, 2, "110");

            TestValue(9.99999, 6, "9.99999");
            TestValue(9.99999, 5, "9.99999");
            TestValue(9.99999, 4, "10.0");
            TestValue(9.99999, 2, "10.0");
            TestValue(9.99999, 1, "10.0");
            TestValue(9.99999, 0, "10");

            TestValue(9.98765, 6, "9.98765");
            TestValue(9.98765, 5, "9.98765");
            TestValue(9.98765, 4, "9.9877");
            TestValue(9.98765, 3, "9.988");
            TestValue(9.98765, 2, "9.99");
            TestValue(9.98765, 1, "10.0");

            TestValue(5.12345, 5, "5.12345", false);
            TestValue(50.12345, 5, "50.12345", false);
            TestValue(500.12345, 5, "500.12345", false);
            TestValue(5000.12345, 5, "5000.12345", false);
            TestValue(50000.12345, 5, "50000.12345", false);
            TestValue(500000.12345, 5, "500000.12345", false);

            TestValue(5.12345, 5, "5.12345", true);
            TestValue(50.12345, 5, "50.1235", true);
            TestValue(500.12345, 5, "500.123", true);
            TestValue(5000.12345, 5, "5000.12", true);
            TestValue(50000.12345, 5, "50000.1", true);
            TestValue(500000.12345, 5, "500000", true);

            TestValue(9.98765, 3, "9.988", true);
            TestValue(99.98765, 3, "99.99", true);
            TestValue(998.98765, 3, "999.0", true);
            TestValue(9987.98765, 3, "9988", true);
            TestValue(99876.98765, 3, "99877", true);
            TestValue(998765.98765, 3, "998766", true);

            TestValue(678.8032741, 9, "678.8032741", true);
            TestValue(678.8032741, 8, "678.803274", true);
            TestValue(678.8032741, 7, "678.80327", true);

            TestValue(1017.9488123, 9, "1017.948812", true);
            TestValue(1017.9488123, 8, "1017.94881", true);
            TestValue(1017.9488123, 7, "1017.9488", true);

            TestValue(10641.5439241, 9, "10641.54392", true);
            TestValue(10641.5439241, 8, "10641.5439", true);
            TestValue(10641.5439241, 7, "10641.544", true);

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

            TestValue(0.00009876, 0, "0", 1E-6);
            TestValue(0.00009876, 1, "0.0", 1E-6);
            TestValue(0.00009876, 2, "0.0", 1E-6);
            TestValue(0.00009876, 3, "0.0", 1E-6);
            TestValue(0.00009876, 4, "0.0001", 1E-6);
            TestValue(0.00009876, 5, "0.0001", 1E-6);
            TestValue(0.00009876, 6, "0.000099", 1E-6);
            TestValue(0.00009876, 7, "0.0000988", 1E-6);
            TestValue(0.00009876, 8, "0.00009876", 1E-6);
            TestValue(0.00009876, 9, "0.00009876", 1E-6);

            TestValue(0.00009876, 0, "0");
            TestValue(0.00009876, 1, "9.9E-05");
            TestValue(0.00009876, 2, "9.88E-05");
            TestValue(0.00009876, 3, "9.876E-05");
            TestValue(0.00009876, 4, "9.876E-05");
            TestValue(0.00009876, 5, "9.876E-05");

            TestValue(0.00004002, 0, "0", 1E-6);
            TestValue(0.00004002, 1, "0.0", 1E-6);
            TestValue(0.00004002, 2, "0.0", 1E-6);
            TestValue(0.00004002, 3, "0.0", 1E-6);
            TestValue(0.00004002, 4, "0.0", 1E-6);
            TestValue(0.00004002, 5, "0.00004", 1E-6);
            TestValue(0.00004002, 6, "0.00004", 1E-6);
            TestValue(0.00004002, 7, "0.00004", 1E-6);
            TestValue(0.00004002, 8, "0.00004002", 1E-6);

            TestValue(0.00004002, 0, "0");
            TestValue(0.00004002, 1, "4.0E-05");
            TestValue(0.00004002, 2, "4.0E-05");
            TestValue(0.00004002, 3, "4.002E-05");
            TestValue(0.00004002, 4, "4.002E-05");

            TestValue(-0.00004002, 0, "0", 1E-6);
            TestValue(-0.00004002, 1, "0.0", 1E-6);
            TestValue(-0.00004002, 2, "0.0", 1E-6);
            TestValue(-0.00004002, 3, "0.0", 1E-6);
            TestValue(-0.00004002, 4, "0.0", 1E-6);
            TestValue(-0.00004002, 5, "-0.00004", 1E-6);
            TestValue(-0.00004002, 6, "-0.00004", 1E-6);
            TestValue(-0.00004002, 7, "-0.00004", 1E-6);
            TestValue(-0.00004002, 8, "-0.00004002", 1E-6);

            TestValue(-0.00004002, 0, "0");
            TestValue(-0.00004002, 1, "-4.0E-05");
            TestValue(-0.00004002, 2, "-4.0E-05");
            TestValue(-0.00004002, 3, "-4.002E-05");
            TestValue(-0.00004002, 4, "-4.002E-05");

            TestValue(4.94065645E-324, 6, "4.940656E-324");
            TestValue(4.94065645E-150, 6, "4.940656E-150");
            TestValue(4.94065645E-101, 6, "4.940656E-101");
            TestValue(4.02735019E-10, 6, "4.02735E-10");

            TestValue(4.0273501E-5, 6, "4.02735E-05");
            TestValue(4.0273501E-4, 6, "4.02735E-04");
            TestValue(4.0273501E-3, 6, "0.004027");
            TestValue(4.0273501E-2, 6, "0.040274");
            TestValue(4.0273501E-1, 6, "0.402735");
            TestValue(0.0134886, 6, "0.013489");

            TestValue(4.0273501E-10, 6, "0.0", true, 1E-15);
            TestValue(0.0134886, 6, "0.013489", true);
            TestValue(7063.79431, 6, "7063.794", true);
            TestValue(6496286.95, 6, "6496287", true);
        }

        private void TestValue(
            double value,
            byte digitsAfterDecimal,
            string resultExpected,
            double thresholdScientific)
        {
            TestValue(value, digitsAfterDecimal, resultExpected, limitDecimalsForLargeValues: false, thresholdScientific: thresholdScientific);
        }

        private void TestValue(
            double value,
            byte digitsAfterDecimal,
            string resultExpected,
            bool limitDecimalsForLargeValues = false,
            double thresholdScientific = 0.001)
        {
            var result = StringUtilities.DblToString(value, digitsAfterDecimal, limitDecimalsForLargeValues, thresholdScientific);
            Console.Write("{0,20}, digits={1,2}: {2,-8}", value, digitsAfterDecimal, result);

            if (limitDecimalsForLargeValues)
            {
                Console.WriteLine(" (limitDecimals=True)");
            }
            else
            {
                Console.WriteLine();
            }

            var expectedResultFound = string.CompareOrdinal(result, resultExpected) == 0;

            if (expectedResultFound)
            {
                return;
            }

            var errMsg = "Result " + result + " did not match expected result (" + resultExpected + ")";

            Console.WriteLine(errMsg);
            Assert.IsTrue(expectedResultFound, errMsg);
        }
    }
}
