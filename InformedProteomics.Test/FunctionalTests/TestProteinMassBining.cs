using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Scoring.GeneratingFunction;
using InformedProteomics.Scoring.TopDown;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    class TestProteinMassBining
    {
        [Test]
        public void TestProteinMassComparerWithBinning()
        {
            var comparer2 = new FilteredProteinMassBinning(new AminoAcidSet(),  50001);

            for (var i = 9999d; i < 10010; i++)
            {
                Console.WriteLine("{0}, {1}",i, comparer2.GetBinNumber(i));
            }

            //var comparer = new ProteinMassBinning(50, 50001, true);
            /*
            Console.WriteLine(Constants.GetBinNumHighPrecision(50000));
            Console.WriteLine(comparer.NumberOfBins);
            Console.WriteLine(comparer2.NumberOfBins);

            var rnd = new Random();

            var mass = 0d;
            for (var i = 0; i < 450; i ++)
            {
                if (i > 0)
                {
                    var j = rnd.Next(aaSet.Length);
                    mass += aaSet[j].Mass;
                }
                if (mass > comparer.MaxMass) break;

                var binNum = Constants.GetBinNumHighPrecision(mass);
                var binNum1 = comparer.GetBinNumber(mass);
                var binNum2 = comparer2.GetBinNumber(mass);

                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", i, mass, binNum, binNum1, binNum2);
            }*/
        }

        [Test]
        public void CompareWithLinearScaling()
        {
            var comparer = new ProteinMassBinning(50, 50001, true);
            var minMass = 0;
            var maxMass = 50000;
            var hBins = new List<int>();

            for (var m = minMass; m <= maxMass; m += 100)
            {
                var binNum1 = Constants.GetBinNumHighPrecision(m);
                var binNum2 = comparer.GetBinNumber(m);

                Console.WriteLine("{0}\t{1}\t{2}", m, binNum1, binNum2);
            }
        }
    }
}
