using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
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
            const string str = "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG";  // Histone H4
            var shuffled = SimpleStringProcessing.Shuffle(str);
            
            var strSorted = String.Concat(str.OrderBy(c => c));
            var shuffledSorted = String.Concat(shuffled.OrderBy(c => c));
            Assert.IsTrue(strSorted.Equals(shuffledSorted));
        }

        [Test]
        public void TestStringMutation()
        {
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

    }
}
