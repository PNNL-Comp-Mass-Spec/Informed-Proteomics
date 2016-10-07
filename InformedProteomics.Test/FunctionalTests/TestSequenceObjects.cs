using System;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestSequenceObjects
    {
        [Test]
        public void TestCompositions()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var comp1 = Modification.Oxidation.Composition;
            var comp2 = new Composition(1, 2, 3, 4, 5, 6,
                                        new[]
                                            {
                                                new Tuple<Atom, short>(Atom.Get("Au"), 2),
                                                new Tuple<Atom, short>(Atom.Get("13C"), 3)
                                            });
            var comp3 = comp1 + comp2;

            Assert.AreEqual(comp1.ToString(), "C(0) H(0) N(0) O(1) S(0)");
            Assert.AreEqual(comp1.ToPlainString(), "O1");
            Assert.AreEqual(comp2.ToString(), "C(1) H(2) N(3) O(4) S(5) P(6) Au(2) 13C(3)");
            Assert.AreEqual(comp2.ToPlainString(), "C1H2N3O4S5P6Au213C3");
            Assert.AreEqual(comp3.ToString(), "C(1) H(2) N(3) O(5) S(5) P(6) Au(2) 13C(3)");
            var comp4 = new Composition(1, 2, 3, 4, 5, 6,
                                        new[]
                                            {
                                                new Tuple<Atom, short>(Atom.Get("Au"), 2),
                                                new Tuple<Atom, short>(Atom.Get("13C"), 3)
                                            });

            // Testing GetHashCode() and Equals()
            Assert.IsTrue(comp2.Equals(comp4));
            Assert.IsTrue(comp2.GetHashCode() == comp4.GetHashCode());
            Assert.IsTrue(comp2.Equals(comp2 + Composition.Zero));
        }
    }
}
