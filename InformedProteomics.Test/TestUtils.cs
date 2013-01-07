using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    class TestUtils
    {
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
            Assert.AreEqual(compositionArr[0].ToString(), "C34H51N7O14S0");
        }
    }
}
