using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using NUnit.Framework;
using SelectivityScore;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestSelectivityScore
    {
        [Test]
        public void TestSelectivityScoreUsingNominalMasses()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var aaSet = AminoAcidSet.GetStandardAminoAcidSetWithCarboamidomethylCys();

            var aaList = new List<AminoAcid>();
            foreach (var residue in AminoAcid.StandardAminoAcidCharacters)
            {
                var aa = aaSet.GetAminoAcid(residue);
                aaList.Add(aa);
                if(residue == 'M') aaList.Add(new ModifiedAminoAcid(aa, Modification.Oxidation));
            }

            var nominalMassArr = aaList.Select(aa => aa.GetNominalMass()).ToArray();
            var probArr = Enumerable.Repeat(0.05, nominalMassArr.Length).ToArray();
            
            var selScoreCalculator = new SelectivityScoreCalculatorUsingNominalMasses(
                nominalMassArr, probArr);
            var hist = selScoreCalculator.GetPValues(848.41, new[] {1259.5989, 1130.5563, 1059.5192, 958.4715});
            for (var t = hist.Length - 1; t >= 0; t--)
            {
                Console.WriteLine("{0}\t{1}", t, hist[t]);
            }
        }

        [Test]
        public void TestBinSize()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const double mass = 907.4763;
            const double minMass = mass - 25.0;
            const double maxMass = mass + 25.0;

            var minBinNum = Constants.GetBinNumHighPrecision(minMass);
            var maxBinNum = Constants.GetBinNumHighPrecision(maxMass);

            Console.Write("{0}\t{1}\t{2}", minBinNum, maxBinNum, (maxBinNum-minBinNum+1));
        }
    }
}
