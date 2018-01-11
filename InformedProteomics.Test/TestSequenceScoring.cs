using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;
using InformedProteomics.TopDown;
using InformedProteomics.TopDown.Scoring;

namespace InformedProteomics.Test
{
    [TestFixture]
    class TestSequenceScoring
    {
        [Test]
        public void TestScoring()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var rawFile = @"\\protoapps\UserData\Jungkap\Joshua\testData\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf" ;
            var resultFile = @"\\protoapps\UserData\Jungkap\Joshua\IdResult\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv";

            if (!File.Exists(rawFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, rawFile);
                return;
            }

            if (!File.Exists(resultFile))
            {
                Console.WriteLine(@"Warning: Skipping test {0} since file not found: {1}", methodName, resultFile);
                return;
            }

            var tsvParser = new TsvFileParser(resultFile);
            var tsvData = tsvParser.GetAllData();
            var ms2ScanNumbers = tsvData["Scan"];

            var run = PbfLcMsRun.GetLcMsRun(rawFile, 0, 0);

            for (int i = 0; i < 1; i++)
            {
                var scanNum = Int32.Parse(ms2ScanNumbers[i]);
                var spectrum = run.GetSpectrum(scanNum) as ProductSpectrum;
                int tsvIndex = ms2ScanNumbers.FindIndex(x => Int32.Parse(x) == scanNum);

                var seqStr = tsvData["Sequence"].ElementAt(tsvIndex).Trim();
                var seqMod = tsvData["Modifications"].ElementAt(tsvIndex).Trim();
                var aaSet = new AminoAcidSet();
                var sequence = Sequence.CreateSequence(seqStr, seqMod, aaSet);
                Console.WriteLine(sequence.Count);
                var score = GetScoreTest(sequence, spectrum);
                Console.WriteLine(scanNum + ":" + score);
            }
        }

        public double GetScoreTest(Sequence sequence, ProductSpectrum spectrum)
        {
            var score = 0d;
            var tol = new Tolerance(10);
            var matchCounter = new CorrMatchedPeakCounter(spectrum,tol,1,20);
            var prefixCompArr = sequence.GetPrefixCompositions().ToArray();
            foreach (var c in prefixCompArr)
            {
                if(c.Equals(Composition.Zero)) Console.WriteLine("found zero");
            }
            var suffixCompArr = sequence.GetSuffixCompositions().ToArray();
            for (int i = 0; i < prefixCompArr.Length; i++)
            {
                score += matchCounter.GetFragmentScore(prefixCompArr[i], suffixCompArr[i]);
            }
            return score;
        }
    }
}
