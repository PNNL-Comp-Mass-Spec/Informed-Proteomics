using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.LikelihoodScoring;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestPrecursorFilter
    {
        [Test]
        public void FilterPrecursorPeaks()
        {
            const double probabilityThreshold = 0.15;
            const double searchWidth = 100;
            const int retentionCount = 6;
            const string rawFile = @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QE_EDRN\raw\342862_EDRN_Serum_02_DDA_1_12Nov13_Samwise_13-07-28.raw";
            const string tsvFile = @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QE_EDRN\tsv\342862_EDRN_Serum_02_DDA_1_12Nov13_Samwise_13-07-28.tsv";
            const string outFileName = @"C:\Users\wilk011\Documents\DataFiles\TestFolder\HCD_EQ_EDRN_FilterTest.txt";

            var precOff = new Dictionary<int, PrecursorOffsets>();

            var charge1 = new PrecursorOffsets(1, probabilityThreshold);
            charge1.AddOffsetsFromFile(@"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\HCD_QCShew_PrecursorOffsetFrequencies_Charge1.txt");
            precOff.Add(1, charge1);
            var charge2 = new PrecursorOffsets(2, probabilityThreshold);
            charge2.AddOffsetsFromFile(@"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\HCD_QCShew_PrecursorOffsetFrequencies_Charge2.txt");
            precOff.Add(2, charge2);
            var charge3 = new PrecursorOffsets(3, probabilityThreshold);
            charge3.AddOffsetsFromFile(@"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\HCD_QCShew_PrecursorOffsetFrequencies_Charge3.txt");
            precOff.Add(3, charge3);
            var charge4 = new PrecursorOffsets(4, probabilityThreshold);
            charge4.AddOffsetsFromFile(@"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\HCD_QCShew_PrecursorOffsetFrequencies_Charge4.txt");
            precOff.Add(4, charge4);

            var filter = new PrecursorFilter(precOff, new Tolerance(10, ToleranceUnit.Ppm), searchWidth, retentionCount);

            var lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0);
            var spectrumMatches = (new SpectrumMatchList(lcms, new TsvFileParser(tsvFile), ActivationMethod.HCD)).Matches;

            var specMatch = spectrumMatches.First();
            var filteredSpecMatch = (filter.FilterMatches(new List<SpectrumMatch> {specMatch})).First();

            var unfilteredPeaks = specMatch.Spectrum.Peaks;
            var filteredPeaks = filteredSpecMatch.Spectrum.Peaks;

            using (var outFile = new StreamWriter(outFileName))
            {
                outFile.WriteLine("M/Z\tFiltered?");
                for (int i = 0; i < unfilteredPeaks.Length; i++)
                {
                    outFile.Write("{0}\t", Math.Round(unfilteredPeaks[i].Mz, 2));
                    if (i < filteredPeaks.Length)
                    {
                        outFile.Write(filteredPeaks[i].Mz.Equals(unfilteredPeaks[i].Mz) ? "No" : "Yes");
                    }
                    else
                    {
                        outFile.Write("Yes");
                    }
                    outFile.WriteLine();
                }
            }
        }
    }
}
