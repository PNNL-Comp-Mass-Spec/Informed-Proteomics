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
    public class TestOffsetFrequencyFunction
    {
        [Test]
        public void PrintOffsetFrequency()
        {
            var spectrumMatches = InitTest();

            using (var debugFile = new StreamWriter(OutputFileName))
            {
                foreach (var spectrumMatch in spectrumMatches)
                {

                    var offsetFrequencyTable = new PrecursorOffsetFrequencyTable(100, spectrumMatch.PrecursorCharge,
                        1.005/spectrumMatch.PrecursorCharge);

                    var ionType = _ionTypes[spectrumMatch.PrecursorCharge - 1];
                    var ion = ionType.GetIon(spectrumMatch.PeptideComposition);
                    var mz = ion.GetMonoIsotopicMz();

                    offsetFrequencyTable.AddMatches(new List<SpectrumMatch> {spectrumMatch});
                    var offsetFrequencies = offsetFrequencyTable.OffsetFrequencies;

                    debugFile.WriteLine("ScanNum\tM/Z\tPrecursor Charge\tIon Type\tPeptide");
                    debugFile.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}",
                        spectrumMatch.ScanNum, mz, spectrumMatch.PrecursorCharge, ionType.Name, spectrumMatch.Peptide);
                    debugFile.WriteLine("Offset\tM/Z");
                    foreach (var offsetFrequency in offsetFrequencies)
                    {
                        if (offsetFrequency.Found > 0)
                            debugFile.WriteLine("{0}\t{1}", offsetFrequency.Offset, offsetFrequency.Offset+mz);
                    }
                }
            }
        }

        // Configuration
//        private const int ScanNum = 30623;
        private const double NoiseFiltration = 0;
        private const ActivationMethod Act = ActivationMethod.HCD;
        private const int MaxCharge = 10;
        private readonly BaseIonType[] _baseIons = { BaseIonType.Y };
        private readonly NeutralLoss[] _neutralLosses = { NeutralLoss.NoLoss };
        private List<IonType> _ionTypes;

        private const string TsvFile = @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\tsv\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.tsv";
        private const string RawFile = @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\raw\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
        private const string OutputFileName = @"C:\Users\wilk011\Documents\DataFiles\TestFolder\ouput_HCD_QCShew_Precursor.txt";
//        private const string DebugFileName = @"C:\Users\wilk011\Documents\DataFiles\TestFolder\debug_HCD_QCShew_Precursor.txt";
        private IEnumerable<SpectrumMatch> InitTest()
        {
            var ionTypeFactory = new IonTypeFactory(_baseIons, _neutralLosses, MaxCharge);

            _ionTypes = ionTypeFactory.GetAllKnownIonTypes().ToList();

            var lcms = LcMsRun.GetLcMsRun(RawFile, MassSpecDataType.XCaliburRun, NoiseFiltration, NoiseFiltration);

            var spectrumMatches = (new SpectrumMatchList(lcms, new TsvFileParser(TsvFile), Act)).Matches;

            return spectrumMatches;
        }
    }
}
