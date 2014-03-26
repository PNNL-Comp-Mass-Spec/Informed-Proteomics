using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
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
        public List<IonProbability> ComputeOffsetFrequencies(SpectrumMatch spectrumMatch)
        {
            var debugFile = File.AppendText(DebugFileName);

            var peptide = spectrumMatch.Peptide;
            var spectrum = spectrumMatch.Spectrum;
            var prefixes = spectrumMatch.Prefixes;
            var suffixes = spectrumMatch.Suffixes;
            suffixes.Reverse();

            var prefixIonTypes = _ionTypes.Where(ionType => ionType.IsPrefixIon).ToList();
            var suffixIonTypes = _ionTypes.Where(ionType => !ionType.IsPrefixIon).ToList();

            var probabilities = new Dictionary<IonType, IonProbability>();

            debugFile.WriteLine("Scan: {0}, Peptide: {1}", spectrumMatch.ScanNum, spectrumMatch.Peptide);
            debugFile.WriteLine("Prefixes");
            debugFile.WriteLine("Segment\tIon\tM/Z\tFound");
            for (int i = 1; i <= prefixes.Count; i++)
            {
                var segment = peptide.Substring(0, i);
                foreach (var ionType in prefixIonTypes)
                {
                    if (!probabilities.ContainsKey(ionType))
                        probabilities.Add(ionType, new IonProbability(ionType.Name));

                    var ion = ionType.GetIon(prefixes[i-1]);
                    double mz = ion.GetMonoIsotopicMz();
                    var present = spectrum.ContainsIon(ion, _tolerance, RelativeIntensityThreshold);
//                    var present = (spectrum.FindPeak(mz, _tolerance) != null);
                    probabilities[ionType].Total++;
                    if (present)
                        probabilities[ionType].Found++;

                    debugFile.WriteLine("{0}\t{1}\t{2}\t{3}", segment, ionType.Name, mz, present);
                }
            }

            debugFile.WriteLine("Suffixes");
            debugFile.WriteLine("Segment\tIon\tM/Z\tFound");
            for (int i = 1; i <= suffixes.Count; i++)
            {
                var segment = peptide.Substring(i, peptide.Length-i);
                foreach (var ionType in suffixIonTypes)
                {
                    if (!probabilities.ContainsKey(ionType))
                        probabilities.Add(ionType, new IonProbability(ionType.Name));

                    var ion = ionType.GetIon(suffixes[i-1]);
                    double mz = ion.GetMonoIsotopicMz();
                    var present = spectrum.ContainsIon(ion, _tolerance, RelativeIntensityThreshold);
//                    var present = (spectrum.FindPeak(mz, _tolerance) != null);
                    probabilities[ionType].Total++;
                    if (present)
                        probabilities[ionType].Found++;

                    debugFile.WriteLine("{0}\t{1}\t{2}\t{3}", segment, ionType.Name, mz, present);
                }
            }

            debugFile.WriteLine();
            debugFile.Close();
            return probabilities.Values.ToList();
        }

        [Test]
        public void IonPresent()
        {
            var spectrumMatchList = InitTest();

            StreamWriter outputFile = File.AppendText(OutputFileName);
            foreach (var spectrumMatch in spectrumMatchList)
            {
                var offsetFrequencies = new OffsetFrequencyTable();
                offsetFrequencies.AddCleavageProbabilities(new List<SpectrumMatch>{spectrumMatch}, _ionTypes, _tolerance,
                    RelativeIntensityThreshold);
                var probabilities = offsetFrequencies.IonProbabilityTable;

                var testProbabilities = ComputeOffsetFrequencies(spectrumMatch);

                Assert.True(probabilities.Count == testProbabilities.Count);

                probabilities.Sort(new CompareIonProbabilityByIon());
                testProbabilities.Sort(new CompareIonProbabilityByIon());

                outputFile.WriteLine("Scan: {0}, Peptide: {1}", spectrumMatch.ScanNum, spectrumMatch.Peptide);
                outputFile.WriteLine("Ion\tFound\tTotal\tTestFound\tTestTotal\tEqual?");
                for (int i = 0; i < probabilities.Count; i++)
                {
                    Assert.True(probabilities[i].IonName == testProbabilities[i].IonName);

                    bool foundEqual = probabilities[i].Found == testProbabilities[i].Found;
                    bool totalEqual = probabilities[i].Total == testProbabilities[i].Total;

                    outputFile.Write(probabilities[i].IonName + "\t");
                    outputFile.Write("{0}\t{1}\t", probabilities[i].Found, probabilities[i].Total);
                    outputFile.Write("{0}\t{1}\t", testProbabilities[i].Found, testProbabilities[i].Total);
                    outputFile.WriteLine(foundEqual && totalEqual);
                }
                outputFile.WriteLine();
            }
            outputFile.Close();
        }

        // Configuration
//        private const int ScanNum = 30623;
        private const double NoiseFiltration = 0;
        private const ActivationMethod Act = ActivationMethod.HCD;
        private const int MaxCharge = 1;
        const double RelativeIntensityThreshold = 0.8;

        private readonly BaseIonType[] _baseIons = { BaseIonType.B, BaseIonType.Y };
        private readonly NeutralLoss[] _neutralLosses = { NeutralLoss.NoLoss };
        private List<IonType> _ionTypes;
        private readonly Tolerance _tolerance = new Tolerance(15, ToleranceUnit.Ppm);

        private const string TsvFile = @"\\protoapps\UserData\Wilkins\ForChris\HCD_QCShew\tsv\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.tsv";
        private const string RawFile = @"\\protoapps\UserData\Wilkins\ForChris\HCD_QCShew\raw\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";
        private const string DebugFileName = @"C:\Users\wilk011\Documents\DataFiles\TestFolder\debug_RIT0.8_NoFiltration_2.txt";
        private const string OutputFileName = @"C:\Users\wilk011\Documents\DataFiles\TestFolder\results_RIT0.8_NoFiltration_2.txt";
        private List<SpectrumMatch> InitTest()
        {
            var ionTypeFactory = new IonTypeFactory(_baseIons, _neutralLosses, MaxCharge);

            _ionTypes = ionTypeFactory.GetAllKnownIonTypes().ToList();

            var lcms = LcMsRun.GetLcMsRun(RawFile, MassSpecDataType.XCaliburRun, NoiseFiltration, NoiseFiltration);

            var spectrumMatches = (new SpectrumMatchList(lcms, new TsvFileParser(TsvFile), Act)).Matches;

            var matches = spectrumMatches.Where(spectrumMatch => !spectrumMatch.Peptide.Contains("+")).ToList();

            foreach (var spectrumMatch in matches)
            {
                var aaset = new AminoAcidSet();

                var composition = aaset.GetComposition(spectrumMatch.Peptide);

                Assert.True(composition.Equals(spectrumMatch.PeptideComposition));
            }

            return matches;
        }
    }
}
