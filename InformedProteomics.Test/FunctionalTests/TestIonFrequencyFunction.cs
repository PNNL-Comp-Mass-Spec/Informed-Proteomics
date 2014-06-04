﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Data;
using InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestIonFrequencyFunction
    {
        public List<Probability<IonType>> ComputeOffsetFrequencies(SpectrumMatch spectrumMatch)
        {
            var debugFile = File.AppendText(DebugFileName);

            var peptide = spectrumMatch.Peptide;
            var sequence = spectrumMatch.Sequence;
            var spectrum = spectrumMatch.Spectrum;
            var prefixes = spectrumMatch.Prefixes;
            var suffixes = spectrumMatch.Suffixes;

            var prefixIonTypes = _ionTypes.Where(ionType => ionType.IsPrefixIon).ToList();
            var suffixIonTypes = _ionTypes.Where(ionType => !ionType.IsPrefixIon).ToList();

            var probabilities = new Dictionary<IonType, Probability<IonType>>();

            debugFile.WriteLine("Scan:\t{0}\tPeptide:\t{1}\tCharge:\t{2}", spectrumMatch.Spectrum.ScanNum, peptide, spectrumMatch.PrecursorCharge);
            debugFile.WriteLine("Prefixes");
            debugFile.WriteLine("Segment\tIon\tM/Z\tFound");

            for (int i = 0; i < prefixes.Count; i++)
            {
                var segment = spectrumMatch.Sequence.GetRange(0, i+1);
                var segmentStr = segment.Aggregate("", (current, aa) => current + aa.Residue);
                foreach (var ionType in prefixIonTypes)
                {
                    if (!probabilities.ContainsKey(ionType))
                        probabilities.Add(ionType, new Probability<IonType>(ionType));

                    var ion = ionType.GetIon(prefixes[i]);
                    double mz = ion.GetMonoIsotopicMz();
                    var present = spectrum.ContainsIon(ion, _tolerance, RelativeIntensityThreshold);
                    probabilities[ionType].Total++;
                    if (present)
                        probabilities[ionType].Found++;

                    debugFile.WriteLine("{0}\t{1}\t{2}\t{3}", segmentStr, ionType.Name, mz, present);
                }
            }

            debugFile.WriteLine("Suffixes");
            debugFile.WriteLine("Segment\tIon\tM/Z\tFound");

            for (int i = 0; i < suffixes.Count; i++)
            {
                var segment = spectrumMatch.Sequence.GetRange(i+1, sequence.Count-(i+1));
                var segmentStr = segment.Aggregate("", (current, aa) => current + aa.Residue);
                foreach (var ionType in suffixIonTypes)
                {
                    if (!probabilities.ContainsKey(ionType))
                        probabilities.Add(ionType, new Probability<IonType>(ionType));

                    var ion = ionType.GetIon(suffixes[i]);
                    double mz = ion.GetMonoIsotopicMz();
                    var present = spectrum.ContainsIon(ion, _tolerance, RelativeIntensityThreshold);
                    probabilities[ionType].Total++;
                    if (present)
                        probabilities[ionType].Found++;

                    debugFile.WriteLine("{0}\t{1}\t{2}\t{3}", segmentStr, ionType.Name, mz, present);
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

            var ionProbabilityTable = new Dictionary<IonType, Probability<IonType>>[MaxPrecCharge];
            for (int i = 0; i < ionProbabilityTable.Count(); i++)
            {
                ionProbabilityTable[i] = new Dictionary<IonType, Probability<IonType>>();
            }
            StreamWriter outputFile = File.AppendText(OutputFileName);
            foreach (var spectrumMatch in spectrumMatchList)
            {
                var charge = spectrumMatch.PrecursorCharge - 1;
                var offsetFrequencies = new IonFrequencyTable(_ionTypes, _tolerance, RelativeIntensityThreshold);
                offsetFrequencies.AddMatches(new List<SpectrumMatch>{spectrumMatch});
                var probabilities = offsetFrequencies.GetProbabilities();

                var testProbabilities = ComputeOffsetFrequencies(spectrumMatch);

                Assert.True(probabilities.Length == testProbabilities.Count);

                Array.Sort(probabilities, new CompareIonProbabilityByIon());
                testProbabilities.Sort(new CompareIonProbabilityByIon());

                outputFile.WriteLine("Scan:\t{0}\tPeptide:\t{1}\tCharge:\t{2}", spectrumMatch.Spectrum.ScanNum, spectrumMatch.Peptide, spectrumMatch.PrecursorCharge);
                outputFile.WriteLine("Ion\tFound\tTotal\tTestFound\tTestTotal\tEqual?");
                for (int i = 0; i < probabilities.Length; i++)
                {
                    Assert.True(probabilities[i].Label.Name == testProbabilities[i].Label.Name);

                    bool foundEqual = probabilities[i].Found.Equals(testProbabilities[i].Found);
                    bool totalEqual = probabilities[i].Total.Equals(testProbabilities[i].Total);

                    var ion = testProbabilities[i].Label;
                    if (!ionProbabilityTable[charge].ContainsKey(ion))
                        ionProbabilityTable[charge].Add(ion, new Probability<IonType>(ion));
                    ionProbabilityTable[charge][ion] += testProbabilities[i];

                    outputFile.Write(probabilities[i].Label + "\t");
                    outputFile.Write("{0}\t{1}\t", probabilities[i].Found, probabilities[i].Total);
                    outputFile.Write("{0}\t{1}\t", testProbabilities[i].Found, testProbabilities[i].Total);
                    outputFile.WriteLine(foundEqual && totalEqual);
                }
                outputFile.WriteLine();
            }
            outputFile.Close();
            for (int i = 0; i < MaxPrecCharge; i++)
            {
                var chargeOutName = IonProbabilityFileName + "_Charge" + (i + 1) +".txt";
                using (var ionProbOut = new StreamWriter(chargeOutName))
                {
                    ionProbOut.WriteLine("Ion\tProbability");
                    foreach (var key in ionProbabilityTable[i].Keys)
                    {
                        ionProbOut.WriteLine(ionProbabilityTable[i][key].Label.Name + "\t" + ionProbabilityTable[i][key].Prob);
                    }
                }
            }
        }

/*        [Test]
        public void TestSequence()
        {
            var matches = InitTest().ToList();
            var match = matches[0];
            var sequence = match.Sequence;
            var prefixes = match.Prefixes;
            var suffixes = match.Suffixes;
            foreach (var aa in sequence)
            {
                Console.WriteLine(aa.Residue+"\t"+aa.GetMass());
            }

            Console.WriteLine("Prefixes");
            for (int i=0; i<prefixes.Count; i++)
            {
                Console.WriteLine((i+1)+"\t"+prefixes[i].GetIsotopeMass(0));
            }
            Console.WriteLine("Suffixes");
            for (int i = 0; i < suffixes.Count; i++)
            {
                Console.WriteLine((i + 1) + "\t" + suffixes[i].GetIsotopeMass(0));
            }
        }*/

        // Configuration
//        private const int ScanNum = 30623;
        private const double NoiseFiltration = 0;
        private const ActivationMethod Act = ActivationMethod.HCD;
        private const int MaxCharge = 1;
        private const int MaxPrecCharge = 4;
        const double RelativeIntensityThreshold = 1.0;

        private readonly BaseIonType[] _baseIons =
        {
            BaseIonType.B, BaseIonType.Y
        };
        private readonly NeutralLoss[] _neutralLosses = { NeutralLoss.NoLoss, NeutralLoss.H2O };
        private List<IonType> _ionTypes;
        private readonly Tolerance _tolerance = new Tolerance(0.5, ToleranceUnit.Th);

        private const string TsvFile = @"\\protoapps\UserData\Wilkins\BottomUp\HCD_QCShew\tsv\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.tsv";
        private const string RawFile = @"\\protoapps\UserData\Wilkins\MSGFPlusTrainingData\CID_LowRes_Tryp.mgf";
        private const string DebugFileName = @"C:\Users\wilk011\Documents\DataFiles\TestFolder\debug_CID_LowRes_Tryp.txt";
        private const string OutputFileName = @"C:\Users\wilk011\Documents\DataFiles\TestFolder\results_CID_LowRes_Tryp.txt";

        private const string IonProbabilityFileName =
            @"C:\Users\wilk011\Documents\DataFiles\TestFolder\IonProbabilities_CID_LowRes_Tryp";
        private IEnumerable<SpectrumMatch> InitTest()
        {
            var ionTypeFactory = new IonTypeFactory(_baseIons, _neutralLosses, MaxCharge);
            _ionTypes = new List<IonType> {ionTypeFactory.GetIonType("b"), ionTypeFactory.GetIonType("y-H2O")};

//            _ionTypes = ionTypeFactory.GetAllKnownIonTypes().ToList();
            

//            var lcms = LcMsRun.GetLcMsRun(RawFile, MassSpecDataType.XCaliburRun, NoiseFiltration, NoiseFiltration);

//            var spectrumMatches = (new SpectrumMatchList(lcms, new TsvFileParser(TsvFile), Act, false, MaxPrecCharge));

            var spectrumMatches = new SpectrumMatchList(RawFile, MaxPrecCharge);
//            spectrumMatches.FilterSpectra();
            return spectrumMatches;
        }
    }

    class CompareIonProbabilityByIon : IComparer<Probability<IonType>>
    {
        public int Compare(Probability<IonType> x, Probability<IonType> y)
        {
            return String.Compare(x.Label.Name, y.Label.Name, StringComparison.Ordinal);
        }
    }
}
