using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text.RegularExpressions;
using System.Xml.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Tests.Base;
using MathNet.Numerics.Statistics;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestUtils
    {
        [Test]
        [Category("Local_Testing")]
        public void TestCountingMs2Scans()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var numSpectra = 0;
            const string rawFileDir = @"H:\Research\TopDownTestData";

            if (!Directory.Exists(rawFileDir))
            {
                Assert.Ignore(@"Skipping test " + methodName + @" since folder not found: " + rawFileDir);
            }

            foreach (var rawFilePath in Directory.GetFiles(rawFileDir, "*.raw"))
            {
                var run = PbfLcMsRun.GetLcMsRun(rawFilePath);
                numSpectra += run.GetScanNumbers(2).Count;
            }
            Console.WriteLine(numSpectra);
        }

        [Test]
        [TestCase(100, 15000, 29, 134946816, 27434, 952705, 4087595)]
        [TestCase(200, 50000, 29, 135077888, 54867, 1043777, 13661895)]
        [TestCase(200, 25000, 29, 135077888, 54867, 912705, 6803514)]
        [TestCase(500, 25000, 29, 135260160, 137168, 730433, 6721213)]
        public void TestFloatBinning(
            double minMass, double maxMass, int numBits,
            int expectedStartBin1, int expectedStartBin2,
            int expectedBinCount1, int expectedBinCount2)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var comparer = new MzComparerWithBinning(numBits);

            var minMassBin1 = comparer.GetBinNumber(minMass);
            var maxMassBin1 = comparer.GetBinNumber(maxMass);

            var minMassBin2 = LcMsMatchMap.GetBinNumber(minMass);
            var maxMassBin2 = LcMsMatchMap.GetBinNumber(maxMass);

            var binCount1 = maxMassBin1 - minMassBin1 + 1;
            var binCount2 = maxMassBin2 - minMassBin2 + 1;

            Console.WriteLine("MzComparer   bin range for {0} to {1} with {2} bits: {3} to {4}", minMass, maxMass, numBits, minMassBin1, maxMassBin1);
            Console.WriteLine("LcMsMatchMap bin range for {0} to {1} with {2} bits: {3} to {4}", minMass, maxMass, numBits, minMassBin2, maxMassBin2);

            Console.WriteLine("MzComparer   binCount: {0}", binCount1);
            Console.WriteLine("LcMsMatchMap binCount: {0}", binCount2);

            Assert.AreEqual(expectedStartBin1, minMassBin1, "BinCount mismatch for MzComparerWithBinning");
            Assert.AreEqual(expectedStartBin2, minMassBin2, "BinCount mismatch for MzComparerWithBinning");

            Assert.AreEqual(expectedBinCount1, binCount1, "BinCount mismatch for MzComparerWithBinning");
            Assert.AreEqual(expectedBinCount2, binCount2, "BinCount mismatch for LcMsMatchMap");

        }

        [Test]
        [TestCase(3, 4, 81, 15)]
        [TestCase(5, 7, 78125, 330)]
        public void TestGeneratingNtoKCombinationsWithRepetition(int n, int k, int expectedCount, int expectedCount2)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var combinations = SimpleMath.GetNtoTheKCombinations(n, k);
            var count = 0;
            foreach (var combination in combinations)
            {
                ++count;
                if (count < 15)
                    Console.WriteLine((count - 1) + ": " + string.Join(",", combination));
                else if (count == 15)
                    Console.WriteLine("...");
            }

            var count2 = SimpleMath.GetCombination(n + k - 1, k);

            Console.WriteLine("Count: " + count);
            Console.WriteLine("Count2: " + count2);

            Assert.AreEqual(expectedCount, count, "Mismatch for number of items returned by GetNtoTheKCombinations");
            Assert.AreEqual(expectedCount2, count2, 0.0001, "Mismatch for SimpleMath.GetCombination");
        }

        [Test]
        [TestCase(5, 3, 35, 35)]
        [TestCase(5, 7, 330, 330)]
        public void TestGeneratingCombinations(int n, int k, int expectedCount, int expectedCount2)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var combinations = SimpleMath.GetCombinationsWithRepetition(n, k);
            var count = 0;
            foreach (var combination in combinations)
            {
                ++count;
                if (count < 15)
                    Console.WriteLine((count - 1) + ": " + string.Join(",", combination));
                else if (count == 15)
                    Console.WriteLine("...");
            }

            var count2 = SimpleMath.GetCombination(n + k - 1, k);

            Console.WriteLine("Count: " + count);
            Console.WriteLine("Count2: " + count2);

            Assert.AreEqual(expectedCount, count, "Mismatch for number of items returned by GetCombinationsWithRepetition");
            Assert.AreEqual(expectedCount2, count2, 0.0001, "Mismatch for SimpleMath.GetCombination");
        }

        [Test]
        [TestCase("+229.163C+57.021GLGGSGTPVDELDK+229.163C+57.021C+57.021QTHDNC+57.021YDQAK+229.163", 968.95333373)]
        public void ParseMsGfString(string msgfPepStr, double expectedMonoMz)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(msgfPepStr);
            var ion = new Ion(sequence.Composition + Composition.H2O, 4);
            var computedMonoMz = ion.GetMonoIsotopicMz();

            Console.WriteLine("{0} {1} {2}", msgfPepStr, ion.Composition, computedMonoMz);

            Assert.AreEqual(expectedMonoMz, computedMonoMz, 0.000001, "Mass mismatch for {0}", msgfPepStr);
        }

        [Test]
        [TestCase("H(250) C(150) N(40) O(50) S(2)", 3475.76909165)]
        [TestCase("H(225) C(150) N(40) O(50) S(2) H(25)", 3475.76909165)]
        [TestCase("H(257) C(150) N(42) O(56) S(4) 13C(12) 15N(3)", 3871.784229)]
        [TestCase("H(225) C(150) N(40) O(50) S(2) Fe(3)", 3618.378283675)]
        [TestCase("H(225) C(150) N(40) O(50) S(2) Fe(2) Fe(1)", 3618.378283675)]
        public void ParseComposition(string empiricalFormula, double expectedMass)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var comp = Composition.Parse(empiricalFormula);
            var computedMass = comp.Mass;
            Console.WriteLine("AveragineMass: {0}", computedMass);

            Assert.AreEqual(expectedMass, computedMass, 0.000001, "Mass mismatch for {0}", empiricalFormula);
        }

        [Test]
        [TestCase("15.9949", "Oxidation", "15.995")]
        [TestCase("229.163", "TMT6plex", "229.163")]
        [TestCase("123456", "", "")]
        public void TestGetModFromMass(string modMassText, string expectedModName, string expectedFormattedMass)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var matchingMods = Modification.GetFromMass(modMassText);
            if (matchingMods == null || matchingMods.Count == 0)
            {
                Console.WriteLine("Mod mass of {0} does not resolve to any mods", modMassText);
                Assert.AreEqual("", expectedModName, "Expected mod to resolve to {0} but did not resolve", expectedModName);
                return;
            }

            var matchingMod = matchingMods.First();
            var modMassFormatted = string.Format("{0:N3}", matchingMod.Mass);

            Console.WriteLine("Mod mass of {0} resolves to {1}, mass ", modMassText, matchingMod.Name, modMassFormatted);

            Assert.AreEqual(expectedFormattedMass, modMassFormatted);
        }

        /// <summary>
        /// Test formatting of modification mass values for various modifications
        /// </summary>
        /// <remarks>Does not use TestCase since we cannot reference a Modification instance at design time</remarks>
        [Test]
        public void TestFormattingByModName()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var expectedMassByMod = new Dictionary<Modification, string> {
                { Modification.Phosphorylation, "79.966"},
                { Modification.Tmt6Plex, "229.163"},
                { Modification.TriMethylation, "42.047"}
            };

            foreach (var mod in expectedMassByMod)
            {

                var modMass = mod.Key.Mass;
                var modMassFormatted = string.Format("{0:N3}", modMass);

                Console.WriteLine("Modification mass of {0,-9} is {1,7}", mod.Key.Name, modMassFormatted);

                Assert.AreEqual(mod.Value, modMassFormatted);
            }

        }

        [Test]
        [TestCase(5105.8432, 7, 8, 730.413448, 0.2473681, 730.8434571)]
        [TestCase(10247.529329, 12, 11, 855.1352797, 0.271738, 855.469731)]
        [TestCase(10247.529329, 14, 11, 733.116993, 0.271738, 733.403666)]
        public void TestAveragine(
            double monoMass, int charge, int expectedPeaks,
            double expectedFirstPeakMz, double expectedFirstPeakIntensity,
            double expectedMzHighestIntensity)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var profile = Averagine.GetTheoreticalIsotopeProfile(monoMass, charge);

            Console.WriteLine("Isotope ions:");

            var highestIntensity = 0.0;
            var highestIntensityMz = 0.0;

            foreach (var p in profile)
            {
                var plotText = new string('+', (int)(p.Intensity * 100 / 5));
                Console.WriteLine("{0,9:N4} {1,8:N5} {2}", p.Mz, p.Intensity, plotText);

                if (p.Intensity > highestIntensity)
                {
                    highestIntensity = p.Intensity;
                    highestIntensityMz = p.Mz;
                }
            }
            Console.WriteLine();

            Assert.AreEqual(expectedPeaks, profile.Count, "Peak count mismatch");
            Assert.AreEqual(expectedFirstPeakMz, profile.First().Mz, 0.00001, "First peak, unexpected m/z");
            Assert.AreEqual(expectedFirstPeakIntensity, profile.First().Intensity, 0.00001, "First peak, unexpected intensity");
            Assert.AreEqual(expectedMzHighestIntensity, highestIntensityMz, 0.00001, "Peak count mismatch");
        }

        [Test]
        public void TestAaSet()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // Note that SequenceLocation.Everywhere is not supported when TargetResidue is '*'
            //var oxEverywhere = new SearchModification(Modification.Oxidation, '*', SequenceLocation.Everywhere, true);

            var acetylNTerm = new SearchModification(Modification.Acetylation, '*', SequenceLocation.PeptideNTerm, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);

            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                acetylNTerm,
                glutathioneC,
                oxM
            };

            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);
            aaSet.Display();

        }

        [Test]
        [TestCase("IRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG",
            "C(185) H(299) N(55) O(53) S(1)", 13, 7, 4171.211301, 4172.214656, 0.370361)]
        [TestCase("METTKPSFQDVLEFVRLFRRKNKLQREIQDVEKKIRDNQKRVLLLDNLSDYIKPGMSVEAIQGIIASMKGDYEDRVDDYIIKNAELSKERRDISKKLKAMGEMKNGEAK",
            "C(557) H(929) N(161) O(171) S(5)", 18, 15, 12769.755127, 12770.758481, 0.004518)]
        public void TestPeptide(
            string sequence, string expectedEmpiricalFormula, int chargeState, int expectedPeakCount,
            double expectedMass, double expectedSecondIsotopeMz, double expectedFirstIsotopomerIntensity)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var aaSet = new AminoAcidSet();
            var composition = aaSet.GetComposition(sequence) + Composition.H2O;

            Console.WriteLine("Empirical formula: " + composition);
            Console.WriteLine("Monoisotopic mass: {0:F4}", composition.Mass);
            Console.WriteLine("Nominal mass:      {0}", composition.NominalMass);

            Assert.AreEqual(expectedEmpiricalFormula, composition.ToString());

            // Remove parentheses and spaces from the expected empirical formula, then compare to composition.ToPlainString
            var expectedPlainString = Regex.Replace(expectedEmpiricalFormula, "[() ]", "");
            Assert.AreEqual(expectedPlainString, composition.ToPlainString());

            Assert.AreEqual(expectedMass, composition.Mass, 0.00001, "Monoisotopic mass mismatch");

            // 2nd isotope
            Console.WriteLine();
            Console.WriteLine("First 3 isotopes");
            Console.WriteLine("{0:5F}" ,composition.GetIsotopeMass(0));
            Console.WriteLine("{0:5F}", composition.GetIsotopeMass(1));
            Console.WriteLine("{0:5F}", composition.GetIsotopeMass(2));

            Assert.AreEqual(expectedSecondIsotopeMz, composition.GetIsotopeMass(1), 0.00001, "Second isotope m/z mismatch");


            Console.WriteLine();
            Console.WriteLine("Isotopomer Envelope:");

            var ionIndex = 0;
            var intensities = composition.GetIsotopomerEnvelopeRelativeIntensities();

            var isotopes = new List<Isotope>();
            foreach (var intensity in intensities)
            {
                isotopes.Add(new Isotope(ionIndex, intensity));
                ionIndex++;
            }

            ShowIsotopes(isotopes);

            Assert.AreEqual(expectedPeakCount, intensities.Length, "Isotopomer envelope count mismatch");
            Assert.AreEqual(expectedFirstIsotopomerIntensity, intensities.First(), 0.00001, "First isotopomer envelope peak intensity mismatch");
            Console.WriteLine();

            Console.WriteLine("Isotope ions for charge {0}: ", chargeState);
            var ion = new Ion(composition + Composition.H2O, chargeState);

            ShowIsotopes(ion.GetIsotopes(0.1));

            Console.WriteLine();
        }

        [Test]
        [TestCase("C(82) H(149) N(23) O(24) S(3)", 3, 5, 1936.030795435, 646.350875, 0.960131)]
        [TestCase("C(82) H(149) N(23) O(24) S(3)", 6, 5, 1936.030795435, 323.679076, 0.960131)]
        [TestCase("C(212) H(350) N(50) O(32) S(1)", 6, 7, 4140.70180111, 691.124243, 0.373737)]
        [TestCase("C(212) H(350) N(50) O(32) Cl(1)", 6, 7, 4143.69858313, 691.623707, 0.373737)]
        public void TestIsoProfile(
            string empiricalFormula, int charge, int expectedCount,
            double expectedMass, double expectedMonoMz, double expectedFirstIsotopeIntensity)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var composition = Composition.Parse(empiricalFormula);
            var ion = new Ion(composition, charge);

            Console.WriteLine("Mass {0:F4} Da; {1:F4} m/z for charge {2}", composition.Mass, ion.GetIsotopeMz(0), charge);

            Assert.AreEqual(expectedMass, composition.Mass, 0.00001, "Mono mass mismatch");
            Assert.AreEqual(expectedMonoMz, ion.GetIsotopeMz(0), 0.00001, "Ion m/z mismatch");

            var isotopes = ion.GetIsotopes(0.1).ToList();
            ShowIsotopes(isotopes);

            Assert.AreEqual(expectedFirstIsotopeIntensity, isotopes.First().Ratio, 0.00001, "First isotope intensity mismatch");

        }


        [Test]
        [TestCase("CCAADDKEACFAVEGPK", 2, 2, "C(78) H(122) N(22) O(29) S(3)", 964.4027858, 136,
            "b,3,463.142800|y,14,1536.699885|b-H2O,10,1310.459854|y2,10,554.260599")]
        [TestCase("IRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG", 3, 2, "C(185) H(299) N(55) O(53) S(1)", 1391.4110436, 296,
            "b,3,456.256508|y,16,1785.975864|b-H2O, 8,1031.515636|y2,12,670.365110")]
        public void TestGetProductIons(
            string peptide, int precursorCharge, int maxCharge,
            string expectedEmpiricalFormula, double expectedMz, int expectedProductIons,
            string expectedProductIonMasses)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var sequence = new Sequence(peptide, aaSet);

            var ionTypeFactory = new IonTypeFactory(
                new[] {BaseIonType.B, BaseIonType.Y},
                new[] {NeutralLoss.NoLoss, NeutralLoss.H2O},
                maxCharge);

            var precursor = sequence.GetPrecursorIon(precursorCharge);
            var monoIsotopicMz = precursor.GetMonoIsotopicMz();

            Console.WriteLine("Precursor Ion: {0}, mass {1:F4}", precursor.Composition, monoIsotopicMz);

            Assert.AreEqual(expectedEmpiricalFormula, precursor.Composition.ToString(), "Empirical formula mismatch");
            Assert.AreEqual(expectedMz, monoIsotopicMz, 0.00001, "Mono m/z mismatch");

            Console.WriteLine("Product ions: ");
            Console.WriteLine("{0,-10} {1,2}  {2,-30} {3,6} {4}", "IonType", "N", "Composition", "Charge", "m/z");

            var productIons = sequence.GetProductIons(ionTypeFactory.GetAllKnownIonTypes());

            // Show the first 3 and the last 3 product ions of each ion type
            var ionCache = new Dictionary<IonType, List<Ion>>();

            foreach (var theoIon in productIons)
            {
                var ionTypeAndIndex = theoIon.Key;
                var ionType = ionTypeAndIndex.Item1;
                // var index = ionTypeAndIndex.Item2;
                var ion = theoIon.Value;

                if (ionCache.TryGetValue(ionType, out var ionsForIonType))
                {
                    ionsForIonType.Add(ion);
                }
                else
                {
                    ionCache.Add(ionType, new List<Ion> { ion });
                }
            }

            foreach (var ionType in (from item in ionCache.Keys orderby item.Name select item))
            {
                var ionList = ionCache[ionType];
                if (ionList.Count <= 6)
                    ShowProductIons(ionType, 0, ionList);
                else
                {
                    ShowProductIons(ionType, 0, ionList.Take(3));
                    Console.WriteLine("...");
                    var startIndex = ionList.Count - 3;
                    ShowProductIons(ionType, startIndex, ionList.Skip(startIndex));
                }
                Console.WriteLine();
            }

            foreach (var expectedIon in expectedProductIonMasses.Split('|'))
            {
                var tokens = expectedIon.Split(',');

                if (tokens.Length != 3)
                {
                    Assert.Fail("Invalid format for expectedProductIon {0}; should be IonTypeName,IonNumber,Mass", expectedIon);
                }

                var ionName = tokens[0].Trim();
                if (!int.TryParse(tokens[1], out var ionIndex))
                {
                }

                if (!double.TryParse(tokens[2], out var expectedIonMass))
                {
                    Assert.Fail("Mass not numeric in expectedProductIon {0}", expectedIon);
                }

                // Look for ionName in ionCache

                var ionTypeToFind = ionTypeFactory.GetIonType(ionName);

                if (!ionCache.TryGetValue(ionTypeToFind, out var observedProductIons))
                {
                    Assert.Fail("Observed product ions do not have ion type " + ionTypeToFind.Name);
                }

                var observedProductIon = observedProductIons[ionIndex];
                Assert.AreEqual(expectedIonMass, observedProductIon.GetMonoIsotopicMz(), 0.00001, "Product ion mass mismatch for {0}{1}", ionName, ionIndex);

            }

            Assert.AreEqual(expectedProductIons, productIons.Count, "Product ion count mismatch");

        }

        [Test]
        [TestCase("_.STR._",      2, 3,       "C(13) H(24) N(6) O(5) S(0)",    "C(13) H(26) N(6) O(11) S(0) P(2)")]
        [TestCase("_.PEPTIDES._", 2, 3,       "C(37) H(56) N(8) O(16) S(0)",   "C(37) H(58) N(8) O(22) S(0) P(2)")]
        [TestCase("_.PEPTIDESTRINGS._", 2, 3, "C(62) H(100) N(18) O(25) S(0)", "C(62) H(102) N(18) O(31) S(0) P(2)")]
        [TestCase("_.PEPTIDESTRINGS._", 4, 5, "C(62) H(100) N(18) O(25) S(0)", "C(62) H(104) N(18) O(37) S(0) P(4)")]
        public void TestSequenceGraph(
            string annotation, int maxModsPerPeptide,
            int expectedCompositionCount, string compositionFirst, string compositionLast)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var phosPhoS = new SearchModification(Modification.Phosphorylation, 'S', SequenceLocation.Everywhere, false);
            var phosPhoT = new SearchModification(Modification.Phosphorylation, 'T', SequenceLocation.Everywhere, false);
            var phosPhoY = new SearchModification(Modification.Phosphorylation, 'Y', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var fixCarbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C', SequenceLocation.Everywhere, true);

            var searchModifications = new List<SearchModification> { phosPhoS, phosPhoT, phosPhoY, oxM, fixCarbamidomethylC };

            var aaSet = new AminoAcidSet(searchModifications, maxModsPerPeptide);

            var pepSeq = annotation.Substring(2, annotation.Length - 4);

            var graph = SequenceGraph.CreateGraph(aaSet, annotation);

            var unmodifiedSequenceComp = graph.GetUnmodifiedSequenceComposition();

            Console.WriteLine("Annotation Compositions for {0} with unmodified formula {1} ", annotation, unmodifiedSequenceComp);

            Assert.AreEqual(unmodifiedSequenceComp, aaSet.GetComposition(pepSeq));

            var compositions = graph.GetSequenceCompositions();

            var index = 0;
            foreach (var composition in compositions)
            {
                Console.WriteLine("{0}: {1}", index + 1, composition);
                index++;
            }

            Assert.AreEqual(expectedCompositionCount, compositions.Length, "Composition count mismatch");
            Assert.AreEqual(compositionFirst, compositions.First().ToString(), "First composition mismatch");
            Assert.AreEqual(compositionLast, compositions.Last().ToString(), "Last composition mismatch");

            //const int seqIndex = 1;
            //Console.WriteLine("Fragment Compositions (" + seqIndex +")");
            //var scoringGraph = graph.GetScoringGraph(seqIndex);
            //foreach (var composition in scoringGraph.GetCompositions())
            //{
            //    Console.WriteLine(composition);
            //}
        }

        /// <summary>
        /// Test computing combinations of modifications
        /// </summary>
        /// <remarks>Does not use TestCase since we cannot reference a Modification instance at design time</remarks>
        [Test]
        public void TestModificationParams()
        {
            const int maxNumDynModsPerSequence = 3;

            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // Keys in the dictionary are the list of modifications to use
            // Values are tuples of the expected result:
            //  ComboCount, MinMass, MaxMass, Result_for_GetModificationCombinationIndex(8, 0), Result_for_GetModificationCombinationIndex(12, 1)
            var modCombos = new Dictionary<List<Modification>, Tuple<int, double, double, int, int>>();

            AppendModCombo(modCombos,
                           new List<Modification> {Modification.Acetylation, Modification.Phosphorylation, Modification.Oxidation},
                           new Tuple<int, double, double, int, int>(20, 0, 239.898993, 14, -1));

            AppendModCombo(modCombos,
                           new List<Modification> { Modification.Acetylation, Modification.Phosphorylation, Modification.Oxidation, Modification.PyroGluQ },
                           new Tuple<int, double, double, int, int>(35, -51.07965, 239.898993, 18, 28));

            foreach (var combo in modCombos)
            {
                var modifications = combo.Key;
                var expectedResults = combo.Value;

                var expectedComboCount = expectedResults.Item1;
                var expectedMassMin = expectedResults.Item2;
                var expectedMassMax = expectedResults.Item3;
                var expectedModCombo1 = expectedResults.Item4;
                var expectedModCombo2 = expectedResults.Item5;

                var modParams = new ModificationParams(modifications.ToArray(), maxNumDynModsPerSequence);

                var comboList = (from item in modifications select item.Name);
                Console.WriteLine("Combinations for " + string.Join(", ", comboList));

                var numCombinations = modParams.NumModificationCombinations;
                var minMass = double.MaxValue;
                var maxMass = double.MinValue;

                for (var modCombIndex = 0; modCombIndex < numCombinations; modCombIndex++)
                {
                    var modCombination = modParams.GetModificationCombination(modCombIndex);
                    var comboMass = modCombination.Composition.Mass;

                    if (modCombIndex < 3 || modCombIndex >= numCombinations - 3)
                    {
                        Console.WriteLine("{0}: Mods {1,-25}, formula {2,-25}, mass {3:F4}",
                                          modCombIndex, modCombination, modCombination.Composition, comboMass);
                    } else if (modCombIndex == 3)
                        Console.WriteLine("...");

                    minMass = Math.Min(minMass, comboMass);
                    maxMass = Math.Max(maxMass, comboMass);
                }

                var modCombo1 = modParams.GetModificationCombinationIndex(8, 0);
                var modCombo2 = modParams.GetModificationCombinationIndex(12, 1);

                Console.WriteLine("Found {0} combos, min mass {1:F4}, max mass {2:F4}", numCombinations, minMass, maxMass);
                Console.WriteLine("CombinationIndex(8, 0) is {0} and CombinationIndex(12, 0) is {1}", modCombo1, modCombo2);

                Assert.AreEqual(expectedComboCount, numCombinations, "Combination count mismatch");

                Assert.AreEqual(expectedMassMin, minMass, 0.00001, "Min mass mismatch");
                Assert.AreEqual(expectedMassMax, maxMass, 0.00001, "Max mass mismatch");

                Assert.AreEqual(expectedModCombo1, modCombo1, "Mismatch for mod combo1 (8,0)");
                Assert.AreEqual(expectedModCombo2, modCombo2, "Mismatch for mod combo1 (12,1)");

                Console.WriteLine();
            }

        }

        [Test]
        public void TestCSharpSyntax()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            //const string path = @"C:\cygwin\home\kims336\Developments\InformedProteomics\InformedProteomics.Test\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            //Console.WriteLine(Path.GetFileName(path));
            //Console.WriteLine(Path.GetFileNameWithoutExtension(path));
            //Console.WriteLine(Path.GetExtension(path));
            //Console.WriteLine(Path.GetDirectoryName(path));

            //const int size = int.MaxValue/4-1;
            //var hugeList = new List<int>(size);
            //Console.WriteLine("Success: " + size + " " + hugeList.Capacity);

            var set = new HashSet<double>();
            Console.WriteLine("Max: " + set.DefaultIfEmpty().Max(n => n*2));
        }

        [Test]
        [TestCase("a.",     1,  -26.987089,  "C(-1) H(1) N(0) O(-1) S(0)")]
        [TestCase("b",      1,    0,         "C(0) H(0) N(0) O(0) S(0)")]
        [TestCase("y",      1,   18.010565,  "C(0) H(2) N(0) O(1) S(0)")]
        [TestCase("b2-H2O", 2,  -18.010565,  "C(0) H(-2) N(0) O(-1) S(0)")]
        [TestCase("c-H2O",  1,   -0.9840156, "C(0) H(1) N(1) O(-1) S(0)")]
        [TestCase("z3",     3,    1.9918406, "C(0) H(0) N(-1) O(1) S(0)")]
        [TestCase("w2",     2,   73.0290,    "C(0) H(0) N(0) O(0) S(0) 73.029")]
        [TestCase("b3-NH3", 3,  -17.0265491, "C(0) H(-3) N(-1) O(0) S(0)")]
        public void TestIonTypeGeneration(string ionTypeName, int expectedCharge, double expectedMass, string expectedComposition)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var printTable = (ionTypeName == "b");
            var matchFound = false;

            if (printTable)
                Console.WriteLine("{0,2} {1,-10} {2,-40} {3,-12} {4,4} {5,6}",
                                  "#", "Name", "Composition", "      Mass", "Charge", "IsPrefixIon");

            var ionTypeFactory = new IonTypeFactory();
            var index = 0;
            foreach (var ionType in ionTypeFactory.GetAllKnownIonTypes())
            {
                var validateIonInfo = (ionType.Name == ionTypeName);

                index++;
                if (printTable)
                {
                    Console.WriteLine("{0,2} {1,-10} {2,-40} {3,12:F4} {4,4} {5,6}",
                                      index, ionType.Name, ionType.OffsetComposition, ionType.Mass, ionType.Charge, ionType.IsPrefixIon);
                } else if (validateIonInfo)
                {
                    Console.WriteLine("Ion {0} with composition {1}, mass {2:F4}, charge {3}+",
                                      ionType.Name, ionType.OffsetComposition, ionType.Mass, ionType.Charge);
                }

                if (!validateIonInfo)
                    continue;

                matchFound = true;

                Assert.AreEqual(expectedCharge, ionType.Charge, "Charge mismatch");
                Assert.AreEqual(expectedMass, ionType.Mass, 0.00001, "Mass mismatch");
                Assert.AreEqual(expectedComposition, ionType.OffsetComposition.ToString(), "Composition mismatch");

                if (!printTable)
                {
                    break;
                }
            }

            Assert.True(matchFound, "Match not found to ion {0}", ionTypeName);

        }

        [Test]
        public void TestIsotopemerProfileByKyowon() // is faster and more accurate than IsotopicDistributionCalculator
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            //C78H120N22O28S3 C150H120N220O28S30
            //var additionalElements = new[]
            //    {
            //        new Tuple<Atom, short>(Atom.Get("P"), 1),
            //        new Tuple<Atom, short>(Atom.Get("13C"), 3),
            //        new Tuple<Atom, short>(Atom.Get("15N"), 1),
            //    };
            //var composition = new Composition(149, 244, 44, 57, 0, additionalElements);
//            var composition = new Composition(83, 136, 22, 24, 1);
            //var composition = new Composition(210, 323, 54, 61, 0);
            var composition = new Composition(419, 699, 119, 129, 1);
            const int charge = 14;
            var ion = new Ion(composition + Composition.H2O, charge);
            var ff = composition.GetIsotopomerEnvelopeRelativeIntensities();
            var isotopeIndex = -1;
            foreach (var ii in ff)
            {
                ++isotopeIndex;
                Console.WriteLine("{0,2}:  {1:F4}  {2:F4}", isotopeIndex, ion.GetIsotopeMz(isotopeIndex), ii);
            }
        }

        [Test]
        public void TestReadingPnnlOmicsXmlFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var xmlFileName = Utils.GetTestFile(methodName, Path.Combine(Utils.DEFAULT_TEST_FILE_FOLDER, "PNNLOmicsElementData.xml"));

            var xdocument = XDocument.Load(xmlFileName.FullName);
            var parameterBaseElement = xdocument.Element("parameters");

            if (parameterBaseElement == null)
            {
                throw new IOException("Problem reading xml file " + xmlFileName + "; Expected element 'parameters' but it was not found");
            }

// ReSharper disable PossibleNullReferenceException
            var elements = parameterBaseElement.Element("ElementIsotopes").Elements("Element");
// ReSharper restore PossibleNullReferenceException
            foreach (var element in elements)
            {
                var symbol = element.Element("Symbol");
                var name = element.Element("Name");
                if (symbol != null && name != null)
                {
                    Console.WriteLine("{0,3}\t{1}", symbol.Value, name.Value);
                }
            }
        }

        [Test]
        public void TestIndexSorting()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var isotopes = new[] {0.8, 0.9, 0.6, 0.3};
            var index = Enumerable.Range(0, isotopes.Length).ToArray();

            Array.Sort(index, (i,j) => isotopes[j].CompareTo(isotopes[i]));

            for (var i = 0; i < isotopes.Length; i++)
            {
                Console.WriteLine("{0}\t{1}", index[i], isotopes[index[i]]);
            }
        }

        [Test]
        public void TestEnum()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const ActivationMethod activation = ActivationMethod.ETD;
            const byte code = (byte) activation;
            Console.WriteLine((byte)activation);
            Console.WriteLine(code);
            Console.WriteLine((ActivationMethod)code);
        }

        [Test]
        public void TestOverflow()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            Console.WriteLine("{0}", Math.Exp(13021));
        }

        [Test]
        public void TestRegEx()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            //const string protAnnotation = "A.HAHLTHQYPAANAQVTAAPQAITLNFSEGVETGFSGAKITGPKNENIKTLPAKRNEQDQKQLIVPLADSLKPGTYTVDWHVVSVDGHKTKGHYTFSVK._";
            //const char delimiter = (char)FastaDatabase.Delimiter;
            //Console.WriteLine(@"^[A-Z" + delimiter + @"]\.[A-Z]+\.[A-Z" + delimiter + @"_]$");
            //Console.WriteLine(Regex.IsMatch(protAnnotation, @"^[A-Z" + delimiter + @"]\.[A-Z]+\.[A-Z" + delimiter + @"_]$"));
            //const string s = "H(12)";
            //Console.WriteLine(Regex.IsMatch(s, @"^\d*[a-zA-Z]+$"));
            //const string plaincompositionStr = "H(12) C(4) 13C(3) N 15N O asdf1";
            const string plaincompositionStr = "H(12)";
            Console.WriteLine(Regex.IsMatch(plaincompositionStr, @"^([A-Z][a-z]?\d*)+$"));
        }

        [Test]
        public void TestParsingComposition()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string compositionStr = "H(230) C(136) N(40) O(46) S 13C(6) 15N(2)";
            var composition = Composition.Parse(compositionStr);
            Console.WriteLine("{0}\t{1}", composition, composition.Mass);
        }

        [Test]
        public void TestIonTypeFactory()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string sequenceStr = "PEPTIDE";
            var aminoAcidSet = new AminoAcidSet();
            var sequence = new Sequence(sequenceStr, aminoAcidSet);
            var compositionOfFirstPrefix = sequence.GetComposition(0, 2);
            Console.WriteLine("{0}\t{1}", compositionOfFirstPrefix, compositionOfFirstPrefix.Mass);

            var ionTypeFactory = new IonTypeFactory(new[] { BaseIonType.B, BaseIonType.Y }, new[] { NeutralLoss.NoLoss }, 10);
            var bIon = ionTypeFactory.GetIonType("b");
            var yIon = ionTypeFactory.GetIonType("y2");
            Console.WriteLine("{0}\t{1}\t{2}", bIon, bIon.OffsetComposition, bIon.OffsetComposition.Mass);
            Console.WriteLine("{0}\t{1}\t{2}", yIon, yIon.OffsetComposition, yIon.OffsetComposition.Mass);

            // Compute mass of y2 = AveragineMass(DE) + OffsetY
            var compositionOfSecondSuffix = sequence.GetComposition(sequenceStr.Length - 2, sequenceStr.Length);
            var y2Ion = yIon.GetIon(compositionOfSecondSuffix);
            Console.WriteLine("m/z of y++: {0}", y2Ion.GetMonoIsotopicMz());
        }

        [Test]
        public void TestDeconvolutedIonTypes()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var comp1 = Composition.H2O;
            var comp2 = new CompositionWithDeltaMass(-1);
            Console.WriteLine(comp1 + comp2);
            Console.WriteLine(comp1 - comp2);

            var ionTypeFactory = IonTypeFactory.GetDeconvolutedIonTypeFactory(new[] {BaseIonType.B, BaseIonType.Y}, new[] { NeutralLoss.NoLoss, NeutralLoss.H2O});
            foreach (var ionType in ionTypeFactory.GetAllKnownIonTypes())
            {
                Console.WriteLine(ionType);
            }
        }

        [Test]
        public void TestLinq()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            var intArr = new List<double> {4.0, 1.0, 2.0, 6.0};
            Console.WriteLine(intArr.Median());
        }

        [Test]
        public void TestRegularExpressions()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string str = "C2HBr365Ag2";
            var matches = Regex.Matches(str, @"[A-Z][a-z]?\d*");
            foreach (var match in matches)
            {
                var element = match.ToString();
                var atom = Regex.Match(element, @"[A-Z][a-z]?");
                var num = element.Substring(atom.Index + atom.Length);
                if (num.Length == 0) num = "1";
                Console.WriteLine("{0} ({1})", atom, num);
            }

            var comp = Composition.ParseFromPlainString(str);
            Console.WriteLine(comp.ToPlainString());
            Console.WriteLine(comp.ToString());
        }

        [Test]
        [Category("Local_Testing")]
        public void TestNumIsoWindows()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            //const string specFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";   // DDA
            const string specFilePath = @"\\protoapps\UserData\Wilkins\BottomUp\DIA_10mz\data\Q_2014_0523_50_10_fmol_uL_10mz.raw"; // DIA

            if (!File.Exists(specFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, specFilePath);
            }

            var run = PbfLcMsRun.GetLcMsRun(specFilePath);
            Console.WriteLine("NumIsoWindows: " + run.GetNumUniqueIsoWindows());
            Console.WriteLine("MinWidth: " + run.GetMinIsolationWindowWidth());
        }

        [Test]
        [TestCase("EAQADAAAEIAEDAAEAEDAGKPK", @"\\protoapps\UserData\Wilkins\BottomUp\DIA_10mz\data\Q_2014_0523_50_10_fmol_uL_10mz.raw", 54407)]
        [TestCase("KYETIDSLQIDDLMNRREVRQPADWQADENGSNDKGNGKGEPAVKVDEVVKSAPAEAELKDADESPVK",
                  @"\\protoapps\UserData\Wilkins\TopDown\Anil\QC_Shew_IntactProtein_new_CID-30CE-4Sep14_Bane_C2Column_3.raw", 2755)]
        [Category("Local_Testing")]
        public void TestPpmErrorCalculation(string seqText, string rawFilePath, int scanNum)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName, rawFilePath);

            var tolerance = new Tolerance(10, ToleranceUnit.Ppm);
            const int maxCharge = 15;
            const double relIntThres = 0.1;

            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            // init
            var sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(seqText);
            var lcms = PbfLcMsRun.GetLcMsRun(rawFilePath);
            var spectrum = lcms.GetSpectrum(scanNum);

            var ionTypeFactory = new IonTypeFactory(maxCharge);
            var iontypes = ionTypeFactory.GetAllKnownIonTypes();

            foreach (var iontype in iontypes)
            {
                var ion = iontype.GetIon(sequence.Composition);
                var obsPeaks = spectrum.GetAllIsotopePeaks(ion, tolerance, relIntThres);
                if (obsPeaks == null) continue;
                var isotopes = ion.GetIsotopes(relIntThres).ToArray();
                for (var i = 0; i < isotopes.Length; i++)
                {
                    if (obsPeaks[i] == null) continue;
                    var obsMz = obsPeaks[i].Mz;
                    var theoMz = ion.GetIsotopeMz(isotopes[i].Index);
                    var ppmError = (obsMz - theoMz)/theoMz*1e6;
                    Assert.True(ppmError <= tolerance.GetValue());
                }
            }
        }

        [Test]
        public void TestUpdateModificationComposition()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string modName1 = "mod1";
            var composition1 = new Composition(0, 1, 0, 2, 0, 0);
            var mod = Modification.RegisterAndGetModification(modName1, composition1);
            Assert.NotNull(mod);

            var getmod = Modification.Get(modName1);
            Assert.NotNull(getmod);
            Assert.AreEqual(mod, getmod);

            const string modName2 = "mod2";
            var composition2 = new Composition(1,0,2,4,1,9);
            var editedMod = Modification.UpdateAndGetModification(modName2, composition2);
            Assert.NotNull(editedMod);

            var getEditedMod = Modification.Get(modName2);
            Assert.NotNull(getEditedMod);
            Assert.AreEqual(editedMod, getEditedMod);
            Assert.AreNotEqual(mod, getEditedMod);
        }

        [Test]
        public void TestUpdateModificationMass()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string modName1 = "mod1";
            const int mass1 = 1000;
            var mod = Modification.RegisterAndGetModification(modName1, mass1);
            Assert.NotNull(mod);

            var getmod = Modification.Get(modName1);
            Assert.NotNull(getmod);
            Assert.AreEqual(mod, getmod);

            const string modName2 = "mod2";
            const int mass2 = 1100;
            var editedMod = Modification.UpdateAndGetModification(modName2, mass2);
            Assert.NotNull(editedMod);

            var getEditedMod = Modification.Get(modName2);
            Assert.NotNull(getEditedMod);
            Assert.AreEqual(editedMod, getEditedMod);
            Assert.AreNotEqual(mod, getEditedMod);
        }

        private void AppendModCombo(
            IDictionary<List<Modification>, Tuple<int, double, double, int, int>> modCombos,
            List<Modification> modList,
            Tuple<int, double, double, int, int> expectedResults)
        {
            modCombos.Add(modList, expectedResults);
        }


        private void ShowIsotopes(IEnumerable<Isotope> isotopes)
        {
            ShowIsotopes(isotopes.ToList());
        }

        private void ShowIsotopes(List<Isotope> isotopes)
        {
            var maxIntensity = (from item in isotopes select item.Ratio).Max();
            if (Math.Abs(maxIntensity) < float.Epsilon)
                maxIntensity = 1.0;

            foreach (var isotope in isotopes)
            {
                var plotBar = new string('+', (int)(isotope.Ratio / maxIntensity * 100 / 5));
                Console.WriteLine("{0,2} {1,8:N5} {2}", isotope.Index, isotope.Ratio, plotBar);
            }

        }

        private void ShowProductIons(IonType ionType, int indexStart, IEnumerable<Ion> productIons)
        {
            var i = 0;
            foreach (var ion in productIons)
            {
                Console.WriteLine("{0,-10} {1,2}  {2,-30} {3,3} {4:F6}", ionType.Name, indexStart + i, ion.Composition, ion.Charge, ion.GetMonoIsotopicMz());
                i++;
            }
        }


    }
}
