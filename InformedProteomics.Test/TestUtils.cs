﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml.Linq;
using DeconTools.Backend.Utilities.IsotopeDistributionCalculation;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Statistics;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestUtils
    {
        [Test]
        public void ParseMsGfString()
        {
            const string msgfPepStr = "+229.163C+57.021GLGGSGTPVDELDK+229.163C+57.021C+57.021QTHDNC+57.021YDQAK+229.163";
            var sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(msgfPepStr);
            var ion = new Ion(sequence.Composition + Composition.H2O, 4);
            Console.WriteLine("{0} {1} {2}", msgfPepStr, ion.Composition, ion.GetMonoIsotopicMz());

            var comp = Composition.Parse("H(257) C(150) N(42) O(56) S(4) 13C(12) 15N(3)");
            Console.WriteLine("AveragineMass: " + comp.Mass);
        }

        [Test]
        public void TestFormatting()
        {
            Console.WriteLine("{0:N3}", Modification.Tmt6Plex.Composition.Mass);
            Console.WriteLine(Modification.GetFromMass("229.163")[0].Name);
        }

        [Test]
        public void TestAveragine()
        {
            const double monoMass = 10247.5293287335;
            const int charge = 14;
            var profile = Averagine.GetTheoreticalIsotopeProfile(monoMass, charge);

            Console.WriteLine("Isotope ions:");
            foreach (var p in profile) Console.WriteLine("{0}\t{1}", p.Mz, p.Intensity);
            Console.WriteLine();
        }

        [Test]
        public void TestAaSet()
        {
            var oxEverywhere = new SearchModification(Modification.Oxidation, '*', SequenceLocation.Everywhere, true);
            var acetylNTerm = new SearchModification(Modification.Acetylation, '*', SequenceLocation.PeptideNTerm,
                false);
            var iTraqCTerm = new SearchModification(Modification.Itraq4Plex, '*', SequenceLocation.PeptideCTerm, true);
            var iTraqK = new SearchModification(Modification.Itraq4Plex, 'K', SequenceLocation.PeptideCTerm, true);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var carbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C',
                SequenceLocation.Everywhere, true);
            var dehydroC = new SearchModification(Modification.Dehydro, 'C', SequenceLocation.Everywhere, false);
            var glutathioneC = new SearchModification(Modification.Glutathione, 'C', SequenceLocation.Everywhere, false);
            var nitrosylC = new SearchModification(Modification.Nitrosyl, 'C', SequenceLocation.Everywhere, false);
            var nethylmaleimideC = new SearchModification(Modification.Nethylmaleimide, 'C', SequenceLocation.Everywhere, false);
            const int numMaxModsPerProtein = 4;
            var searchModifications = new List<SearchModification>
            {
                oxEverywhere
                //iTraqK,
                //iTraqCTerm
                //acetylNTerm
                //carbamidomethylC,
                //dehydroC,
                //glutathioneC,
                //nitrosylC,
                //nethylmaleimideC,
                //oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPerProtein);
            aaSet.Display();
        }

        [Test]
        public void TestPeptide()
        {
            //const string sequence = "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG";  // Histone H4
            const string sequence = "IRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG";  // Histone H4
            //const string sequence = "MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVTDELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGIVVDYVLEFDVPDELIVDRIVGRRVHAASGRVYHVKFNPPKVEGKDDVTGEDLTTRKDDQEETVRKRLVEYHQMTAPLIGYYQKEAEAGNTKYAKVDGTQAVADVRAALEKILG";
            //const string sequence = "MNKTQLIDVIAEKAELSKTQAKAALESTLAAITESLKEGDAVQLVGFGTFKVNHRAERTGRNPQTGKEIKIAAANVPAFVSGKALKDAVK";
            //const string sequence =
            //    "METTKPSFQDVLEFVRLFRRKNKLQREIQDVEKKIRDNQKRVLLLDNLSDYIKPGMSVEAIQGIIASMKGDYEDRVDDYIIKNAELSKERRDISKKLKAMGEMKNGEAK";
            var aaSet = new AminoAcidSet();
            var composition = aaSet.GetComposition(sequence) + Composition.H2O;

            Console.WriteLine(composition);
            Console.WriteLine(composition.Mass);
            Console.WriteLine(composition.NominalMass);
            // 2nd isotope
            Console.WriteLine(composition.GetIsotopeMass(0));
            Console.WriteLine(composition.GetIsotopeMass(1));
            Console.WriteLine(composition.GetIsotopeMass(2));
            //Assert.AreEqual(composition.ToPlainString(), "C34H51N7O14");

            Console.WriteLine("Isotopomer Envelope:");
            foreach (var e in composition.GetIsotopomerEnvelopeRelativeIntensities()) Console.WriteLine(e);
            Console.WriteLine();

            Console.WriteLine("Isotope ions:");
            var ion = new Ion(composition + Composition.H2O, 13);
            foreach (var p in ion.GetIsotopes(0.1)) Console.WriteLine("{0}\t{1}", ion.GetIsotopeMz(p.Index), p.Ratio);
            Console.WriteLine();
        }

        [Test]
        public void TestIsoProfile()
        {
            var composition = Composition.Parse("C(82) H(149) N(23) O(24) S(3)");
            const int charge = 3;
            var ion = new Ion(composition, charge);
            foreach (var isotope in ion.GetIsotopes(0.1))
            {
                Console.WriteLine("{0}\t{1}\t{2}", isotope.Index, ion.GetIsotopeMz(isotope.Index), isotope.Ratio);
            }
            Console.WriteLine(composition.Mass);
        }

        [Test]
        public void TestGetProductIons()
        {
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var sequence = new Sequence("CCAADDKEACFAVEGPK", aaSet);

            var ionTypeFactory = new IonTypeFactory(
                new[] {BaseIonType.B, BaseIonType.Y}, 
                new[] {NeutralLoss.NoLoss, NeutralLoss.H2O}, 
                maxCharge: 2);

            Console.WriteLine("Precursor Ion: {0}\t{1}", sequence.GetPrecursorIon(2).Composition, sequence.GetPrecursorIon(2).GetMonoIsotopicMz());
            Console.WriteLine("Product ions: ");
            var productIons = sequence.GetProductIons(ionTypeFactory.GetAllKnownIonTypes());
            foreach (var theoIon in productIons)
            {
                var ionTypeAndIndex = theoIon.Key;
                var ionType = ionTypeAndIndex.Item1;
                var index = ionTypeAndIndex.Item2;
                var ion = theoIon.Value;
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", ionType.Name, index, ion.Composition, ion.Charge, ion.GetMonoIsotopicMz());
            }
        }

        [Test]
        public void TestSequenceGraph()
        { 
            var phosPhoS = new SearchModification(Modification.Phosphorylation, 'S', SequenceLocation.Everywhere, false);
            var phosPhoT = new SearchModification(Modification.Phosphorylation, 'T', SequenceLocation.Everywhere, false);
            var phosPhoY = new SearchModification(Modification.Phosphorylation, 'Y', SequenceLocation.Everywhere, false);
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var fixCarbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C', SequenceLocation.Everywhere, true);

            var searchModifications = new List<SearchModification> { phosPhoS, phosPhoT, phosPhoY, oxM, fixCarbamidomethylC };
            //var searchModifications = new List<SearchModification> { phosPhoT, fixCarbamidomethylC };
            const int numMaxModsPepPeptide = 2;

            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPepPeptide);
            const string annotation = "_.STR._";
            var pepSeq = annotation.Substring(2, annotation.Length - 4);
            Console.WriteLine(aaSet.GetComposition(pepSeq));
            var graph = SequenceGraph.CreateGraph(aaSet, annotation);
            Console.WriteLine(graph.GetUnmodifiedSequenceComposition());
            Assert.AreEqual(graph.GetUnmodifiedSequenceComposition(), aaSet.GetComposition(pepSeq));

            Console.WriteLine("Annotation Compositions:");
            var index = -1;
            foreach (var composition in graph.GetSequenceCompositions())
            {
                Console.WriteLine(++index+": "+composition);
            }

            //const int seqIndex = 1;
            //Console.WriteLine("Fragment Compositions (" + seqIndex +")");
            //var scoringGraph = graph.GetScoringGraph(seqIndex);
            //foreach (var composition in scoringGraph.GetCompositions())
            //{
            //    Console.WriteLine(composition);
            //}
        }



        [Test]
        public void TestModificationParams()
        {
            var modifications = new[] { Modification.Acetylation, Modification.Phosphorylation, Modification.Oxidation}; //, Modification.PyroGluQ };
            var modParams = new ModificationParams(modifications, 3);
            int numCombinations = modParams.NumModificationCombinations;
            for (int modCombIndex = 0; modCombIndex < numCombinations; modCombIndex++)
            {
                var modCombination = modParams.GetModificationCombination(modCombIndex);
                Console.WriteLine("{0}: {1} {2} {3}", modCombIndex, modCombination, modCombination.Composition, modCombination.Composition.Mass);
            }

            Console.WriteLine(modParams.GetModificationCombinationIndex(8, 0));
            Console.WriteLine(modParams.GetModificationCombinationIndex(19, 1));
        }

        [Test]
        public void TestGeneratingCombinations()
        {
            const int n = 5;
            const int k = 3;
            var combinations = SimpleMath.GetCombinationsWithRepetition(n, k);
            int count = 0;
            foreach (var combination in combinations)
            {
                ++count;
                Console.WriteLine((count-1)+": "+string.Join(",", combination));
            }
            Console.WriteLine("Count: " + count);
            Console.WriteLine("Count2: " + SimpleMath.GetCombination(n + k - 1, k));
        }

        [Test]
        public void TestCSharpSyntax()
        {
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
        public void TestIonTypeGeneration()
        {
            var ionTypeFactory = new IonTypeFactory();
            int index = 0;
            foreach (var ionType in ionTypeFactory.GetAllKnownIonTypes())
            {
                Console.WriteLine(++index + ": " + ionType);
            }
            var yIon = ionTypeFactory.GetIonType("y2-H2O");
            Console.WriteLine(yIon.GetMz(0));
        }

        [Test]
        public void TestIsotopomerProfile()
        {
            //const string molFormula = "C78H120N22O28S3";    // CCAADDKEACFAVEGPK
            const string molFormula = "C83H136N22O24S1"; 
//            const string molFormula = "C4195H6470N1164O1213S34";

            var isoCalc = IsotopicDistributionCalculator.Instance;
            var profile = isoCalc.GetIsotopePattern(molFormula);

            var sb = new StringBuilder();

            foreach (var peak in profile.Peaklist)
            {
                sb.Append(peak.XValue);
                sb.Append("\t");
                sb.Append(peak.Height);
                sb.Append("\n");
            }

            Console.Write(sb.ToString());
        }

        [Test]
        public void TestIsotopemerProfileByKyowon() // is faster and more accurate than IsotopicDistributionCalculator
        {
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
                Console.WriteLine("{0}: {1}\t{2}", isotopeIndex, ion.GetIsotopeMz(isotopeIndex), ii);
            }
        }

        [Test]
        public void TestTimeToComputeIsotopomerProfiles()
        {
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\HumanPeptides.txt";

            int numPeptides = 0;
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            var isoCalc = IsotopicDistributionCalculator.Instance;
            foreach (var annotation in File.ReadLines(dbFilePath))
            {
                ++numPeptides;
                var peptide = annotation.Substring(2, annotation.Length - 4);
                var composition = aaSet.GetComposition(peptide);
                var molFormula = composition.ToPlainString();
                isoCalc.GetIsotopePattern(molFormula);
                composition.GetIsotopomerEnvelopeRelativeIntensities();
            }

            Console.WriteLine("NumPeptides: " + numPeptides);
            sw.Stop();
            var sec = sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }

        [Test]
        public void TestReadingPnnlOmicsXmlFile()
        {
            const string xmlFileName = @"..\..\..\PNNLOmicsElementData.xml";
            var xdocument = XDocument.Load(xmlFileName);
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
                    Console.WriteLine("{0} {1}", symbol.Value, name.Value);
                }
            }
        }

        [Test]
        public void TestIndexSorting()
        {
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
            const ActivationMethod activation = ActivationMethod.ETD;
            const byte code = (byte) activation;
            Console.WriteLine((byte)activation);
            Console.WriteLine(code);
            Console.WriteLine((ActivationMethod)code);
        }

        [Test]
        public void TestOverflow()
        {
            Console.WriteLine("{0}", Math.Exp(13021));
        }

        [Test]
        public void TestRegEx()
        {
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
            const string compositionStr = "H(230) C(136) N(40) O(46) S 13C(6) 15N(2)";
            var composition = Composition.Parse(compositionStr);
            Console.WriteLine("{0}\t{1}", composition, composition.Mass);
        }

        [Test]
        public void TestIonTypeFactory()
        {
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
            var intArr = new List<double> {4.0, 1.0, 2.0, 6.0};
            Console.WriteLine(intArr.Median());
        }

        [Test]
        public void TestRegularExpressions()
        {
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
        public void TestNumIsoWindows()
        {
            //const string specFilePath = @"C:\cygwin\home\kims336\Data\QCShewQE\QC_Shew_13_04_A_17Feb14_Samwise_13-07-28.raw";   // DDA
            const string specFilePath = @"\\protoapps\UserData\Wilkins\BottomUp\DIA_10mz\data\Q_2014_0523_50_10_fmol_uL_10mz.raw"; // DIA
            var run = PbfLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);
            Console.WriteLine("NumIsoWindows: " + run.GetNumUniqueIsoWindows());
            Console.WriteLine("MinWidth: " + run.GetMinIsolationWindowWidth());
        }

        [TestCase("EAQADAAAEIAEDAAEAEDAGKPK", @"\\protoapps\UserData\Wilkins\BottomUp\DIA_10mz\data\Q_2014_0523_50_10_fmol_uL_10mz.raw", 54407)]
        [TestCase("KYETIDSLQIDDLMNRREVRQPADWQADENGSNDKGNGKGEPAVKVDEVVKSAPAEAELKDADESPVK",
                  @"\\protoapps\UserData\Wilkins\TopDown\Anil\QC_Shew_IntactProtein_new_CID-30CE-4Sep14_Bane_C2Column_3.raw", 2755)]
        public void TestPpmErrorCalculation(string seqText, string rawFilePath, int scanNum)
        {
            var tolerance = new Tolerance(10, ToleranceUnit.Ppm);
            const int maxCharge = 15;
            const double relIntThres = 0.1;

            // init
            var sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(seqText);
            var lcms = PbfLcMsRun.GetLcMsRun(rawFilePath, MassSpecDataType.XCaliburRun);
            var spectrum = lcms.GetSpectrum(scanNum);

            var ionTypeFactory = new IonTypeFactory(maxCharge);
            var iontypes = ionTypeFactory.GetAllKnownIonTypes();

            foreach (var iontype in iontypes)
            {
                var ion = iontype.GetIon(sequence.Composition);
                var obsPeaks = spectrum.GetAllIsotopePeaks(ion, tolerance, relIntThres);
                if (obsPeaks == null) continue;
                var isotopes = ion.GetIsotopes(relIntThres).ToArray();
                for (int i = 0; i < isotopes.Length; i++)
                {
                    if (obsPeaks[i] == null) continue;
                    var obsMz = obsPeaks[i].Mz;
                    var theoMz = ion.GetIsotopeMz(isotopes[i].Index);
                    var ppmError = (obsMz - theoMz)/theoMz*1e6;
                    Assert.True(ppmError <= tolerance.GetValue());
                }
            }
        }
    }
}
