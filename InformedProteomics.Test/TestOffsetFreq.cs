using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    public class IonProbability
    {
        public int Found { get; set; }
        public int Total { get; set; }
        public IonProbability(int f, int t) { Found = f; Total = t; }
    }

    [TestFixture]
    public class OffsetTable
    {
        private const double FdrThreshold = 0.01;
        private const double EValueThreshold = 0.1;
        private const double PepQThreshold = 0.01;
        private const int PrecursorCharge = 2;
        private const int TotalCharges = 2;
        private const int NumMutations = 3;
        private const string PrecChargeHeader = "Charge";
        private const string PeptideHeader = "Peptide";
        private const string ScanHeader = "ScanNum";
        private const string EvalueHeader = "EValue";
        private const string FDRHeader = "FDR";
        private static IEnumerable<Tuple<string, Spectrum>> CleanScans(string txtFileName, string rawFileName)
        {
            var tsvParser = new TsvFileParser(txtFileName);

            var scans = tsvParser.GetData(ScanHeader);
            var peptides = tsvParser.GetPeptides(PepQThreshold);
            var fdrs = tsvParser.GetData(FDRHeader);
            var evalues = tsvParser.GetData(EvalueHeader);
            var charges = tsvParser.GetData(PrecChargeHeader);
            IEnumerable<string> minValues = fdrs;
            var lcms = LcMsRun.GetLcMsRun(rawFileName, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

            var clean = new List<Tuple<string, Spectrum>>();
            using (var chargei = charges.GetEnumerator())
            using (var scani = scans.GetEnumerator())
            using (var peptidei = peptides.GetEnumerator())
            {
                while (scani.MoveNext() && peptidei.MoveNext() && chargei.MoveNext())
                {
                    var spec = lcms.GetSpectrum(Convert.ToInt32(scani.Current));
                    int precCharge = Convert.ToInt32(chargei.Current);
                    if (precCharge != PrecursorCharge) continue;
                    clean.Add(new Tuple<string, Spectrum>(peptidei.Current, spec));
                }
            }
            return clean;
        }

        private static Dictionary<string, IonProbability> GetOffsetCounts(IEnumerable<Tuple<string, Spectrum>> cleanScans,
                                        string[] ionTypes, bool useDecoy, ActivationMethod act, IonTypeFactory ionTypeFactory)
        {
            var probabilities = new Dictionary<string, IonProbability>();
            foreach (string ionTypeStr in ionTypes)
            {
                if (!probabilities.ContainsKey(ionTypeStr))
                    probabilities.Add(ionTypeStr, new IonProbability(0, 0));
            }

            var pepDict = new Dictionary<string, int>();

            foreach (var node in cleanScans)
            {
                var protein = node.Item1;
                var spectrum = node.Item2;
                if (pepDict.ContainsKey(protein)) continue;
                else
                    pepDict.Add(protein, 0);
                if (useDecoy)
                {
                    var shuffled = SimpleStringProcessing.Shuffle(protein);
                    protein = SimpleStringProcessing.Mutate(shuffled, NumMutations);
                }
                var sequence = Sequence.GetSequenceFromMsGfPlusPeptideStr(protein);
                var spec = spectrum as ProductSpectrum;
                if (spec == null) continue;
                for (int i = 0; i < protein.Length - 1; i++)
                {
                    if (spec.ActivationMethod == act)
                    {
                        foreach (var ionTypeStr in ionTypes)
                        {
                            Composition sequenceComposition = null;
                            if (ionTypeStr[0] == 'a' || ionTypeStr[0] == 'b' || ionTypeStr[0] == 'c')
                                sequenceComposition = sequence.GetComposition(0, i);
                            else if (ionTypeStr[0] == 'x' || ionTypeStr[0] == 'y' || ionTypeStr[0] == 'z')
                                sequenceComposition = sequence.GetComposition(protein.Length - i, protein.Length);
                            else
                                throw new FormatException();

                            var ionType = ionTypeFactory.GetIonType(ionTypeStr);
                            var ion = ionType.GetIon(sequenceComposition);
                            ion.Composition.ComputeApproximateIsotopomerEnvelop();

                            probabilities[ionTypeStr].Total++;
                            if (spectrum.ContainsIon(ion, new Tolerance(15, ToleranceUnit.Ppm), 0.8))
                                probabilities[ionTypeStr].Found++;
                        }
                    }
                }
            }
            return probabilities;
        }

        private static void WriteProbFile(string outFile, Dictionary<string, IonProbability> offsetCounts, string[] ionTypes)
        {
            using (var file = new StreamWriter(outFile))
            {
                for (int i = 0; i < ionTypes.Length; i++)
                {
                    file.WriteLine("{0}\t{1}", ionTypes[i], Math.Round((double)(offsetCounts[ionTypes[i]].Found) / (offsetCounts[ionTypes[i]].Total), 5));
                }
            }
        }

        private static void WriteOffsetCountsFile(string outFile, Dictionary<string, IonProbability> offsetCounts, string[] ionTypes)
        {
            using (var file = new StreamWriter(outFile))
            {
                for (int i = 0; i < ionTypes.Length; i++)
                {
                    file.Write("{0}Found\t{0}Total", ionTypes[i]);
                    if (i != ionTypes.Length - 1)
                        file.Write("\t");
                }
                file.WriteLine();
                for (int i = 0; i < ionTypes.Length; i++)
                {
                    file.Write("{0}\t{1}", offsetCounts[ionTypes[i]].Found, offsetCounts[ionTypes[i]].Total);
                    if (i != ionTypes.Length - 1)
                        file.Write("\t");
                }
            }
        }

        [Test]
        public static void OffsetFreq()
        {
            const string fileList = @"C:\Users\wilk011\Documents\DataFiles\HCD_QE_TMT10.txt";
            const string preRes = @"\\protoapps\UserData\Sangtae\ForChris\HCD_QE_TMT10_Charles\tsv\";
            const string preRaw = @"\\protoapps\UserData\Sangtae\ForChris\HCD_QE_TMT10_Charles\raw\";
            string outPre = @"C:\Users\wilk011\Documents\DataFiles\BottomUp Offset\HCD_QE_TMT10_Charles\Charge "+PrecursorCharge+@"\";
            const string outSuff = ".out";
            const bool useDecoy = false;
            string finalOutput = "HCD_QE_TMT10_Charles.Charge" + PrecursorCharge + ".txt";
            string[] ionTypes = { "a", "a-H2O", "a-NH3", 
                                  "a2", "a2-H2O", "a2-NH3", 
                                  "b", "b-H2O", "b-NH3",
                                  "b2", "b2-H2O", "b2-NH3",
                                  "c", "c-H2O",
                                  "c2", "c2-H2O",
                                  "x", "x-H2O", "x-NH3",
                                  "x2", "x2-H2O", "x2-NH3",
                                  "y", "y-H2O", "y-NH3",
                                  "y2", "y2-H2O", "y2-NH3",
                                  "z", "z-H2O", "z-NH3",
                                  "z2", "z2-H2O", "z2-NH3" };
            var losses = new[] { NeutralLoss.NoLoss, NeutralLoss.NH3, NeutralLoss.H2O };
            const ActivationMethod act = ActivationMethod.HCD;

            var ionTypeFactory = new IonTypeFactory(new[] { BaseIonType.A, BaseIonType.B, BaseIonType.C, 
                                                            BaseIonType.X, BaseIonType.Y, BaseIonType.Z },
                                                    losses, TotalCharges);

            var fileNameParser = new TsvFileParser(fileList);

            var txtFiles = fileNameParser.GetData("text");
            var rawFiles = fileNameParser.GetData("raw");

            using (var txtFileIt = txtFiles.GetEnumerator())
            using (var rawFileIt = rawFiles.GetEnumerator())
            {
                while (txtFileIt.MoveNext() && rawFileIt.MoveNext())
                {
                    string textFile = preRes + txtFileIt.Current;
                    string rawFile = preRaw + rawFileIt.Current;
                    Console.WriteLine(rawFile);
                    var scans = CleanScans(textFile, rawFile);
                    var cleanScans = scans as Tuple<string, Spectrum>[] ?? scans.ToArray();
                    
                    var offsetCounts = GetOffsetCounts(cleanScans, ionTypes, false, act, ionTypeFactory);
                    var outFile = outPre + rawFileIt.Current + ".Charge" + PrecursorCharge + outSuff;
                    WriteOffsetCountsFile(outFile, offsetCounts, ionTypes);
                    WriteProbFile(outFile + ".prob", offsetCounts, ionTypes);
                    if (useDecoy)
                    {
                        var decoyOffsetCounts = GetOffsetCounts(cleanScans, ionTypes, true, act, ionTypeFactory);
                        var decoyOutFile = outPre + rawFileIt.Current + ".Charge" + PrecursorCharge + ".decoy" + outSuff;
                        WriteOffsetCountsFile(decoyOutFile, decoyOffsetCounts, ionTypes);
                        WriteProbFile(decoyOutFile + ".prob", decoyOffsetCounts, ionTypes);
                    }
                }
            }

            // Consolidate files
            var found = new int[ionTypes.Length];
            var total = new int[ionTypes.Length];
            foreach (var rawFile in rawFiles)
            {
                string fileName = outPre + rawFile + ".Charge" + PrecursorCharge + outSuff;
                var outReader = new TsvFileParser(fileName);
                for (int i = 0; i < ionTypes.Length; i++)
                {
                    var foundStr = outReader.GetData(ionTypes[i] + "Found")[0];
                    var totalStr = outReader.GetData(ionTypes[i] + "Total")[0];
                    found[i] = Convert.ToInt32(foundStr);
                    total[i] = Convert.ToInt32(totalStr);
                }
            }
            using (var finalOutputFile = new StreamWriter(outPre + finalOutput))
            {
                for (int i = 0; i < ionTypes.Length; i++)
                {
                    finalOutputFile.WriteLine("{0}\t{1}", ionTypes[i], Math.Round((double)(found[i]) / (total[i]), 5));
                }
            }
        }
    }
}
