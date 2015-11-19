using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.DIA.Search;
using NUnit.Framework;

namespace InformedProteomics.Test.Obsolete
{
    [TestFixture]
    public class TestEdrn
    {
        [Ignore("File Missing, test obsolete, or long test")]
        [Test]
        public void GenerateVennDiagrams()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            // DIA
            const string dir = @"H:\Research\EDRN\Ic\DIA_Replicate";

            const string rep1 = dir + @"\EDRN_Serum_07_DIA_1_01_13Nov13_Samwise_13-07-28_IcTda.tsv";
            //const string rep2 = dir + @"\EDRN_Serum_07_DIA_1_02_13Nov13_Samwise_13-07-28_IcTda.tsv";
            //const string rep3 = dir + @"\EDRN_Serum_07_DIA_1_03_13Nov13_Samwise_13-07-28_IcTda.tsv";
            const string rep4 = dir + @"\EDRN_Serum_07_DIA_1_04_13Nov13_Samwise_13-07-28_IcTda.tsv";
            //const string rep5 = dir + @"\EDRN_Serum_07_DIA_1_05_18Nov13_Samwise_13-07-28_IcTda.tsv";

            const string resultPath1 = rep1;
            const string resultPath2 = rep4;

            var result1 = new TsvFileParser(resultPath1);
            var result2 = new TsvFileParser(resultPath2);

            const double pepQValueThreshold = 0.01;
            var vennDiagram = new VennDiagram<string>(result1.GetPeptidesAboveQValueThreshold(pepQValueThreshold),
                                                      result2.GetPeptidesAboveQValueThreshold(pepQValueThreshold));
            Console.WriteLine("{0}\t{1}\t{2}",
                              vennDiagram.Set1Only.Count, // + vennDiagram.Intersection.Count,
                              vennDiagram.Intersection.Count,
                              vennDiagram.Set2Only.Count //+ vennDiagram.Intersection.Count
                              );
        }


        [Ignore("File Missing, test obsolete, or long test")]
        [Test]
        public void RunIpa()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string resultPath = @"D:\Research\Data\EDRN\Replicates_Frac7\343513_EDRN_Serum_07_DIA_1_05_18Nov13_Samwise_13-07-28.tsv";
            const string specFilePath = @"D:\Research\Data\EDRN\RawFiles\DIA_Replicate\343513_EDRN_Serum_07_DIA_1_05_18Nov13_Samwise_13-07-28.raw";
            const string outputFilePath = @"D:\Research\Data\EDRN\Replicates_Frac7\Rep5_Summary.tsv";

            var postProcessor = new MsGfPostProcessor(specFilePath, resultPath, new Tolerance(20), new Tolerance(10));
            var numId = postProcessor.PostProcessing(outputFilePath);

            Console.WriteLine("NumId: {0}", numId);
        }

        [Ignore("File Missing, test obsolete, or long test")]
        [Test]
        public void GenerateSpecCount()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string resultFileDir = @"D:\Research\Data\EDRN\DDA\NTT1_NoMod\";
            var resultFiles = Directory.GetFiles(resultFileDir, "*.tsv");

            var fractionToIds = new Dictionary<int, List<MsGfMatch>>();
            foreach (var resultFile in resultFiles)
            {
                var fraction = Convert.ToInt32(resultFile.Substring(resultFile.IndexOf("_Serum", StringComparison.Ordinal) + 7, 2));
                var msgfResults = new MsGfResults(resultFile);
                fractionToIds.Add(fraction, msgfResults.GetMatchesAtPsmFdr(0.01));
            }

            const string spikedInPeptideFile = @"D:\Research\Data\EDRN\SpikedPeptides.txt";
            var spikedInPeptides = File.ReadAllLines(spikedInPeptideFile);

            Console.WriteLine("Peptide\t{0}", string.Join("\t", Enumerable.Range(1, 12)));
            foreach (var spikedInPeptide in spikedInPeptides)
            {
                Console.Write(spikedInPeptide);
                for (var fraction = 1; fraction <= 12; fraction++)
                {
                    var matches = fractionToIds[fraction];
                    var count = matches.Count(m => m.Peptide.Replace("C+57.021", "C").Replace("K+8.014", "K").Replace("R+10.008", "R").Equals(spikedInPeptide));
                    Console.Write("\t{0}", count);
                }
                Console.WriteLine();
            }
        }

        [Ignore("File Missing, test obsolete, or long test")]
        [Test]
        public void SummarizeDda()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string spikedInPeptideFile = @"D:\Research\Data\EDRN\SpikedPeptides.txt";
            var spikedInPeptides = File.ReadAllLines(spikedInPeptideFile);

            //var peptideToFractions = new Dictionary<string, List<int>>();

            const string resultFileDir = @"D:\Research\Data\EDRN\DDA\Heavy\";

            var numPeptides = 0;
            var numIdentifiedPeptides = 0;
            foreach (var spikedInPeptide in spikedInPeptides)
            {
                numPeptides++;
                var fractionList = new List<int>();
                var resultFiles = Directory.GetFiles(resultFileDir, "*.prm.txt");
                foreach (var resultFile in resultFiles)
                {
                    foreach (var line in File.ReadLines(resultFile))
                    {
                        if (line.StartsWith("Peptide")) continue;
                        var token = line.Split('\t');
                        var peptide = token[0];

                        var fraction = Convert.ToInt32(resultFile.Substring(resultFile.IndexOf("_Serum", System.StringComparison.Ordinal) + 7, 2));
                        if (peptide.Equals(spikedInPeptide)) fractionList.Add(fraction);
                    }
                }
                if (fractionList.Any()) numIdentifiedPeptides++;
                Console.WriteLine("{0}\t{1}", spikedInPeptide, string.Join(",", fractionList));
            }
            Console.WriteLine("NumPeptides: {0}", numPeptides);
            Console.WriteLine("NumIdentifiedPeptides: {0}", numIdentifiedPeptides);
        }

        [Ignore("File Missing, test obsolete, or long test")]
        [Test]
        public void GeneratePrmInfo()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string resultFileDir = @"D:\Research\Data\EDRN\DDA\Heavy\";

            var resultFiles = Directory.GetFiles(resultFileDir, "*.tsv");
            foreach (var resultFile in resultFiles)
            {
                var outputPath = Path.ChangeExtension(resultFile, "prm.txt");
                //Console.WriteLine("{0} -> {1}", Path.GetFileName(resultFile), Path.GetFileName(outputPath));
                GeneratePrmInfo(resultFile, outputPath);
            }
        }

        [Ignore("File Missing, test obsolete, or long test")]
        [Test]
        public void GeneratePrmInfo(string resultFilePath, string outputFilePath)
        {
            Console.Write("Processing {0}", Path.GetFileName(resultFilePath));
            Console.Out.Flush();

            var rawFilePath =
                @"D:\Research\Data\EDRN\DDA\raw\" + Path.GetFileNameWithoutExtension(resultFilePath) + ".raw";
            var reader = new XCaliburReader(rawFilePath);
            var run = InMemoryLcMsRun.GetLcMsRun(rawFilePath);

            var tolerance = new Tolerance(10, ToleranceUnit.Ppm);

            const string spikedInPeptideFile = @"D:\Research\Data\EDRN\SpikedPeptides.txt";
            var spikedInPeptides = File.ReadAllLines(spikedInPeptideFile);
            var spikedInPepSet = new HashSet<string>();
            foreach (var p in spikedInPeptides)
            {
                spikedInPepSet.Add(p);
            }
//            const string resultFilePath = @"D:\Research\Data\EDRN\DDA\Frac7_NTT2.tsv";
            //const string resultFilePath = @"D:\Research\Data\EDRN\DDA\Heavy\342865_EDRN_Serum_07_DDA_1_12Nov13_Samwise_13-07-28.tsv";
//            const string resultFilePath = @"D:\Research\Data\EDRN\DDA\NTT1_NoMod\342865_EDRN_Serum_07_DDA_1_12Nov13_Samwise_13-07-28.tsv";
            const double qValueThreshold = 0.01;

            var pepSet = new HashSet<string>();
            MsGfPlusHeaderInformation headerInfo = null;

            //var prefix = new HashSet<string>();
            //var suffix = new HashSet<string>();
            var numPeptides = 0;

            var prevScanNum = -1;
            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine("Peptide\tCharge\tMonoMz\tMostAbundantMz\tMs2ScanNum\tRtMs2\tRtApex\tRtStart\tRtEnd\tSpecEValue\tPepQValue");
                foreach (var line in File.ReadLines(resultFilePath))
                {
                    if (line.StartsWith("#"))
                    {
                        headerInfo = new MsGfPlusHeaderInformation(line);
                        continue;
                    }

                    var match = new MsGfMatch(line, headerInfo);

                    if (match.ScanNum == prevScanNum) continue;
                    prevScanNum = match.ScanNum;

                    if (!match.IsValid || match.Protein.StartsWith(FastaDatabase.DecoyProteinPrefix)) continue;
                    if (match.PepQValue > qValueThreshold) continue;
                    var peptide = match.Peptide.Replace("C+57.021", "C").Replace("K+8.014", "K").Replace("R+10.008", "R");

                    if (pepSet.Contains(peptide)) continue;
                    pepSet.Add(peptide);

                    if (spikedInPepSet.Contains(peptide))
                    {
                        var ion = new Ion(match.Formula, match.Charge);
                        var mostAbundantIonMz = ion.GetMostAbundantIsotopeMz();
                        var xic = run.GetPrecursorExtractedIonChromatogram(mostAbundantIonMz, tolerance, match.ScanNum);
                        if (xic.Count == 0) continue;
                        var minScan = xic.Min().ScanNum;
                        var maxScan = xic.Max().ScanNum;
                        writer.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}",
                            peptide,
                            match.Charge,
                            ion.GetMonoIsotopicMz(),
                            mostAbundantIonMz,
                            match.ScanNum,
                            reader.RtFromScanNum(match.ScanNum),
                            reader.RtFromScanNum(xic.GetApexScanNum()),    // Rt apex
                            reader.RtFromScanNum(minScan),                 // Rt start
                            reader.RtFromScanNum(maxScan),                 // Rt end
                            match.SpecEValue,
                            match.PepQValue);
                        ++numPeptides;
                    }
                    //else
                    //{
                    //    foreach (var spikedInPeptide in spikedInPeptides)
                    //    {
                    //        if (spikedInPeptide.StartsWith(peptide)) prefix.Add(spikedInPeptide + "\t" + peptide + "\t" + match.ScanNum);
                    //        else if (spikedInPeptide.EndsWith(peptide)) suffix.Add(spikedInPeptide + "\t" + peptide + "\t" + match.ScanNum);
                    //    }
                    //}
                }                
            }

            //Console.WriteLine("*********Prefix");
            //foreach(var p in prefix) Console.WriteLine(p);

            //Console.WriteLine("*********Suffix");
            //foreach (var p in suffix) Console.WriteLine(p);

            Console.WriteLine("\t{0}", numPeptides);
        }

        [Ignore("File Missing, test obsolete, or long test")]
        [Test]
        public void ComputeSpikedInPeptideMzHist()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string pepListFile = @"C:\cygwin\home\kims336\Data\DIA\SpikedPeptides.txt";

            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var charges = new[] {2};

            var hist = new int[4];

            var sum = 0;

            Console.WriteLine("Peptide\tCharge\tMz");
            foreach (var line in File.ReadLines(pepListFile))
            {
                if (line.Length == 0) continue;
                var peptide = line;
                var composition = aaSet.GetComposition(peptide) + Composition.H2O;

                foreach (var charge in charges)
                {
                    var precursorIon = new Ion(composition, charge);
                    var precursorIonMz = precursorIon.GetMonoIsotopicMz();

                    if (precursorIonMz < 400 || precursorIonMz >= 900) continue;
                    var histIndex = (int)((precursorIonMz - 400)/125);
                    hist[histIndex]++;

                    Console.WriteLine("{0}\t{1}\t{2}\t{3}", peptide, charge, precursorIonMz, histIndex);

                    sum++;
                }
            }

            Console.WriteLine("\nRange\tNum\tRatio");
            for (var i = 0; i < hist.Length; i++)
            {
                Console.WriteLine("{0}-{1}\t{2}\t{3}", 400+i*125, 525+i*125, hist[i], hist[i] / (float)sum);
            }

        }
    }
}
