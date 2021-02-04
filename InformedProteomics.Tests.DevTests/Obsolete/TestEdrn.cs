using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.MathAndStats;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.DevTests.Obsolete
{
    [TestFixture]
    public class TestEdrn
    {
        [Ignore("File Missing, test obsolete, or long test")]
        [Test]
        [Category("Local_Testing")]
        public void GenerateVennDiagrams()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

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
        [Category("Local_Testing")]
        public void SummarizeDda()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

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
                        if (line.StartsWith("Peptide"))
                        {
                            continue;
                        }

                        var token = line.Split('\t');
                        var peptide = token[0];

                        var fraction = Convert.ToInt32(resultFile.Substring(resultFile.IndexOf("_Serum", System.StringComparison.Ordinal) + 7, 2));
                        if (peptide.Equals(spikedInPeptide))
                        {
                            fractionList.Add(fraction);
                        }
                    }
                }
                if (fractionList.Count > 0)
                {
                    numIdentifiedPeptides++;
                }

                Console.WriteLine("{0}\t{1}", spikedInPeptide, string.Join(",", fractionList));
            }
            Console.WriteLine("NumPeptides: {0}", numPeptides);
            Console.WriteLine("NumIdentifiedPeptides: {0}", numIdentifiedPeptides);
        }

        [Ignore("File Missing, test obsolete, or long test")]
        [Test]
        [Category("Local_Testing")]
        public void ComputeSpikedInPeptideMzHist()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string pepListFile = @"C:\cygwin\home\kims336\Data\DIA\SpikedPeptides.txt";

            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var charges = new[] {2};

            var hist = new int[4];

            var sum = 0;

            Console.WriteLine("Peptide\tCharge\tMz");
            foreach (var line in File.ReadLines(pepListFile))
            {
                if (line.Length == 0)
                {
                    continue;
                }

                var peptide = line;
                var composition = aaSet.GetComposition(peptide) + Composition.H2O;

                foreach (var charge in charges)
                {
                    var precursorIon = new Ion(composition, charge);
                    var precursorIonMz = precursorIon.GetMonoIsotopicMz();

                    if (precursorIonMz < 400 || precursorIonMz >= 900)
                    {
                        continue;
                    }

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
