using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.DevTests
{
    [TestFixture]
    public class TestIprg
    {
        [Ignore("File Missing, test obsolete, or long test")]
        [Test]
        [Category("Local_Testing")]
        public void TestReadingExcelFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            //const string resultFile = @"H:\Research\IPRG2015\Submissions\Submission_32080.xlsx";
        }

        [Test]
        [Category("Local_Testing")]
        public void ProcessIprg2015PreStudy()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string dir = @"H:\Research\IPRG2015";

            const string databaseFilePath = dir + @"\database\yeast6proteaprotein.fasta";
            if (!File.Exists(databaseFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, databaseFilePath);
            }

            var database = new FastaDatabase(databaseFilePath);
            database.Read();

            const string jobFilePath = dir + @"\Jobs.tsv";
            if (!File.Exists(jobFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, jobFilePath);
            }

            var jobParser = new TsvFileParser(jobFilePath);
            var jobs = jobParser.GetData("Jobs").Select(j => Convert.ToInt32(j)).ToArray();
            var experiments = jobParser.GetData("Experiments").Select(e => e.Split('_')[2]).ToArray();

            //const string resultFilePath = dir + @"\AMT_Proteins_NA.tsv";
            //const string outputFilePath = dir + @"\AMT_Proteins.tsv";

            const string resultFilePath = dir + @"\AMT_Peptides_NA.tsv";
            const string outputFilePath = dir + @"\AMT_Peptides.tsv";

            var parser = new TsvFileParser(resultFilePath);
            var headers = parser.GetHeaders();
            var jobColNum = new int[jobs.Length];
            for (var i = 0; i < jobs.Length; i++)
            {
                for (var j = 0; j < headers.Count; j++)
                {
                    if (headers[j].Contains(jobs[i].ToString()))
                    {
                        jobColNum[i] = j;
                        break;
                    }
                }
            }

            for (var i = 0; i < jobs.Length; i++)
            {
                Console.WriteLine("{0}\t{1}\t{2}", jobs[i], jobColNum[i], experiments[i]);
            }

            using var writer = new StreamWriter(outputFilePath);

            var peptides = parser.GetData("Peptide");   // Peptides
            var proteins = parser.GetData("Reference"); // Proteins
            var abundances = new string[jobs.Length][];
            for (var i = 0; i < jobs.Length; i++)
            {
                abundances[i] = parser.GetData(headers[jobColNum[i]]).ToArray();
            }

            if (peptides != null)
            {
                writer.Write("Peptide\t");
            }

            writer.Write("Protein\tLength");
            for (var i = 0; i < jobs.Length; i++)
            {
                writer.Write("\t" + experiments[i]);
            }
            writer.WriteLine("\tSpikeIn");
            for (var i = 0; i < proteins.Count; i++)
            {
                var protein = proteins[i];
                if (protein.StartsWith("XXX") || protein.StartsWith("Contaminant"))
                {
                    continue;
                }

                var length = database.GetProteinLength(protein);
                //if (length <= 0)
                //{
                //    Console.WriteLine("Shit!");
                //    return;
                //}
                if (peptides != null)
                {
                    writer.Write(peptides[i] + "\t");
                }

                writer.Write(protein + "\t" + length);
                for (var j = 0; j < jobs.Length; j++)
                {
                    writer.Write("\t" + abundances[j][i]);
                }
                writer.WriteLine("\t" + (protein.StartsWith("STANDARD") ? 1 : 0));
            }
        }

        [Test]
        [Category("Local_Testing")]
        public void AddNaToTable()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string dir = @"H:\Research\IPRG2015";
            const string resultFilePath = dir + @"\AMT_Peptides_Missing.tsv";

            if (!File.Exists(resultFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, resultFilePath);
            }

            const string outputFilePath = dir + @"\AMT_Peptides_NA.tsv";

            using var writer = new StreamWriter(outputFilePath);

            foreach (var line in File.ReadLines(resultFilePath))
            {
                var token = line.Split('\t');
                writer.WriteLine(string.Join("\t", token.Select(t => t.Length == 0 ? "NA" : (Double.TryParse(t, out var result) ? (result * 1E6).ToString() : t))));
            }
        }

        [Test]
        [Category("Local_Testing")]
        public void AddProteinLengths()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string databaseFilePath = @"H:\Research\IPRG2015\database\yeast6proteaprotein.fasta";
            if (!File.Exists(databaseFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, databaseFilePath);
            }

            var database = new FastaDatabase(databaseFilePath);
            database.Read();

            const string resultPath = @"H:\Research\IPRG2015\AMT_Peptides_NA.tsv";
            if (!File.Exists(resultPath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, resultPath);
            }

            const string outputFilePath = @"H:\Research\IPRG2015\AMT_Peptides.tsv";

            using var writer = new StreamWriter(outputFilePath);

            foreach (var line in File.ReadLines(resultPath))
            {
                var data = line.Split(null);
                if (data.Length != 14)
                {
                    continue;
                }

                var peptide = data[0];
                if (peptide.Equals("Peptide"))
                {
                    writer.WriteLine("Peptide\tProtein\tLength\t{0}", string.Join("\t", data.Skip(2)));
                    continue;
                }
                var protein = data[1];
                var length = database.GetProteinLength(protein);
                writer.WriteLine("{0}\t{1}\t{2}\t{3}", peptide, protein, length, string.Join("\t", data.Skip(2)));
            }
        }

        [Test]
        [Category("Local_Testing")]
        public void GetProteinAccessions()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string uniprotAccession = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}";
            var uniProtPattern = new Regex(uniprotAccession);
            const string databaseFilePath = @"H:\Research\IPRG2015\Henry_results\iPRG2015.TargDecoy.fasta";
            if (!File.Exists(databaseFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, databaseFilePath);
            }

            var database = new FastaDatabase(databaseFilePath);
            database.Read();
            var nameToAccession = new Dictionary<string, string>();
            foreach (var proteinName in database.GetProteinNames())
            {
                var start = proteinName.IndexOf('|');
                var end = proteinName.LastIndexOf('|');
                //var accession = proteinName.Substring(start + 1, end - start - 1);
                var name = proteinName.Substring(end + 1);
                if (proteinName.StartsWith("DECOY"))
                {
                    name += "-DECOY";
                }
                //                Console.WriteLine(name + " -> " +accession);
                Assert.IsTrue(uniProtPattern.IsMatch(proteinName));
                nameToAccession.Add(name, proteinName);
                //                Console.WriteLine(name);
            }

            const string resultPath = @"H:\Research\IPRG2015\Henry_results\ProteinNames.txt";
            if (!File.Exists(resultPath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, resultPath);
            }

            foreach (var line in File.ReadLines(resultPath))
            {
                if (line.Length == 0)
                {
                    continue;
                }

                var name = line;
                //                if (name.Contains(";"))
                //                {
                //                }
                name = name.Split()[0];
                if (name.Contains('|'))
                {
                    name = name.Substring(name.LastIndexOf('|') + 1);
                }

                if (nameToAccession.TryGetValue(name, out var proteinName))
                {
                    Console.WriteLine(proteinName);
                }
                else
                {
                    Console.WriteLine(name);
                    Assert.IsTrue(false);
                }
            }
        }

        [Test]
        [Category("Local_Testing")]
        public void CreateTargetList()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string databaseFilePath = @"D:\Research\Data\IPRG2014\database\SpikedInPeptides.fasta";
            if (!File.Exists(databaseFilePath))
            {
                Assert.Ignore("Skipping test {0} since file not found: {1}", methodName, databaseFilePath);
            }

            var database = new FastaDatabase(databaseFilePath);
            database.Read();
            var indexedDatabase = new IndexedDatabase(database);
            var numTargets = 0;

            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);

            Console.WriteLine("Peptide\tFormula\tProtein");
            foreach (var annotationAndOffset in indexedDatabase.AnnotationsAndOffsets(6, 30, 1, 1, Enzyme.Trypsin))
            {
                var annotation = annotationAndOffset.Annotation;
                var peptide = annotation.Substring(2, annotation.Length - 4);
                var offset = annotationAndOffset.Offset;

                Console.WriteLine("{0}\t{1}\t{2}", peptide, (aaSet.GetComposition(peptide) + Composition.H2O).ToPlainString(), database.GetProteinName(offset));
                numTargets++;
            }
            Console.WriteLine("NumTargets: {0}", numTargets);
        }
    }
}
