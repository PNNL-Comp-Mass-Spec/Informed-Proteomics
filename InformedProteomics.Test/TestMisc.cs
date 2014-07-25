using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.Utils;
using InformedProteomics.DIA.Search;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestMisc
    {
        [Test]
        public void TestDivisionByInfinity()
        {
            Console.WriteLine(10f/float.NegativeInfinity);
        }

        [Test]
        public void ProcessIprg2015PreStudy()
        {
            const string dir = @"H:\Research\IPRG2015";

            const string databaseFilePath = dir + @"\yeast6mix.fasta";
            var database = new FastaDatabase(databaseFilePath);
            database.Read();

            const string jobFilePath = dir + @"\Jobs.tsv";
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
                    if (headers[j].Contains("" + jobs[i]))
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

            using (var writer = new StreamWriter(outputFilePath))
            {
                var peptides = parser.GetData("Peptide");   // Peptides
                var proteins = parser.GetData("Reference");     // Proteins
                var abundances = new string[jobs.Length][];
                for (var i = 0; i < jobs.Length; i++)
                {
                    abundances[i] = parser.GetData(headers[jobColNum[i]]).ToArray();
                }

                if(peptides != null) writer.Write("Peptide\t");
                writer.Write("Protein\tLength");
                for (var i = 0; i < jobs.Length; i++)
                {
                    writer.Write("\t"+experiments[i]);
                }
                writer.WriteLine("\tSpikeIn");
                for (var i = 0; i < proteins.Count; i++)
                {
                    var protein = proteins[i];
                    if (protein.StartsWith("XXX") || protein.StartsWith("Contaminant")) continue;
                    var length = database.GetProteinLength(protein);
                    //if (length <= 0)
                    //{
                    //    Console.WriteLine("Shit!");
                    //    return;
                    //}
                    if(peptides != null) writer.Write(peptides[i]+"\t");
                    writer.Write(protein+"\t"+length);
                    for (var j = 0; j < jobs.Length; j++)
                    {
                        writer.Write("\t"+abundances[j][i]);
                    }
                    writer.WriteLine("\t"+ (protein.StartsWith("STANDARD") ? 1 : 0));
                }
            }
        }

        [Test]
        public void AddNaToTable()
        {
            const string dir = @"H:\Research\IPRG2015";
            const string resultFilePath = dir + @"\AMT_Peptides.tsv";

            const string outputFilePath = dir + @"\AMT_Peptides_NA.tsv";
            using (var writer = new StreamWriter(outputFilePath))
            {
                foreach (var line in File.ReadLines(resultFilePath))
                {
                    var token = line.Split('\t');
                    writer.WriteLine(string.Join("\t", token.Select(t => t.Length == 0 ? "NA" : t)));
                }
            }
        }

        [Test]
        public void AddProteinLengths()
        {
            const string databaseFilePath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\database\ID_002216_235ACCEA.fasta";
            var database = new FastaDatabase(databaseFilePath);
            database.Read();

            const string proteinListPath = @"C:\cygwin\home\kims336\Data\TopDownQCShew\MSAlign\Proteins.txt";
            foreach (var protein in File.ReadLines(proteinListPath))
            {
                var proteinId = protein.Split(null)[0];
                var length = database.GetProteinLength(proteinId);
                Console.WriteLine("{0}\t{1}", proteinId, length);
            }
        }

        [Test]
        public void GenerateAbrfSpecCountAllProteins()
        {
            const string dir = @"H:\Research\IPRG2015";
            const double qValueThreshold = 0.01;
            //var names = new[] { "ENO1_YEAST", "ADH1_YEAST", "CYC_BOVIN", "ALBU_BOVIN" };
            //var accessions = new[] { "P00924", "P00330", "P62894", "P02769" };

            const string resultDir = dir + @"\NTT1";
            var msgfResultFiles = Directory.GetFiles(resultDir, "*.tsv").ToArray();

            var specCount = new Dictionary<string, int[]>();  // protein name => array of counts

            for (var i = 0; i < msgfResultFiles.Length; i++ )
            {
                var msgfResultFile = msgfResultFiles[i];

                MsGfPlusHeaderInformation headerInfo = null;

                var prevScanNum = -1;
                foreach (var line in File.ReadLines(msgfResultFile))
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
                    if (match.QValue > qValueThreshold) continue;

                    var proteins = match.Protein.Split(';');
                    foreach (var protein in proteins)
                    {
                        var proteinName = protein.Substring(0, protein.LastIndexOf("(pre=", StringComparison.Ordinal));
                        int[] countArr;
                        if (!specCount.TryGetValue(proteinName, out countArr)) specCount[proteinName] = new int[msgfResultFiles.Length];
                        specCount[proteinName][i]++;
                    }
                }
            }

            // Writing
            const string databaseFilePath = dir + @"\database\iPRG2015.fasta";
            var database = new FastaDatabase(databaseFilePath);
            database.Read();

//            var spikeInAccessions = new[] { "STANDARD_Alpha-Casein", "STANDARD_Beta-Lactoglobulin", "STANDARD_Carbonic-Anhydrase", "P02769"};

            const string outputFilePath = dir + @"\SpecCountAllProteins.tsv";
            using (var writer = new StreamWriter(outputFilePath))
            {
                var fileIds = msgfResultFiles.Select(f => f.Substring(f.IndexOf("_sample", StringComparison.Ordinal) + 1, 
                    f.LastIndexOf('.') - f.IndexOf("_sample", StringComparison.Ordinal) - 1));
                writer.WriteLine("Protein\tLength\t" + string.Join("\t", fileIds) + "\tSpikeIn");
                foreach (var entry in specCount)
                {
                    var proteinId = entry.Key;
                    var length = database.GetProteinLength(proteinId);
                    Assert.True(length > 0);
                    var counts = entry.Value;
                    Assert.True(counts.Length == msgfResultFiles.Length);
                    var spikeIn = 0;
                    //if (spikeInAccessions.Any(spikeInAccession => proteinId.StartsWith("sp|" + spikeInAccession)))
                    if(proteinId.StartsWith("sp|"))
                    {
                        spikeIn = 1;
                    }
                    writer.WriteLine("{0}\t{1}\t{2}\t{3}", proteinId, length, string.Join("\t", counts), spikeIn);
                }
            }
        }

        [Test]
        public void TestAbrfSpecCount()
        {
            const string dir = @"D:\Research\Data\IPRG2014";
            const double qValueThreshold = 0.01;
            var names = new[] {"ENO1_YEAST", "ADH1_YEAST", "CYC_BOVIN", "ALBU_BOVIN"};
            var accessions = new[] { "P00924", "P00330", "P62894", "P02769"};

            const string databaseFilePath = dir + @"\database\E_coli_K12_uniprot_reviewed_2013-01-31.revCat.fasta";
            const string resultDir = dir + @"\10ppm_TI0_NTT1";
            
            Console.WriteLine("Run\tTotal PSM\t"+string.Join("\t", names));

            foreach (var msgfResultFile in Directory.GetFiles(resultDir, "*.tsv"))
            {
                var fileId = msgfResultFile.Substring(msgfResultFile.LastIndexOf('_') + 1, 3);
                Console.Write(fileId);

                var totalPsm = 0;
                var specCount = new int[accessions.Length];

                MsGfPlusHeaderInformation headerInfo = null;

                var prevScanNum = -1;
                foreach (var line in File.ReadLines(msgfResultFile))
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
                    if (match.QValue > qValueThreshold) continue;

                    totalPsm++;

                    for (var i = 0; i < accessions.Length; i++)
                    {
                        if (match.Protein.StartsWith("sp|" + accessions[i]))
                        {
                            specCount[i]++;
                        }
                    }
                }

                Console.Write("\t"+totalPsm);
                for (var i = 0; i < accessions.Length; i++)
                {
                    Console.Write("\t{0}", specCount[i]);
                }
                Console.WriteLine();
            }
        }

        [Test]
        public void CreateTargetList()
        {
            const string databaseFilePath = @"D:\Research\Data\IPRG2014\database\SpikedInPeptides.fasta";
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

                Console.WriteLine("{0}\t{1}\t{2}", peptide, (aaSet.GetComposition(peptide)+Composition.H2O).ToPlainString(), database.GetProteinName(offset));
                numTargets++;
            }
            Console.WriteLine("NumTargets: {0}", numTargets);
        }

        [Test]
        public void TestPathUtils()
        {
            const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1.raw";
            Console.WriteLine(Path.ChangeExtension(rawFilePath, "_Target.tsv"));
            Console.WriteLine(Path.GetDirectoryName(rawFilePath));
            Console.WriteLine(Path.GetDirectoryName(rawFilePath) + Path.DirectorySeparatorChar + Path.GetFileNameWithoutExtension(rawFilePath)+"_IcTarget.tsv");

            var outputDir = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\L1_1_Mode2\Synocho_L1_1_IcTarget.tsv";
            if (outputDir[outputDir.Length - 1] == Path.DirectorySeparatorChar) outputDir = outputDir.Remove(outputDir.Length - 1);
            if (!Directory.Exists(outputDir))
            {
                if (!File.GetAttributes(outputDir).HasFlag(FileAttributes.Directory))
                {
                    throw new Exception(outputDir + " is not a directory!");
                } 
                Directory.CreateDirectory(outputDir);
            }
            Console.WriteLine(outputDir);

        }
    }
}
