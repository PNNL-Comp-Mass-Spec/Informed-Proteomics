using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Database;
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
        public void AddNaToTable()
        {
            const string dir = @"D:\Research\Data\IPRG2014";
            const string resultFilePath = dir + @"\AMTAllPeptidesMissingValues.tsv";
            foreach (var line in File.ReadLines(resultFilePath))
            {
                var token = line.Split('\t');
                Console.WriteLine(string.Join("\t", token.Select(t => t.Length == 0 ? "NA" : t)));
            }
        }

        [Test]
        public void AddIPrgProteinLengths()
        {
            const string dir = @"D:\Research\Data\IPRG2014";
            const string databaseFilePath = dir + @"\database\E_coli_K12_uniprot_reviewed_2013-01-31.fasta";
            var database = new FastaDatabase(databaseFilePath);
            database.Read();

            const string proteinListPath = dir + @"\ProteinListIntensity.txt";
            foreach (var proteinId in File.ReadLines(proteinListPath))
            {
                var length = database.GetProteinLength(proteinId);
                Console.WriteLine("{0}\t{1}", proteinId, length);
            }
        }

        [Test]
        public void GenerateAbrfSpecCountAllProteins()
        {
            const string dir = @"D:\Research\Data\IPRG2014";
            const double qValueThreshold = 0.01;
            //var names = new[] { "ENO1_YEAST", "ADH1_YEAST", "CYC_BOVIN", "ALBU_BOVIN" };
            //var accessions = new[] { "P00924", "P00330", "P62894", "P02769" };

            const string resultDir = dir + @"\10ppm_TI0_NTT1";
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
                        var proteinName = protein.Substring(0,
                                                            protein.LastIndexOf("(pre=", System.StringComparison.Ordinal));
                        int[] countArr;
                        if (!specCount.TryGetValue(proteinName, out countArr)) specCount[proteinName] = new int[msgfResultFiles.Length];
                        specCount[proteinName][i]++;
                    }
                }
            }

            // Writing
            const string databaseFilePath = dir + @"\database\E_coli_K12_uniprot_reviewed_2013-01-31.fasta";
            var database = new FastaDatabase(databaseFilePath);
            database.Read();

            var spikeInAccessions = new[] { "P00924", "P00330", "P62894", "P02769"};

            const string outputFilePath = dir + @"\SpecCountAllProteins.tsv";
            using (var writer = new StreamWriter(outputFilePath))
            {
                var fileIds = msgfResultFiles.Select(f => f.Substring(f.LastIndexOf('_') + 1, 3));
                writer.WriteLine("Protein\tLength\t" + string.Join("\t", fileIds) + "\tSpikeIn");
                foreach (var entry in specCount)
                {
                    var proteinId = entry.Key;
                    var length = database.GetProteinLength(proteinId);
                    Assert.True(length > 0);
                    var counts = entry.Value;
                    Assert.True(counts.Length == msgfResultFiles.Length);
                    var spikeIn = 0;
                    if (spikeInAccessions.Any(spikeInAccession => proteinId.StartsWith("sp|" + spikeInAccession)))
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
    }
}
