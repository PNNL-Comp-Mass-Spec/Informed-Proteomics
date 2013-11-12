using System;
using System.IO;
using InformedProteomics.Backend.Database;
using InformedProteomics.DIA.Search;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestMisc
    {
        private const double QValueThreshold = 0.01;
        [Test]
        public void TestAbrfSpecCount()
        {
            var names = new[] {"ENO1_YEAST", "ADH1_YEAST", "CYC_BOVIN", "ALBU_BOVIN"};
            var accessions = new[] { "P00924", "P00330", "P62894", "P02769"};

            const string dir = @"D:\Research\Data\IPRG2014\10ppm_TI0_NTT1";
            
            Console.WriteLine("Run\tTotal PSM\t"+string.Join("\t", names));
            
            foreach (var msgfResultFile in Directory.GetFiles(dir, "*.tsv"))
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
                    if (match.QValue > QValueThreshold) continue;

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
    }
}
