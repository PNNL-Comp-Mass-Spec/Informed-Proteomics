using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.DevTests
{
    [TestFixture]
    public class TestMisc
    {
        [Test]
        public void RemovePepFdrFromFile()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string henryResultPath = @"H:\Research\IPRG2015\Henry_results\tsv";
            if (!Directory.Exists(henryResultPath))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, henryResultPath);
            }

            foreach (var resultFile in Directory.GetFiles(henryResultPath, "*_FDR.tsv"))
            {
                var fileName = Path.GetFileName(resultFile);
                if (fileName == null) continue;
                var sample = fileName.Substring(0, 2);
                Console.WriteLine("Processing {0}", sample);

                var outputFilePath = Path.Combine(henryResultPath,  "ID_" + sample + ".tsv");
                Console.WriteLine("Writing to {0}", outputFilePath);
                using (var writer = new StreamWriter(outputFilePath))
                {
                    foreach (var line in File.ReadLines(resultFile))
                    {
                        writer.WriteLine(line.Substring(0, line.LastIndexOf('\t')));
                    }
                }
            }
        }

        [Test]
        public void CreatePeptideAbundanceTableWithSkyline()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            // Reading Henry's results
            var pepKeySet = new HashSet<string>();
            var resultDic = new Dictionary<string, Tuple<double, double>>();
            const string henryResultPath = @"H:\Research\IPRG2015\Henry_results\tsv";
            if (!Directory.Exists(henryResultPath))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, henryResultPath);
            }

            var aaSet = new AminoAcidSet();
            foreach (var resultFile in Directory.GetFiles(henryResultPath, "*.tsv"))
            {
                var fileName = Path.GetFileName(resultFile);
                if (fileName == null) continue;
                var sample = fileName.Substring(0, 2);
                Console.WriteLine("Processing {0}", sample);
                var tsvReader = new TsvFileParser(resultFile);
                var peptides = tsvReader.GetData("Peptide").ToArray();
                var charge = tsvReader.GetData("Charge").Select(c => Convert.ToInt32(c)).ToArray();
                var prob = tsvReader.GetData("Prob").Select(Convert.ToDouble).ToArray();
                var qValue = tsvReader.GetData("QValue").Select(Convert.ToDouble).ToArray();
                for (var i = 0; i < tsvReader.NumData; i++)
                {
                    var peptide = peptides[i];
                    var nominalMass = GetNominalMass(aaSet, peptide);
                    var key = sample + ":" + GetPeptide(peptides[i]) + ":" + nominalMass + ":" + charge[i];
                    var pepKey = GetPeptide(peptides[i]) + ":" + nominalMass;
                    pepKeySet.Add(pepKey);
                    Tuple<double, double> existingScores;
                    if (resultDic.TryGetValue(key, out existingScores))
                    {
                        if (prob[i] > existingScores.Item1)
                        {
                            resultDic[key] = new Tuple<double, double>(prob[i], qValue[i]);
                        }
                    }
                    else
                    {
                        resultDic.Add(key, new Tuple<double, double>(prob[i], qValue[i]));
                    }
                }
            }

            const string skylineFilePath = @"H:\Research\IPRG2015\MySkyline\TransitionResults.csv";
            if (!File.Exists(skylineFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, skylineFilePath);
            }

            var skylineTable = new TsvFileParser(skylineFilePath, ',');

            const string outputFilePath = @"H:\Research\IPRG2015\MySkyline\SkylineTransitionResultsWithScores3.tsv";
            using (var writer = new StreamWriter(outputFilePath))
            {
                var peptides = skylineTable.GetData("Peptide Sequence").ToArray();
                var samples = skylineTable.GetData("Replicate Name").Select(s => "" + s[0] + s[2]).ToArray();
                var charges = skylineTable.GetData("Precursor Charge").Select(c => Convert.ToInt32(c)).ToArray();
                var precursorMzs = skylineTable.GetData("Precursor Mz").Select(Convert.ToDouble).ToArray();

                writer.WriteLine("{0}\tProbability\tQValue", string.Join("\t", skylineTable.GetHeaders().Take(skylineTable.GetHeaders().Count-2)));
                for (var i = 0; i < skylineTable.NumData; i++)
                {
                    var precursorMz = precursorMzs[i];
                    var charge = charges[i];
                    var nominalMass = (int)Math.Round(((precursorMz - Constants.Proton)*charge - Composition.H2O.Mass)*
                                      Constants.RescalingConstant);
                    var pepKey = peptides[i] + ":" + nominalMass;
                    if (!pepKeySet.Contains(pepKey))
                    {
                        //Console.WriteLine("Removing {0}", pepKey);
                        continue;
                    }
                    var key = samples[i] + ":" + peptides[i] + ":" + nominalMass + ":" + charge;
                    double? prob = null, qValue = null;
                    Tuple<double, double> scores;
                    if (resultDic.TryGetValue(key, out scores))
                    {
                        prob = scores.Item1;
                        qValue = scores.Item2;
                    }
                    var skylineData = skylineTable.GetRows()[i].Split(',');
                    for (var j = 0; j < skylineData.Length - 2; j++)
                    {
                        if(j != 2) writer.Write(skylineData[j]+"\t");
                        else writer.Write("" + skylineData[j][0] + skylineData[j][2]+"\t");
                    }
                    writer.WriteLine("{0}\t{1}",
                        prob != null ? prob.ToString() : "NA",
                        qValue != null ? qValue.ToString() : "NA");
                }
            }
            Console.WriteLine("Done");
        }

        private int GetNominalMass(AminoAcidSet aaSet, string peptide)
        {
            var nominalMass = 0;
            StringBuilder buf = null;
            var curNominalMass = 0;
            foreach (var c in peptide)
            {
                if (char.IsUpper(c))
                {
                    curNominalMass = aaSet.GetAminoAcid(c).GetNominalMass();
                    nominalMass += curNominalMass;
                }
                else if (c == '[') buf = new StringBuilder();
                else if (char.IsNumber(c) && buf != null) buf.Append(c);
                else if (c == ']' && buf != null) nominalMass += Convert.ToInt32(buf.ToString()) - curNominalMass;
            }
            return nominalMass;
        }

        private string GetPeptide(string henryPeptide)
        {
            var buf = new StringBuilder();
            foreach(var c in henryPeptide) if (char.IsUpper(c)) buf.Append(c);
            return buf.ToString();
        }

        [Test]
        public void TestPathUtils()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            const string rawFilePath = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\Synocho_D1_1.raw";
            if (!File.Exists(rawFilePath))
            {
                Assert.Ignore(@"Skipping test {0} since file not found: {1}", methodName, rawFilePath);
            }

            Console.WriteLine(MassSpecDataReaderFactory.RemoveExtension(rawFilePath) + "_Target.tsv");
            Console.WriteLine(Path.GetDirectoryName(rawFilePath));
            Console.WriteLine(Path.Combine(Path.GetDirectoryName(rawFilePath), Path.GetFileNameWithoutExtension(rawFilePath) + "_IcTarget.tsv"));

            var outputDir = @"C:\cygwin\home\kims336\Data\TopDownJia\raw\L1_1_Mode2\Synocho_L1_1_IcTarget.tsv";

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
