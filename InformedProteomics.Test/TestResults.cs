using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestResults
    {
        [Test]
        public void SummarizeAnilResults()
        {
            const string resultFolder = @"H:\Research\Anil\Oct28";
            var actMethods = new[] { ActivationMethod.CID, ActivationMethod.ETD, ActivationMethod.HCD };

            Console.WriteLine("Data\tCID\t\tETD\t\tHCD\t");
            Console.WriteLine("\tNumId\tMaxMass\tNumId\tMaxMass\tNumId\tMaxMass");
            foreach (var rawFile in Directory.GetFiles(resultFolder, "*.raw"))
            {
                var datasetName = Path.GetFileNameWithoutExtension(rawFile);

                var resultFile = Path.GetDirectoryName(rawFile) + Path.DirectorySeparatorChar + datasetName + "_IcTda.tsv";
                var numId = new Dictionary<ActivationMethod, int>();
                var maxMass = new Dictionary<ActivationMethod, double>();

                foreach (var actMethod in actMethods)
                {
                    numId[actMethod] = 0;
                    maxMass[actMethod] = 0.0;
                }

                var run = PbfLcMsRun.GetLcMsRun(rawFile);
                var tsvParser = new TsvFileParser(resultFile);
                var qValues = tsvParser.GetData("QValue").Select(Convert.ToDouble).ToArray();
                var scanNums = tsvParser.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
                var masses = tsvParser.GetData("Mass").Select(Convert.ToDouble).ToArray();

                for (var i = 0; i < qValues.Length; i++)
                {
                    if (qValues[i] > 0.01) break;
                    var scanNum = scanNums[i];
                    var spec = run.GetSpectrum(scanNum) as ProductSpectrum;
                    Assert.True(spec != null);
                    ++numId[spec.ActivationMethod];

                    var mass = masses[i];
                    if (mass > maxMass[spec.ActivationMethod]) maxMass[spec.ActivationMethod] = mass;
                }
                Console.Write(datasetName);
                foreach (var actMethod in actMethods)
                {
                    Console.Write("\t" + numId[actMethod]);
                    Console.Write("\t" + maxMass[actMethod]);
                }
                Console.WriteLine();
            }
        }

        [Test]
        public void CountIdentifiedPeptides()
        {
            //const string targetResultPath = @"H:\Research\Jarret\10fmol_10mz\NoMod_NTT2\Q_2014_0523_50_10_fmol_uL_10mz_IcTarget.tsv";
            //const string decoyResultPath = @"H:\Research\Jarret\10fmol_10mz\NoMod_NTT2\Q_2014_0523_50_10_fmol_uL_10mz_IcDecoy.tsv";
            //const string path = @"H:\Research\Jarret\10mz\0amol\Q_2014_0523_10_0_amol_uL_10mz";
            const string path = @"H:\Research\Jarret\20mz\NTT2\Q_2014_0523_12_0_amol_uL_20mz";
            //const string path = @"H:\Research\Jarret\DDA\NoMod_NTT2\Q_2014_0523_1_0_amol_uL_DDA";
            const string targetResultPath = path + "_IcTarget.tsv";
            const string decoyResultPath = path + "_IcDecoy.tsv";

            var fdrCalculator = new FdrCalculator(targetResultPath, decoyResultPath, false);         
            if (fdrCalculator.HasError())
            {
                throw new Exception(@"Error computing FDR: " + fdrCalculator.ErrorMessage);
            }

            Console.WriteLine("NumPSMs: {0}", fdrCalculator.NumPsms);
            Console.WriteLine("NumPeptides: {0}", fdrCalculator.NumPeptides);
            Console.WriteLine("Done");
        }

        [Test]
        public void GenerateVennDiagrams()
        {
            // DIA
            const string dir = @"H:\Research\DDAPlus\NTT2";

            const string dda1 = dir + @"\20140701_yeast_DDA_01_IcTda.tsv";
            const string dda2 = dir + @"\20140701_yeast_DDA_02_2_IcTda.tsv";
            const string ddaPlus1 = dir + @"\20140701_yeast_DDAp_binCharge_01_IcTda.tsv";
            const string ddaPlus2 = dir + @"\20140701_yeast_DDAp_binCharge_02_IcTda.tsv";

            const string resultPath1 = ddaPlus1;
            const string resultPath2 = ddaPlus2;

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

        [Test]
        public void TestIcrTools()
        {
            const string icrToolsPath = @"H:\Research\Yufeng\TopDownYufeng\ICRTools\yufeng_column_test2_Isos.csv";
            const double fitScoreThreshold = 0.15;
            var icrToolsparser = new TsvFileParser(icrToolsPath, ',');
            var monoMassArr = icrToolsparser.GetData("monoisotopic_mw").Select(Convert.ToDouble).ToArray();
            var scanArray = icrToolsparser.GetData("scan_num").Select(s => Convert.ToInt32(s)).ToArray();
            var fitArray = icrToolsparser.GetData("fit").Select(Convert.ToDouble).ToArray();
            var peakList = fitArray.Where(fit => fit < fitScoreThreshold).Select((fit, i) => new Peak(monoMassArr[i], scanArray[i])).OrderBy(p => p.Mz).ToList();
            Console.WriteLine("FitScoreThreshold: {0}", fitScoreThreshold); 
            Console.WriteLine("NumFeatures: {0}", peakList.Count);
            //var peakList = monoMassArr.Select((mass, i) => new Peak(mass, scanArray[i])).OrderBy(p => p.Mz).ToList();
            var massTolerance = new Tolerance(10);
            const double scanTolerance = 10.0;

            const string resultFilePath = @"H:\Research\Yufeng\TopDownYufeng\M1_V31\yufeng_column_test2_IcTda.tsv";
            var parser = new TsvFileParser(resultFilePath);
            var qValues = parser.GetData("QValue").Select(Convert.ToDouble).ToArray();
            var masses = parser.GetData("Mass").Select(Convert.ToDouble).ToArray();
            var scans = parser.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
            var compositions = parser.GetData("Composition").ToArray();
            var numTotal = 0;
            var numMatch = 0;
            var compSet = new HashSet<string>();
            var matchedCompSet = new HashSet<string>();
            for (var i = 0; i < qValues.Length; i++)
            {
                if (qValues[i] > 0.01) break;
                var mass = masses[i];
                var scan = scans[i];
                var matches = PeakListUtils.FindAllPeaks(peakList, mass, massTolerance)
                    .Where(p => p.Intensity >= scan - scanTolerance && p.Intensity <= scan + scanTolerance);
                if (matches.Any())
                {
                    ++numMatch;
                    matchedCompSet.Add(compositions[i]);
                }
                ++numTotal;
                compSet.Add(compositions[i]);
            }
            Console.WriteLine("PrSMs: {0} / {1} = {2}", numMatch, numTotal, numMatch/(float)numTotal);
            Console.WriteLine("Compositions: {0} / {1} = {2}", matchedCompSet.Count, compSet.Count, matchedCompSet.Count / (float)compSet.Count);
        }
    }
}
