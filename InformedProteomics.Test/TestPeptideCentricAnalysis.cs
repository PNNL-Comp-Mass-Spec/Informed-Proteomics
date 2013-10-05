using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.DIA.Search;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestPeptideCentricAnalysis
    {
        [Test]
        public void TestQExactiveDdaDataPostProcessing()
        {
            //const string range = "650to775";
            //const string resultPath = @"D:\Research\Data\UW\QExactive\DIA_Results\DIA_" + range + ".tsv";
            //const string specFilePath = @"D:\Research\Data\UW\QExactive\82593_lv_mcx_DIA_5mz_" + range + ".raw";
            //const string outputFilePath = @"D:\Research\Data\UW\QExactive\DIA_" + range + "_Summary.tsv";

            const string resultPath = @"D:\Research\Data\UW\QExactive\DDA_Results\DDA_All.tsv";
            const string specFile = @"D:\Research\Data\UW\QExactive\82593_lv_mcx_DDA_NoCharge.raw";
            const string outputFilePath = @"D:\Research\Data\UW\QExactive\DDA_All_Summary.tsv";

            var postProcessor = new MsGfPostProcessor(specFile, resultPath, new Tolerance(20), new Tolerance(10));
            var numId = postProcessor.PostProcessing(outputFilePath);

            Console.WriteLine("NumId: {0}", numId);
        }

        [Test]
        public void TestQExactiveDiaDataPostProcessing()
        {
            //const string range = "650to775";
            //const string resultPath = @"D:\Research\Data\UW\QExactive\DIA_Results\DIA_" + range + ".tsv";
            //const string specFilePath = @"D:\Research\Data\UW\QExactive\82593_lv_mcx_DIA_5mz_" + range + ".raw";
            //const string outputFilePath = @"D:\Research\Data\UW\QExactive\DIA_" + range + "_Summary.tsv";

            const string resultPath = @"D:\Research\Data\UW\QExactive\DIA_Results\DIA_All.tsv";
            var specFiles = Directory.GetFiles(@"D:\Research\Data\UW\QExactive\", "*_DIA_*.raw");
            const string outputFilePath = @"D:\Research\Data\UW\QExactive\DIA_All_Summary.tsv";

            var postProcessor = new MsGfPostProcessor(specFiles, resultPath, new Tolerance(20), new Tolerance(10));
            var numId = postProcessor.PostProcessing(outputFilePath);

            Console.WriteLine("NumId: {0}", numId);
        }

        [Test]
        public void TestFusionDiaDataPostProcessing()
        {
            const string resultPath = @"D:\Research\Data\UW\Fusion\DIA_Results\DIA_C_1_To_4.tsv";
            const string specFilePath = @"D:\Research\Data\UW\Fusion\WT_D_DIA_130412091220.raw";
            const string outputFilePath = @"D:\Research\Data\UW\Fusion\DIA_Summary.tsv";

            var postProcessor = new MsGfPostProcessor(specFilePath, resultPath, new Tolerance(5), new Tolerance(3));
            var numId = postProcessor.PostProcessing(outputFilePath);
            
            Console.WriteLine("NumId: {0}", numId);
        }

        [Test]
        public void TestFusionDdaDataPostProcessing()
        {
            const string resultPath = @"D:\Research\Data\UW\Fusion\DDA_Results\DDA_C_1_To_4.tsv";
            const string specFilePath = @"D:\Research\Data\UW\Fusion\WT_D_DDA_130412065618.raw";
            const string outputFilePath = @"D:\Research\Data\UW\Fusion\DDA_Summary.tsv";

            var postProcessor = new MsGfPostProcessor(specFilePath, resultPath, new Tolerance(5), new Tolerance(3));
            var numId = postProcessor.PostProcessing(outputFilePath);

            Console.WriteLine("NumId: {0}", numId);
        }

        [Test]
        public void GenerateVennDiagrams()
        {
            // Fusion
            const string fusionMsgfResult = @"D:\Research\Data\UW\Fusion\MSGFPlusResults\TI2.tsv";
            const string fusionDdaResult = @"D:\Research\Data\UW\Fusion\DDA_Summary.tsv";
            const string fusionDiaResult = @"D:\Research\Data\UW\Fusion\DIA_Summary.tsv";

            // Q-Exactive
            const string qeDdaResult = @"D:\Research\Data\UW\QExactive\DDA_Results\82593_lv_mcx_DDA.tsv";
            const string qeDiaResult = @"D:\Research\Data\UW\QExactive\DIA_All_Summary.tsv";

            const string resultPath1 = fusionDiaResult;
            const string resultPath2 = fusionMsgfResult;

            var result1 = new TsvFileParser(resultPath1);
            var result2 = new TsvFileParser(resultPath2);

            const double pepQValueThreshold = 0.01;
            var vennDiagram = new VennDiagram<string>(result1.GetPeptides(pepQValueThreshold),
                                                      result2.GetPeptides(pepQValueThreshold));
            Console.WriteLine("{0}\t{1}\t{2}", vennDiagram.Set1Only.Count, vennDiagram.Intersection.Count, vennDiagram.Set2Only.Count);
        }

        [Test]
        public void TestMacCossFusionData()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            const string specFilePath = @"D:\Research\Data\UW\Fusion\WT_D_DDA_130412065618.raw";

            const string dbFilePath = @"D:\Research\Data\UW\Database\Yeast_SGD_withContam.revCat.fasta";
            var targetDecoyDb = new FastaDatabase(dbFilePath);
            var indexedDb = new IndexedDatabase(targetDecoyDb);
            var sequences = indexedDb.SequencesAsStrings(minLength: 6, maxLength: 40, numTolerableTermini: 1,
                                                         numMissedCleavages: 2, enzyme: Enzyme.Trypsin);

            var aminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);

            var run = new LcMsRun(new XCaliburReader(specFilePath));

            var analyzer = new PeptideCentricAnalysis(run, sequences, aminoAcidSet);


            const string outputFilePath = @"D:\Research\Data\UW\Fusion\WT_D_DIA_130412091220.pca.txt";
            analyzer.IntermediateSearch(outputFilePath);

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);
        }

        [Test]
        public void TestXicGen()
        {
            const string specFilePath = @"D:\Research\Data\UW\Fusion\WT_D_DDA_130412065618.raw";
            var run = new LcMsRun(new XCaliburReader(specFilePath));
            // Test
            var tolerance = new Tolerance(30);

            const string peptide = "AIANGQVDGFPTQEECR";
            const int targetScanNum = 37633;
            const int charge = 2;

            //const string peptide = "IVDTNGAGDAFAGGFMAGLTK";
            //const int targetScanNum = 67513;
            //const int charge = 3;

            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var precursorIon = new Ion(aaSet.GetComposition(peptide) + Composition.H2O, charge);

            Console.WriteLine("Theoretical isotopomer profile:");
            foreach(var p in precursorIon.GetIsotopes(0.1)) Console.WriteLine("{0}\t{1}", precursorIon.GetIsotopeMz(p.Item1), p.Item2);

            var xicArr = new Dictionary<int, Xic>();
            var basePeakIndex = precursorIon.Composition.GetMostAbundantIsotopeZeroBasedIndex();
            for (var i = -1; i < 3; i++)
            {
                xicArr[i] = run.GetExtractedIonChromatogram(precursorIon.GetIsotopeMz(i), tolerance, targetScanNum);
            }

            for (var i = -1; i < 3; i++)
            {
                Console.WriteLine("\nIndex: {0}", i);
                Console.WriteLine("m/z: {0}", precursorIon.GetIsotopeMz(i));
                Console.WriteLine("#XicPeaks: {0}", xicArr[i].Count);
                Console.WriteLine("Intensity: {0}", xicArr[i].GetSumIntensities()/xicArr[basePeakIndex].GetSumIntensities());
                Console.WriteLine("Correlation: {0}", xicArr[i].GetCorrelation(xicArr[basePeakIndex]));
            }
        }

        [Test]
        public void TestFusionDdaData()
        {
            // Parameters
            const double relativeIntensityThreshold = 0.7;
            const double precursorTolerancePpm = 20;
            const double isotopeRatioTolerance = 2;
            const double correlationThreshold = 0.3;
            const double fdrThreshold = 0.01;

            const string specFilePath = @"D:\Research\Data\UW\Fusion\WT_D_DDA_130412065618.raw";
            var run = new LcMsRun(new XCaliburReader(specFilePath));

            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            var tolerance = new Tolerance(precursorTolerancePpm);
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);

            const string resultFilePath = @"D:\Research\Data\UW\Fusion\oldResult\WT_D_DDA_130412065618_10ppm_TI2_SGD_Decoy.tsv";
            var numTargets = 0;
            var numValidTargets = 0;
            var numDecoys = 0;
            var numValidDecoys = 0;

            foreach(var line in File.ReadLines(resultFilePath))
            {
                if (line.StartsWith("#")) continue;
                var token = line.Split('\t');
                if (token.Length != 16) continue;

                var qValue = Convert.ToDouble(token[14]);
                if (qValue > fdrThreshold) continue;

                var peptide = token[8].Replace("C+57.021", "C");
                var scanNum = Convert.ToInt32(token[2]);
                var charge = Convert.ToInt32(token[7]);
                var protein = token[9];
                var isDecoy = protein.StartsWith("XXX_");
                if (isDecoy) numDecoys++;
                else numTargets++;

                var precursorIon = new Ion(aaSet.GetComposition(peptide) + Composition.H2O, charge);
                var basePeakIndex = precursorIon.Composition.GetMostAbundantIsotopeZeroBasedIndex();
                var baseXic = run.GetExtractedIonChromatogram(precursorIon.GetBaseIsotopeMz(), tolerance, scanNum);
                var baseIntensity = baseXic.GetSumIntensities();

                var isValid = true;
                foreach (var isotope in precursorIon.GetIsotopes(relativeIntensityThreshold))
                {
                    if (isotope.Item1 == basePeakIndex) continue;
                    var isotopeMz = precursorIon.GetIsotopeMz(isotope.Item1);
                    var xic = run.GetExtractedIonChromatogram(isotopeMz, tolerance, scanNum);

                    if (xic.Count == 0)
                    {
                        isValid = false;
                        break;
                    }

                    //if (xic.Count > 0)
                    //{
                    //    var isotopeRatio = xic.GetSumIntensities() / baseIntensity / isotope.Item2;
                    //    var correlation = xic.GetCorrelation(baseXic);
                        
                    //    if (isotopeRatio > 0.8 && isotopeRatio < 1.2
                    //        && correlation > 0.8)
                    //    {
                    //        isValid = true;
                    //    }
                    //}

                    // Check if isotope ratio is within tolerance
                    //if (isotopeRatio > isotopeRatioTolerance || isotopeRatio < 1 / isotopeRatioTolerance)
                    //{
                    //    isValid = false;
                    //    //Console.WriteLine("Off ratio\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", isDecoy, peptide, scanNum, charge, precursorIon.GetMz(), isotopeMz, isotopeRatio);
                    //    break;
                    //}

                    // Check if correlation is high
                    //if (correlation < correlationThreshold)
                    //{
                    //    isValid = false;
                    //    //Console.WriteLine("Low correlation\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", isDecoy, peptide, scanNum, charge, precursorIon.GetMz(), isotopeMz, correlation);
                    //    break;
                    //}
                }

                if (isValid && !isDecoy) numValidTargets++;
                else if (isValid) numValidDecoys++;

                //Console.WriteLine("{0}\t{1}\t{2}", peptide, scanNum, charge);
            }
            Console.WriteLine("#Targets: {0}", numTargets);
            Console.WriteLine("#ValidTargets: {0}\t{1}", numValidTargets, numValidTargets/(double)numTargets);
            Console.WriteLine("#Decoys: {0}", numDecoys);
            Console.WriteLine("#ValidDecoys: {0}\t{1}", numValidDecoys, numValidDecoys / (double)numDecoys);

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"TimeForPrecursorValidation {0:f4} sec", sec);
        }

        [Test]
        public void AnalyizeFusionDdaData()
        {
            // Parameters
            const double relativeIntensityThreshold = 0.7;
            const double precursorTolerancePpm = 20;

            const string specFilePath = @"D:\Research\Data\UW\Fusion\WT_D_DDA_130412065618.raw";
            var run = new LcMsRun(new XCaliburReader(specFilePath));
            const double fdrThreshold = 0.01;

            var tolerance = new Tolerance(precursorTolerancePpm);
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);

            const string resultFilePath = @"D:\Research\Data\UW\Fusion\WT_D_DDA_130412065618_10ppm_TI2_SGD_Decoy.tsv";

            Console.WriteLine("IsDecoy\tPeptide\tScanNum\tCharge\tSpecEValue\tQValue\tPrecursorMz" +
                              "\tTheo0\tTheo1\tTheo2\tTheo3" +
                              "\tObs0\tCorr0\tObs1\tCorr1\tObs2\tCorr2\tObs3\tCorr3\tObs-1\tCorr-1\tObs0.5\tCorr0.5");
            foreach (var line in File.ReadLines(resultFilePath))
            {
                if (line.StartsWith("#")) continue;
                var token = line.Split('\t');
                if (token.Length != 16) continue;

                var qValue = Convert.ToDouble(token[14]);
                if (qValue > fdrThreshold) continue;

                var peptide = token[8].Replace("C+57.021", "C");
                var scanNum = Convert.ToInt32(token[2]);
                var charge = Convert.ToInt32(token[7]);
                var specEValue = Convert.ToDouble(token[12]);

                var protein = token[9];
                var isDecoy = protein.StartsWith("XXX_");

                var precursorIon = new Ion(aaSet.GetComposition(peptide) + Composition.H2O, charge);
                var baseXic = run.GetExtractedIonChromatogram(precursorIon.GetBaseIsotopeMz(), tolerance, scanNum);
                var baseIntensity = baseXic.GetSumIntensities();

                Console.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", (isDecoy ? 1 : 0), peptide, scanNum, charge, specEValue, qValue, precursorIon.GetMz());

                var isotopeIndices = new double[] {0, 1, 2, 3, -1, 0.5};
                var theoIsotopes = precursorIon.GetIsotopes(0.01);
                var numIsotopes = 0;
                foreach (var theoIsotope in theoIsotopes)
                {
                   Console.Write("\t"+theoIsotope.Item2);
                    if (++numIsotopes == 4) break;
                }

                foreach (var isotopeIndex in isotopeIndices)
                {
                    var isotopeMz = precursorIon.GetIsotopeMz(isotopeIndex);
                    var xic = run.GetExtractedIonChromatogram(isotopeMz, tolerance, scanNum);
                    var relativeIntensity = xic.GetSumIntensities() / baseIntensity;
                    var correlation = xic.GetCorrelation(baseXic);
                    Console.Write("\t{0}\t{1}", relativeIntensity, correlation);
                }
                Console.WriteLine();
            }
        }
    }
}
