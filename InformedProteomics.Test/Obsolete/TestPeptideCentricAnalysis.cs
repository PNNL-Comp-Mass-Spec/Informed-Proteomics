using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.DIA.Search;
using NUnit.Framework;

namespace InformedProteomics.Test.Obsolete
{
    [TestFixture]
    internal class TestPeptideCentricAnalysis
    {
        [Test]
        public void CompareRt()
        {
            // Q-Exactive
            //const string qeDdaResult = @"D:\Research\Data\UW\QExactive\DDA_All_Summary.tsv";
            //const string qeDiaResult = @"D:\Research\Data\UW\QExactive\DIA_All_Summary.tsv";

            // Fusion
            const string qeDdaResult = @"D:\Research\Data\UW\Fusion\DDA_Summary.tsv";
            const string qeDiaResult = @"D:\Research\Data\UW\Fusion\DIA_Summary.tsv";

            const string specFileDda = @"D:\Research\Data\UW\QExactive\82593_lv_mcx_DDA.raw";
            var ddaReader = new XCaliburReader(specFileDda);

            var specFileToReader = new Dictionary<string, XCaliburReader>();
            var specFilesDia = Directory.GetFiles(@"D:\Research\Data\UW\QExactive\", "*_DIA_*.raw");
            foreach (var specFile in specFilesDia)
            {
                var specFileNoExt = Path.GetFileNameWithoutExtension(specFile);
                if (specFileNoExt == null) continue;
                var reader = new XCaliburReader(specFile);
                specFileToReader.Add(specFileNoExt, reader);
            }

            const string resultPath1 = qeDdaResult;
            const string resultPath2 = qeDiaResult;

            var result1 = new TsvFileParser(resultPath1);
            var result2 = new TsvFileParser(resultPath2);

            const double pepQValueThreshold = 0.01;
            var vennDiagram = new VennDiagram<string>(result1.GetPeptides(pepQValueThreshold),
                                                      result2.GetPeptides(pepQValueThreshold));

            var intersectionPeptides = vennDiagram.Intersection;

            var result1Peptides = result1.GetData("Peptide");
            var result1ScanNums = result1.GetData("ScanNum");

            var result2Peptides = result2.GetData("Peptide");
            var result2ScanNums = result2.GetData("ScanNum");
            var result2SpecFile = result2.GetData("#SpecFile");

            Console.WriteLine("Peptide\tScanNum1\tScanNum2\tRt1\tRt2");
            foreach (var peptide in intersectionPeptides)
            {
                var index1 = result1Peptides.IndexOf(peptide);
                var index2 = result2Peptides.IndexOf(peptide);

                var scanNum1 = Convert.ToInt32(result1ScanNums[index1]);
                var scanNum2 = Convert.ToInt32(result2ScanNums[index2]);

                var diaFile = Path.GetFileNameWithoutExtension(result2SpecFile[index2]);

                var reader1 = ddaReader;
                var reader2 = specFileToReader[diaFile];

                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", peptide.Replace("C+57.021", "C"), scanNum1, scanNum2, reader1.RtFromScanNum(scanNum1), reader2.RtFromScanNum(scanNum2));
            }
        }

        [Test]
        public void CompareRtFusion()
        {
            // Fusion
            const string qeDdaResult = @"D:\Research\Data\UW\Fusion\DDA_Summary.tsv";
            const string qeDiaResult = @"D:\Research\Data\UW\Fusion\DIA_Summary.tsv";

            const string specFileDda = @"D:\Research\Data\UW\Fusion\WT_D_DDA_130412065618.raw";
            var ddaReader = new XCaliburReader(specFileDda);

            const string specFileDia = @"D:\Research\Data\UW\Fusion\WT_D_DIA_130412091220.raw";
            var diaReader = new XCaliburReader(specFileDia);

            const string resultPath1 = qeDdaResult;
            const string resultPath2 = qeDiaResult;

            var result1 = new TsvFileParser(resultPath1);
            var result2 = new TsvFileParser(resultPath2);

            const double pepQValueThreshold = 0.01;
            var vennDiagram = new VennDiagram<string>(result1.GetPeptides(pepQValueThreshold),
                                                      result2.GetPeptides(pepQValueThreshold));

            var intersectionPeptides = vennDiagram.Intersection;

            var result1Peptides = result1.GetData("Peptide");
            var result1ScanNums = result1.GetData("ScanNum");

            var result2Peptides = result2.GetData("Peptide");
            var result2ScanNums = result2.GetData("ScanNum");

            Console.WriteLine("Peptide\tScanNum1\tScanNum2\tRt1\tRt2");
            foreach (var peptide in intersectionPeptides)
            {
                var index1 = result1Peptides.IndexOf(peptide);
                var index2 = result2Peptides.IndexOf(peptide);

                var scanNum1 = Convert.ToInt32(result1ScanNums[index1]);
                var scanNum2 = Convert.ToInt32(result2ScanNums[index2]);

                var reader1 = ddaReader;
                var reader2 = diaReader;

                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", peptide.Replace("C+57.021", "C"), scanNum1, scanNum2, reader1.RtFromScanNum(scanNum1), reader2.RtFromScanNum(scanNum2));
            }
        }

        [Test]
        public void RunPeptideCentricAnalysis()
        {
            //TestQExactiveDdaDataPostProcessing();
            //TestQExactiveDiaDataPostProcessing();
            //TestQExactiveDiaDataPostProcessingNoEdgeNtt2();
            //TestFusionDiaDataPostProcessing();
            //TestFusionDdaDataPostProcessing();
            TestQExactiveDdaDataPostProcessingPerFile();
        }

        [Test]
        public void TestQExactiveDdaDataPostProcessing()
        {
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
            const string resultPath = @"D:\Research\Data\UW\QExactive\DIA_Results\DIA_All.tsv";
            var specFiles = Directory.GetFiles(@"D:\Research\Data\UW\QExactive\", "*_DIA_*.raw");
            const string outputFilePath = @"D:\Research\Data\UW\QExactive\DIA_All_Summary.tsv";

            var postProcessor = new MsGfPostProcessor(specFiles, resultPath, new Tolerance(20), new Tolerance(10));
            var numId = postProcessor.PostProcessing(outputFilePath); 

            Console.WriteLine("NumId: {0}", numId);
        }

        [Test]
        public void TestQExactiveDiaDataPostProcessingNoEdgeNtt2()
        {
            const string resultPath = @"D:\Research\Data\UW\QExactive\NoEdge\DIA_All.tsv";
            var specFiles = Directory.GetFiles(@"D:\Research\Data\UW\QExactive\", "*_DIA_*.raw");
            const string outputFilePath = @"D:\Research\Data\UW\QExactive\DIA_All_Edge_Summary.tsv";

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
        public void TestQExactiveDdaDataPostProcessingPerFile()
        {
            var ranges = new[] { "400to525", "525to650", "650to775", "775to900" };
            foreach (var range in ranges)
            {
                Console.WriteLine("Range: {0}", range);
                TestQExactiveDdaDataPostProcessingPerFile(range);
            }
        }

        private void TestQExactiveDdaDataPostProcessingPerFile(string range)
        {
            var resultPath = @"D:\Research\Data\UW\QExactive\DIA_Results\DIA_" + range + ".tsv";
            var specFilePath = @"D:\Research\Data\UW\QExactive\82593_lv_mcx_DIA_5mz_" + range + ".raw";
            var outputFilePath = @"D:\Research\Data\UW\QExactive\DIA_" + range + "_Summary.tsv";

            var postProcessor = new MsGfPostProcessor(specFilePath, resultPath, new Tolerance(20), new Tolerance(10));
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
            const string qeDdaMsGfResult = @"D:\Research\Data\UW\QExactive\82593_lv_mcx_DDA.tsv";
            const string qeDdaResult = @"D:\Research\Data\UW\QExactive\DDA_All_Summary.tsv";
            const string qeDiaResult = @"D:\Research\Data\UW\QExactive\DIA_All_Summary.tsv";

            // PE-MMR
            const string pemmrQeDia = @"D:\Research\Data\UW\PEMMR\QE\allFDR.tsv";
            const string pemmrFusionDda = @"D:\Research\Data\UW\PEMMR\Fusion\WT_D_DDA_130412065618_MX_PEMMR.tsv";
            const string pemmrFusionDia = @"D:\Research\Data\UW\PEMMR\Fusion\WT_D_DIA_130412091220_MX_PEMMR.tsv";

            const string resultPath1 = pemmrFusionDda;
            const string resultPath2 = fusionDdaResult;

            var result1 = new TsvFileParser(resultPath1);
            var result2 = new TsvFileParser(resultPath2);

            const double pepQValueThreshold = 0.01;
            var vennDiagram = new VennDiagram<string>(result1.GetPeptides(pepQValueThreshold),
                                                      result2.GetPeptides(pepQValueThreshold));
            Console.WriteLine("{0}\t{1}\t{2}",
                              vennDiagram.Set1Only.Count + vennDiagram.Intersection.Count,
                              vennDiagram.Intersection.Count,
                              vennDiagram.Set2Only.Count + vennDiagram.Intersection.Count);
        }

        [Test]
        public void ProcessPemmrData()
        {
            const string resultPath = @"D:\Research\Data\PEMMR\Ox\iTRAQ_N33T34_10ug_100cm_300min_C2_061213_All.tsv";
            const string outputFilePath = @"D:\Research\Data\PEMMR\Ox\IPA_Summary.tsv";
            const string specFilePath = @"D:\Research\Data\PEMMR\Spectra\iTRAQ_N33T34_10ug_100cm_300min_C2_061213.raw"; ;

            var postProcessor = new MsGfPostProcessor(specFilePath, resultPath, new Tolerance(5), new Tolerance(5));
            var numId = postProcessor.PostProcessing(outputFilePath);

            Console.WriteLine("NumId: {0}", numId);
        }

        [Test]
        public void GenerateVennDiagramsPeMmr()
        {
            // No PE-MMR
            const string noPeMmr = @"D:\Research\Data\PEMMR\iTRAQ_N33T34_10ug_100cm_300min_C2_061213.tsv";

            // PE-MMR Scan based FDR
            const string scanBasedPeMmr = @"D:\Research\Data\PEMMR\NewSpectra\iTRAQ_N33T34_10ug_100cm_300min_C2_061213_MX_PEMMR_UMCID_ScanFDR.tsv";

            // UMC based FDR
            const string umcBasedPeMmr = @"D:\Research\Data\PEMMR\NewSpectra\iTRAQ_N33T34_10ug_100cm_300min_C2_061213_MX_PEMMR_UMCID_UMCFDR.tsv";

            // IPA
            const string ipa = @"D:\Research\Data\PEMMR\Ox\IPA_Summary_TargetOnly.tsv";

            const string resultPath1 = umcBasedPeMmr;
            const string resultPath2 = ipa;

            var result1 = new TsvFileParser(resultPath1);
            var result2 = new TsvFileParser(resultPath2);

            const double pepQValueThreshold = 0.01;
            var vennDiagram = new VennDiagram<string>(result1.GetPeptides(pepQValueThreshold),
                                                      result2.GetPeptides(pepQValueThreshold));
            Console.WriteLine("{0}\t{1}\t{2}",
                              vennDiagram.Set1Only.Count + vennDiagram.Intersection.Count,
                              vennDiagram.Intersection.Count,
                              vennDiagram.Set2Only.Count + vennDiagram.Intersection.Count);
            Console.WriteLine("{0}\t{1}\t{2}",
                              vennDiagram.Set1Only.Count,
                              vennDiagram.Intersection.Count,
                              vennDiagram.Set2Only.Count);

            foreach (var peptide in vennDiagram.Set2Only)
            {
                Console.WriteLine(peptide);
                var peptides = result2.GetData("Peptide");
            }
        }

        [Test]
        public void ProcessMhcData()
        {
            const string resultPath = @"D:\Research\Data\ImmunoPeptidomics\Benchmarking\IPA\carone_C1309_All.tsv";
            const string outputFilePath = @"D:\Research\Data\ImmunoPeptidomics\Benchmarking\IPA\IPA_Summary.tsv";
            var specFiles = Directory.GetFiles(@"D:\Research\Data\ImmunoPeptidomics\Benchmarking\raw", "*.raw");

            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);

            var searchModifications = new List<SearchModification>
            {
                oxM
            };
            var aaSet = new AminoAcidSet(searchModifications, 2);

            var postProcessor = new MsGfPostProcessor(specFiles, resultPath, new Tolerance(5), new Tolerance(3));
            var numId = postProcessor.PostProcessing(outputFilePath);

            Console.WriteLine("NumId: {0}", numId);
        }

        [Test]
        public void TestSpecEValueCalibration()
        {
            const int targetCharge = 4;
            //const string resultFilePath = @"D:\Research\Data\UW\QExactive\DIA_Results\DIA_All.tsv";
            const string resultFilePath = @"D:\Research\Data\UW\QExactive\DIA_All_Summary - Copy.tsv";
            var histSpecEValue = new float[2,500];
            const int specEValueBinningConstant = 10;
            
            var histMsGfScore = new float[2,1000];
            const int minMsGfScore = -100;
            var numOccs = new[] {0, 0};
            MsGfPlusHeaderInformation headerInfo = null;
            var prevScanNum = new[] {-1, -1};
            
            foreach (var line in File.ReadLines(resultFilePath))
            {
                if (line.StartsWith("#"))
                {
                    if(headerInfo == null) headerInfo = new MsGfPlusHeaderInformation(line);
                    continue;
                }

                if (headerInfo == null) continue;

                var token = line.Split('\t');
                if (token.Length != headerInfo.NumColumns) continue;

                var protein = token[headerInfo.ProteinColNum];
                var isDecoy = protein.StartsWith("XXX");

                var targetIndex = isDecoy ? 1 : 0;

                var scanNum = Convert.ToInt32(token[headerInfo.ScanNumColNum]);
                if (scanNum == prevScanNum[targetIndex]) continue;
                prevScanNum[targetIndex] = scanNum;

                var specFile = token[headerInfo.SpecFileColNum];
                var precursorMz = Convert.ToDouble(token[headerInfo.PrecursorColNum]);
                var charge = Convert.ToInt32(token[headerInfo.ChargeColNum]);

                //if (charge == 1) continue;
                if (targetCharge > 0 && charge != targetCharge) continue;
                var specEValue = Convert.ToDouble(token[headerInfo.SpecEValueColNum]);
                var specEValueBin = (int)Math.Round(-Math.Log10(specEValue) * specEValueBinningConstant);
                histSpecEValue[targetIndex,specEValueBin]++;

                var msgfScore = Convert.ToInt32(token[headerInfo.MsgfScoreColNum]);
                histMsGfScore[targetIndex, msgfScore - minMsGfScore]++;

                numOccs[targetIndex]++;
            }

            for (var i = 0; i < histMsGfScore.GetLength(0); i++)
            {
                for (var j = 0; j < histMsGfScore.GetLength(1); j++)
                {
                    histMsGfScore[i,j] /= numOccs[i];
                }
            }

            var specEValueArr = new float[histSpecEValue.GetLength(1)];
            for (var i = 0; i < histSpecEValue.GetLength(0); i++)
            {
                for (var j = 0; j < histSpecEValue.GetLength(1); j++)
                {
                    if(i==0) specEValueArr[j] = (float)j / specEValueBinningConstant;
                    histSpecEValue[i, j] /= numOccs[i];
                }
            }
            
            Console.WriteLine("specEValueX <- c({0});", string.Join(",", specEValueArr));
            Console.WriteLine("specEValueHistTarget <- c({0});", string.Join(",", Slice(histSpecEValue, 0)));
            Console.WriteLine("specEValueHistDecoy <- c({0});", string.Join(",", Slice(histSpecEValue, 1)));
            Console.WriteLine("#NumTargetOccs: {0}", numOccs[0]);
            Console.WriteLine("#NumDecoyOccs: {0}", numOccs[1]);
        }

        public static IEnumerable<float> Slice(float[,] histArray, int targetIndex)
        {
            for (var i = 0; i < histArray.GetLength(1); i++)
            {
                yield return histArray[targetIndex, i];
            }
        }

        [Test]
        public void TestXicGen()
        {
            const string specFilePath = @"D:\Research\Data\UW\Fusion\WT_D_DDA_130412065618.raw";
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);

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
            foreach(var p in precursorIon.GetIsotopes(0.1)) Console.WriteLine("{0}\t{1}", precursorIon.GetIsotopeMz(p.Index), p.Ratio);

            var xicArr = new Dictionary<int, Xic>();
            var basePeakIndex = precursorIon.Composition.GetMostAbundantIsotopeZeroBasedIndex();
            for (var i = -1; i < 3; i++)
            {
                xicArr[i] = run.GetPrecursorExtractedIonChromatogram(precursorIon.GetIsotopeMz(i), tolerance, targetScanNum);
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
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);

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
                var baseXic = run.GetPrecursorExtractedIonChromatogram(precursorIon.GetMostAbundantIsotopeMz(), tolerance, scanNum);
                var baseIntensity = baseXic.GetSumIntensities();

                var isValid = true;
                foreach (var isotope in precursorIon.GetIsotopes(relativeIntensityThreshold))
                {
                    if (isotope.Index == basePeakIndex) continue;
                    var isotopeMz = precursorIon.GetIsotopeMz(isotope.Index);
                    var xic = run.GetPrecursorExtractedIonChromatogram(isotopeMz, tolerance, scanNum);

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
                    //    //Console.WriteLine("Off ratio\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", isDecoy, peptide, scanNum, charge, precursorIon.GetMonoIsotopicMz(), isotopeMz, isotopeRatio);
                    //    break;
                    //}

                    // Check if correlation is high
                    //if (correlation < correlationThreshold)
                    //{
                    //    isValid = false;
                    //    //Console.WriteLine("Low correlation\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", isDecoy, peptide, scanNum, charge, precursorIon.GetMonoIsotopicMz(), isotopeMz, correlation);
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
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun);
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
                var baseXic = run.GetPrecursorExtractedIonChromatogram(precursorIon.GetMostAbundantIsotopeMz(), tolerance, scanNum);
                var baseIntensity = baseXic.GetSumIntensities();

                Console.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", (isDecoy ? 1 : 0), peptide, scanNum, charge, specEValue, qValue, precursorIon.GetMonoIsotopicMz());

                var isotopeIndices = new double[] {0, 1, 2, 3, -1, 0.5};
                var theoIsotopes = precursorIon.GetIsotopes(0.01);
                var numIsotopes = 0;
                foreach (var theoIsotope in theoIsotopes)
                {
                   Console.Write("\t"+theoIsotope.Ratio);
                    if (++numIsotopes == 4) break;
                }

                foreach (var isotopeIndex in isotopeIndices)
                {
                    var isotopeMz = precursorIon.GetIsotopeMz(isotopeIndex);
                    var xic = run.GetPrecursorExtractedIonChromatogram(isotopeMz, tolerance, scanNum);
                    var relativeIntensity = xic.GetSumIntensities() / baseIntensity;
                    var correlation = xic.GetCorrelation(baseXic);
                    Console.Write("\t{0}\t{1}", relativeIntensity, correlation);
                }
                Console.WriteLine();
            }
        }
    }
}
