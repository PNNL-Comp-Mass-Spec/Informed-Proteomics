using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.IMS;
using InformedProteomics.Backend.IMSScoring;
using NUnit.Framework;
using UIMFLibrary;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestIcIms
    {
        [Test]
        public void TestSarcSearch()
        {
            const string uimfFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo_Collision_Energy_Collapsed.UIMF";

            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);

            var imsData = new ImsDataCached(uimfFilePath);

            const string paramFile = @"..\..\..\TestFiles\HCD_train.mgf_para.txt";
            var imsScorerFactory = new ImsScorerFactory(paramFile);

            Console.WriteLine("Finished reading the UIMF file.");

            var dbFilePaths = new string[] 
                {
                    @"C:\cygwin\home\kims336\Data\IMS_Sarc\HumanPeptides.txt",
                    @"C:\cygwin\home\kims336\Data\IMS_Sarc\HumanPeptides_Reverse.txt"
                };

            var isDecoy = 0;

            using (var writer = new StreamWriter(@"C:\cygwin\home\kims336\Data\IMS_Sarc\ICResults.txt"))
            {
                writer.WriteLine(
                    "Peptide\tPrecursorMz\tPrecursorCharge\tLCBegin\tLCEnd\tIMSBegin\tIMSEnd\tApexLC\tApexIMS\tSumIntensities\tNumPoints\tScore\tIsDecoy");
                foreach (var dbFilePath in dbFilePaths)
                {
                    var lineNum = 0;
                    foreach (var annotation in File.ReadLines(dbFilePath))
                    {
                        ++lineNum;
                        if (lineNum%10000 == 1)
                        {
                            Console.WriteLine("Processing {0} peptides...", lineNum);
                        }
                        var peptide = annotation.Substring(2, annotation.Length - 4);

                        var seqGraph = new SequenceGraph(aaSet, peptide);
                        foreach (var scoringGraph in seqGraph.GetScoringGraphs())
                        {
                            scoringGraph.RegisterImsData(imsData, imsScorerFactory);
                            for (var precursorCharge = 1; precursorCharge <= 5; precursorCharge++)
                            {
                                double precursorMz = scoringGraph.GetPrecursorIon(precursorCharge).GetMz();
                                if (precursorMz > imsData.MaxPrecursorMz || precursorMz < imsData.MinPrecursorMz)
                                    continue;
                                var best = scoringGraph.GetBestFeatureAndScore(precursorCharge);
                                var feature = best.Item1;
                                if (best.Item1 != null && best.Item2 > 0)
                                {
                                    writer.WriteLine(
                                        "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}"
                                        , annotation
                                        , precursorMz
                                        , precursorCharge
                                        , feature.ScanLcStart
                                        , feature.ScanLcStart + feature.ScanLcLength - 1
                                        , feature.ScanImsStart
                                        , feature.ScanImsStart + feature.ScanImsLength - 1
                                        , feature.ScanLcStart + feature.ScanLcRepOffset
                                        , feature.ScanImsStart + feature.ScanImsRepOffset
                                        , feature.SumIntensities
                                        , feature.NumPoints
                                        , best.Item2
                                        , isDecoy
                                        );
                                }
                            }
                        }
                    }
                    ++isDecoy;
                }
            }
        }

        [Test]
        public void TestBsaSearch()
        {
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var fixCarbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C', SequenceLocation.Everywhere, true);

            var searchModifications = new List<SearchModification> { fixCarbamidomethylC };
            const int numMaxModsPepPeptide = 2;
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPepPeptide);

            const string uimfFilePath =
                @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            var imsData = new ImsDataCached(uimfFilePath);

            const string paramFile = @"..\..\..\TestFiles\HCD_train.mgf_para.txt";
            var imsScorerFactory = new ImsScorerFactory(paramFile);


            var dbFilePaths = new string[] 
                {
                    @"..\..\..\TestFiles\BSAPeptides_ST.txt",
                    @"..\..\..\TestFiles\BSA_ReversePeptides_ST.txt"
                };

            var isDecoy = 0;

            //using (var writer = new StreamWriter(@"C:\cygwin\home\kims336\Data\IMS_Sarc\ICResults.txt"))
            using (var writer = new StreamWriter(@"..\..\..\TestFiles\IC_BSA_Results.txt"))
            {
                writer.WriteLine(
                    "Peptide\tPrecursorMz\tPrecursorCharge\tLCBegin\tLCEnd\tIMSBegin\tIMSEnd\tApexLC\tApexIMS\tSumIntensities\tNumPoints\tScore\tIsDecoy");
                foreach (var dbFilePath in dbFilePaths)
                {
                    foreach (var annotation in File.ReadLines(dbFilePath))
                    {
                        var peptide = annotation.Substring(2, annotation.Length - 4);

                        //var scoringGraph = seqGraph.GetScoringGraph(0);
                        var seqGraph = new SequenceGraph(aaSet, peptide);
                        foreach (var scoringGraph in seqGraph.GetScoringGraphs())
                        {
                            scoringGraph.RegisterImsData(imsData, imsScorerFactory);
                            for (var precursorCharge = 1; precursorCharge <= 5; precursorCharge++)
                            {
                                var precursorMz = scoringGraph.GetPrecursorIon(precursorCharge).GetMz();
                                if (precursorMz > imsData.MaxPrecursorMz || precursorMz < imsData.MinPrecursorMz)
                                    continue;
                                var best = scoringGraph.GetBestFeatureAndScore(precursorCharge);
                                var feature = best.Item1;
                                if (best.Item1 != null && best.Item2 > 0)
                                {
                                    writer.WriteLine(
                                        "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}"
                                        , annotation
                                        , precursorMz
                                        , precursorCharge
                                        , feature.ScanLcStart
                                        , feature.ScanLcStart + feature.ScanLcLength - 1
                                        , feature.ScanImsStart
                                        , feature.ScanImsStart + feature.ScanImsLength - 1
                                        , feature.ScanLcStart + feature.ScanLcRepOffset
                                        , feature.ScanImsStart + feature.ScanImsRepOffset
                                        , feature.SumIntensities
                                        , feature.NumPoints
                                        , best.Item2
                                        , isDecoy
                                        );
                                }
                            }
                        }
                    }
                    ++isDecoy;
                }
            }
        }

        [Test]
        public void TestDbSearch()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            var imsData = new ImsDataCached(uimfFilePath);

            const string paramFile = @"..\..\..\TestFiles\HCD_train.mgf_para.txt";
            var imsScorerFactory = new ImsScorerFactory(paramFile);
        }

        [Test]
        public void TestImsScoring()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            //var imsData = new ImsDataCached(uimfFilePath);
            var imsData = new ImsDataCached(uimfFilePath, 400.0, 2000.0, 10.0, 2500.0,
                new Tolerance(25, DataReader.ToleranceType.PPM), new Tolerance(25, DataReader.ToleranceType.PPM));
            const string paramFile = @"..\..\..\TestFiles\HCD_train.mgf_para.txt";
            var imsScorerFactory = new ImsScorerFactory(paramFile);

            //const string targetPeptide = "CCAADDKEACFAVEGPK";
            //const string targetPeptide = "ECCHGDLLECADDRADLAK";
            const string targetPeptide = "CHGDLLECADDRADLAK";
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);

            var seqGraph = new SequenceGraph(aaSet, targetPeptide);
            var scoringGraph = seqGraph.GetScoringGraph(0);
            //foreach (var composition in scoringGraph.GetCompositions())
            //{
            //    Console.WriteLine(composition);
            //}
            scoringGraph.RegisterImsData(imsData, imsScorerFactory);
            for (var precursorCharge = 3; precursorCharge <= 3; precursorCharge++)
            {
                var best = scoringGraph.GetBestFeatureAndScore(precursorCharge);
                Console.WriteLine("PrecursorMz: " + scoringGraph.GetPrecursorIon(precursorCharge).GetMz());
                Console.WriteLine("Charge: " + precursorCharge);
                Console.WriteLine("Feature: " + best.Item1);
                Console.WriteLine("Score: " + best.Item2);
                Console.WriteLine();
            }
        }

        [Test]
        public void TestImsFeatureFinding()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            var imsData = new ImsDataCached(uimfFilePath);

            //Console.WriteLine("Generating precursor features (MinMz: " + imsData.MinPrecursorMz + " MaxMz: " + imsData.MaxPrecursorMz + ")");
            //int numPrecursorFeatures = imsData.CreatePrecursorFeatures();
            //Console.WriteLine("TotalNumPrecursorFeatures: " + numPrecursorFeatures);            

            const string targetPeptide = "CCAADDKEACFAVEGPK";//EYANQFMWEYSTNYGQAPLSLLVSYTK  CCHGDLLECADDRADLAK
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var precursorComposition = aaSet.GetComposition(targetPeptide);
            for (int charge = 3; charge <= 3; charge++)
            {
                var precursorIon = new Ion(precursorComposition+Composition.H2O, charge);
                double precursorMz = precursorIon.GetMz();
                FeatureSet precursorFeatures = imsData.GetPrecursorFeatures(precursorMz);
                Console.WriteLine("Precursor: {0}, Charge: {1}\n", precursorMz, charge);
                foreach (Feature precursorFeature in precursorFeatures)
                {
                    Console.WriteLine(precursorFeature);
                    Console.WriteLine("LC Apex profile");
                    Console.WriteLine(string.Join("\t", precursorFeature.LcApexPeakProfile));
                    Console.WriteLine("IMS Apex profile");
                    Console.WriteLine(string.Join("\t", precursorFeature.ImsApexPeakProfile));
                    Console.WriteLine();
                }
            }

            //Console.WriteLine("Generating fragment features (MinMz: " + imsData.MinFragmentMz + " MaxMz: " + imsData.MaxFragmentMz + ")");
            //int numFragmentFeatures = imsData.CreateFragmentFeatures();
            //Console.WriteLine("TotalNumFragmentFeatures: " + numFragmentFeatures);

            const string targetFragment = "AADDKEACFAVEGPK"; //"YGQAPLSLLVSYTK";
            var fragmentComposition = aaSet.GetComposition(targetFragment);
                for (int charge = 2; charge <= 2; charge++)
            {
                var y15C2 = new Ion(fragmentComposition + Composition.H2O, 2);
                double fragmentMz = y15C2.GetMz();
                FeatureSet fragmentFeatures = imsData.GetFragmentFeatures(fragmentMz);
                Console.WriteLine("Fragment: {0}, Charge: {1}", fragmentMz, charge);
                foreach (Feature fragmentFeature in fragmentFeatures)
                {
                    Console.WriteLine(fragmentFeature);
                    Console.WriteLine("LC Apex profile");
                    Console.WriteLine(string.Join("\t", fragmentFeature.LcApexPeakProfile));
                    Console.WriteLine("IMS Apex profile");
                    Console.WriteLine(string.Join("\t", fragmentFeature.ImsApexPeakProfile));
                    Console.WriteLine();
                }
            }
        }

        [Test]
        public void TestGeneratingAllPrecursorXicsForEachMzBin()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";

            var imsData = new ImsData(uimfFilePath);
            var precursorXicMap = new Dictionary<int, FeatureSet>();

            int totalNumFeatures = 0;

            const double minMz = 400.0;
            const double maxMz = 2500.0;

            int minTargetBin = imsData.GetBinFromMz(minMz);
            int maxTargetBin = imsData.GetBinFromMz(maxMz);
            if (maxTargetBin >= imsData.GetNumberOfBins())
                maxTargetBin = imsData.GetNumberOfBins()-1;

            //Console.WriteLine("TargetMz: " + theoMz);
            Console.WriteLine("MinMz: " + minMz + " MaxMz: " + maxMz);
            Console.WriteLine("MinTargetBin: " + minTargetBin + " MaxTargetBin: " + maxTargetBin);
            for (int targetBin = minTargetBin; targetBin <= maxTargetBin; targetBin++)
            {
                double mz = imsData.GetMzFromBin(targetBin);
                FeatureSet featureSet = imsData.GetFeatures(targetBin, DataReader.FrameType.MS1);
                int numFeatures = featureSet.GetFeatures().Count();
                Console.WriteLine("Bin: " + targetBin + "\tMZ: " + mz + "\tNumFeatures: " + numFeatures);
                if (featureSet.GetFeatures().Any())
                    precursorXicMap.Add(targetBin, featureSet);
                totalNumFeatures += numFeatures;
            }
            Console.WriteLine("TotalNumPrecursorFeatures: " + totalNumFeatures);


        }

        [Test]
        public void TestGeneratingAllProductXicsForEachMzBin()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";

            var imsData = new ImsData(uimfFilePath);
            var fragmentXicMap = new Dictionary<int, FeatureSet>();

            int totalNumFragmentFeatures = 0;

            const double minMz = 0.0;
            const double maxMz = 2500.0;

            int minFragmentTargetMz = imsData.GetBinFromMz(minMz);
            int maxFragmentTargetBin = imsData.GetBinFromMz(maxMz);
            if (maxFragmentTargetBin >= imsData.GetNumberOfBins())
                maxFragmentTargetBin = imsData.GetNumberOfBins() - 1;

            Console.WriteLine("MinMz: " + minMz + " MaxMz: " + maxMz);
            Console.WriteLine("MinTargetBin: " + minFragmentTargetMz + " MaxTargetBin: " + maxFragmentTargetBin);
            for (int targetBin = minFragmentTargetMz; targetBin <= maxFragmentTargetBin; targetBin++)
            {
                double mz = imsData.GetMzFromBin(targetBin);
                FeatureSet featureSet = imsData.GetFeatures(targetBin, DataReader.FrameType.MS2);
                int numFeatures = featureSet.GetFeatures().Count();
                Console.WriteLine("Bin: " + targetBin + "\tMZ: " + mz + "\tNumFeatures: " + numFeatures);
                if (featureSet.GetFeatures().Any())
                    fragmentXicMap.Add(targetBin, featureSet);
                totalNumFragmentFeatures += numFeatures;
            }
            Console.WriteLine("TotalNumFragmentFeatures: " + totalNumFragmentFeatures);
        }

        [Test]
        public void TestSimple()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            var uimfReader = new DataReader(uimfFilePath);         
            Console.WriteLine("NumFrames: " + uimfReader.GetGlobalParameters().NumFrames);
            Console.WriteLine("NumScans: " + uimfReader.GetFrameParameters(1).Scans);
        }

        [Test]
        public void TestGeneratingSpectrum()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            var dataReader = new DataReader(uimfFilePath);
            const int startFrameNumber = 579;
            const int endFrameNumber = 581;
            const int startScanNumber = 124;
            const int endScanNumber = 135;

            double[] mzArray;
            int[] intensityArray;
            dataReader.GetSpectrum(startFrameNumber, endFrameNumber, DataReader.FrameType.MS2, startScanNumber,
                                   endScanNumber, out mzArray, out intensityArray);

            var size = mzArray.Length;
            using (var writer = new StreamWriter(@"C:\cygwin\home\kims336\Data\IMS_BSA\testSpec.mgf"))
            {
                writer.WriteLine("BEGIN IONS");
                writer.WriteLine("TITLE=test");
                writer.WriteLine("PEPMASS=582.3447");
                writer.WriteLine("CHARGE=2+");
                const int windowSize = 10;
                for (var i = 0; i < size; i++)
                {
                    bool isPeak = true;
                    int rank = 0;
                    for (var j = i - windowSize; j < i + windowSize; j++)
                    {
                        if (j >= 0 && j < intensityArray.Length)
                        {
                            if (intensityArray[i] < intensityArray[j])
                            {
                                isPeak = false;
                                break;
                            }
                            if (mzArray[j] > mzArray[i] - 0.2 && mzArray[j] < mzArray[i] + 0.2)
                            {
                                if (intensityArray[j] <= intensityArray[i]) rank++;
                            }
                        }
                    }
                    //if (rank <= 1)
                    {
                        writer.WriteLine(mzArray[i] + "\t" + intensityArray[i]);
                    }
                }
                writer.WriteLine("END IONS");
            }
        }
    }
}
