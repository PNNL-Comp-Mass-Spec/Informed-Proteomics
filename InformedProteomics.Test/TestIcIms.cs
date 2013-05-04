using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.IMS;
using NUnit.Framework;
using UIMFLibrary;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestIcIms
    {
        [Test]
        public void TestImsFeatureFinding()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            var imsData = new ImsDataCached(uimfFilePath);

            //Console.WriteLine("Generating precursor features (MinMz: " + imsData.MinPrecursorMz + " MaxMz: " + imsData.MaxPrecursorMz + ")");
            //int numPrecursorFeatures = imsData.CreatePrecursorFeatures();
            //Console.WriteLine("TotalNumPrecursorFeatures: " + numPrecursorFeatures);            

            const string targetPeptide = "CCHGDLLECADDRADLAK";//EYANQFMWEYSTNYGQAPLSLLVSYTK CCAADDKEACFAVEGPK
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

            const string targetFragment = "YGQAPLSLLVSYTK";//"AADDKEACFAVEGPK";
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
    }
}
