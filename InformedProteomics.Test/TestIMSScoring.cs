using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.IMS;
using InformedProteomics.Backend.IMSScoring;
using InformedProteomics.Backend.IMSTraining;
using InformedProteomics.Backend.Scoring;
using NUnit.Framework;
using UIMFLibrary;
using IonType = InformedProteomics.Backend.Data.Spectrometry.IonType;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestImsScoring
    {
        [Test]
        public void TestTrainingUsingMgfFile()
        {
            const string mgf = @"C:\Users\kwj\Dropbox\Training\HCD_train.mgf";
            var tolerance = new Tolerance(30, DataReader.ToleranceType.PPM);
            var spectra = MgfParser.Parse(mgf);
            TrainerUsingMgfFile.Train(@"..\..\..\TestFiles\HCD_train.mgf_para.txt", spectra, tolerance, 3);
        }

        [Test]
        public void TestMsMsSpectrum()
        {
            const string mgf = @"C:\Users\kwj\Dropbox\Training\HCD_train.mgf";
            var spectra = MgfParser.Parse(mgf, 5);
            var tolerance = new Tolerance(20, DataReader.ToleranceType.PPM);
            foreach (var spectrum in spectra)
            {
                var annotation = spectrum.Annotation;
                var ionTypes = new List<IonType>(new IonTypeFactory().GetAllKnownIonTypes());
                for (var k = 1; k < annotation.Count; k++)
                {
                    var peaks = spectrum.GetExplainedPeaks(annotation, k, ionTypes, tolerance);
                    for (var i = 0; i < ionTypes.Count; i++)
                    {
                        if (peaks[i].Intensity > 0)
                            Console.WriteLine(ionTypes[i] + "\t" + peaks[i].Mz + "\t" + peaks[i].Intensity);
                    }
                }
            }

        }

        [Test]
        public void TestScoring()
        {
            const string uimfFilePath =
                @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            const string paramFile = @"..\..\..\TestFiles\HCD_train.mgf_para.txt";
            var imsData = new ImsDataCached(uimfFilePath);
            var imsScorerFactory = new ImsScorerFactory(paramFile);
            
            const string fasta = @"..\..\..\TestFiles\BSA.fasta";
            //const string fasta = @"..\..\..\TestFiles\CCAADDKEACFAVEGPK.fasta";
            var targetDist = new int[1000];
            var decoyDist = new int[1000];
            var targetScores = new Dictionary<double, List<string>>();
            var decoyScores = new Dictionary<double, List<string>>();

            for (var q = 0; q < 2; q++)
            {
                var dist = targetDist;
                var num = 0;
                var pepIndex = 0;
                if (q != 0) dist = decoyDist;
                foreach (var targetPeptide in Misc.GetPeptidesFromFasta(fasta, true, 2, q != 0))
                    // stupid function made by kyowon.
                {
                    Console.WriteLine("{0}: {1}", ++pepIndex, targetPeptide);
                    var pep = targetPeptide;// CACSRKNQVK"GNYKNAYYLLEPAYFYPHR";// "CCAADDKEACFAVEGPK"//targetPeptide; "QLSACKLRQK";
                    //if (pep.Length < 5) continue;
                    var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
                    var precursorComposition = aaSet.GetComposition(pep);
                    var sequence = new Sequence(precursorComposition + Composition.H2O, pep, aaSet);
                    var maxScore = double.NegativeInfinity;
                    for (var charge = 2; charge <= 3; charge++)
                    {
                        var precursorIon = new Ion(precursorComposition + Composition.H2O, charge);
                        var imsScorer = imsScorerFactory.GetImsScorer(imsData, precursorIon);

                        for (var i = 0; i <= 0; i++)
                        {
                            var precursorMz = precursorIon.GetIsotopeMz(i + precursorIon.Composition.GetMostAbundantIsotopeZeroBasedIndex());
                            if (precursorMz > imsData.MaxPrecursorMz || precursorMz < imsData.MinPrecursorMz) continue;
                            var precursorFeatures = imsData.GetPrecursorFeatures(precursorMz);
                            //Console.WriteLine("Precursor: {0}, Charge: {1}\n", precursorMz, charge + "\t" + precursorComposition);
                            foreach (var precursorFeature in precursorFeatures)
                            {
                                //Console.WriteLine("Precursor Feature: " + precursorFeature + "\n");
                                var score = imsScorer.GetPrecursorScore(precursorFeature);
                                //Console.WriteLine("Feature: " + precursorFeature);
                                //Console.WriteLine("Precursor score: " + score);
                                for (var cutNumber = 1; cutNumber < pep.Length; cutNumber++)
                                {
                                    var cutScore = imsScorer.GetCutScore(pep[cutNumber - 1], pep[cutNumber], sequence.GetComposition(0, cutNumber), precursorFeature);
                                    //Console.WriteLine("{0} {1} {2} {3}", pep[cutNumber-1], pep[cutNumber], sequence.GetComposition(0, cutNumber), cutScore);
                                    //var cutNodeScore = imsScorer.GetCutNodeScore(pep[cutNumber - 1], pep[cutNumber], sequence.GetComposition(0, cutNumber), precursorFeature);
                                    //var cutRatioScore = imsScorer.GetCutRatioScore(pep[cutNumber - 1], pep[cutNumber], sequence.GetComposition(0, cutNumber), precursorFeature);
                                    //var cutLCScore = imsScorer.GetCutLCScore(pep[cutNumber - 1], pep[cutNumber], sequence.GetComposition(0, cutNumber), precursorFeature);
                                    //var cutIMSScore = imsScorer.GetCutIMSScore(pep[cutNumber - 1], pep[cutNumber], sequence.GetComposition(0, cutNumber), precursorFeature);
                                    //Console.WriteLine(cutNumber + "\t" + cutNodeScore + "\t" + cutRatioScore + "\t" + cutLCScore + "\t" + cutIMSScore + "\t" +  cutScore);
                                    score += cutScore;
                                }
                                maxScore = Math.Max(maxScore, score);
                                //Console.WriteLine("Score = " + score + "\n");
                            }
                        }
                    }
                    if (q == 0)
                    {
                        if(!targetScores.ContainsKey(maxScore)) targetScores[maxScore] = new List<string>();
                        targetScores[maxScore].Add(pep);
                    }
                    else
                    {
                        if (!decoyScores.ContainsKey(maxScore)) decoyScores[maxScore] = new List<string>();
                        decoyScores[maxScore].Add(pep);
                    }
                    if (maxScore > 0)
                    {
                        Console.WriteLine((q==0? "T" : "D") + " " + num++ + " Max Score of this peptide " + pep + " is " + maxScore);
                    }
                    //break;
                    var scoreIndex = (int)Math.Min(targetDist.Length - 1, Math.Max(0, maxScore + 50));
                    dist[scoreIndex] = dist[scoreIndex] + 1;

                }
            }

            var decoyScoreList = decoyScores.Keys.ToList();

            Console.WriteLine();
            decoyScoreList.Sort();
            var threshold = decoyScoreList[decoyScores.Count - 1];
            var numTarget = 0;
            foreach (var score in targetScores.Keys)
            {
                if (score > threshold)
                {
                    foreach (var peptide in targetScores[score])
                    {
                        Console.WriteLine(peptide);
                    }
                    numTarget+=targetScores[score].Count;
                }
            }
            Console.WriteLine("Threshold is "+threshold + "\nNumber of target is "+numTarget);


           
            var twriter = new StreamWriter(@"..\..\..\TestFiles\target.m");
            var dwriter = new StreamWriter(@"..\..\..\TestFiles\decoy.m");
            
            twriter.Write("t=[");
            for (var j = 0; j < targetDist.Length;j++ )
            {
                twriter.WriteLine(j+"\t"+targetDist[j]);
            }
            twriter.WriteLine("];");
            twriter.Close();

            dwriter.Write("d=[");
            for (var j = 0; j < decoyDist.Length; j++)
            {
                dwriter.WriteLine(j + "\t" + decoyDist[j]);
            }
            dwriter.WriteLine("];");
            dwriter.Close();
        }
    }
}
