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
using Feature = InformedProteomics.Backend.IMS.Feature;
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
            //var targetDist = new int[1000];
            //var decoyDist = new int[1000];

            var highestScorePerFeature = new Dictionary<Feature, Tuple<double, bool>>();

            const string targetTxt = @"..\..\..\TestFiles\BSA_ST.txt";
            const string decoyTxt = @"..\..\..\TestFiles\BSA_ST_Rev.txt";

            var targetMatches = new List<Tuple<double, string, Feature>>(); // score, peptide, feature
            var decoyMatches = new List<Tuple<double, string, Feature>>(); // score, peptide, feature

            for (var q = 0; q < 2; q++)
            {
                //var dist = targetDist;
                var num = 0;
                var pepIndex = 0;
                //if (q != 0) dist = decoyDist;
                var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
                    
                foreach (var targetPeptide in  Misc.GetPeptidesFromTxt(q==0? targetTxt : decoyTxt))//Misc.GetPeptidesFromFasta(fasta, false, 2, q != 0))
                    // stupid function made by kyowon.
                {
                    //Console.WriteLine("{0}: {1}", ++pepIndex, targetPeptide);
                    var pep = targetPeptide;// "CCAADDKEACFAVEGPK";// "LVNELTEFAK";// targetPeptide;// CACSRKNQVK"GNYKNAYYLLEPAYFYPHR";// "CCAADDKEACFAVEGPK"//targetPeptide; "QLSACKLRQK";
                    var precursorComposition = aaSet.GetComposition(pep);
                    var sequence = new Sequence(precursorComposition + Composition.H2O, pep, aaSet);
                    var maxScore = double.NegativeInfinity;
                    var maxPortionOfExplainedFrag = 0.0;

                    Feature maxFeature = null;
                    for (var charge = 1; charge <= 5; charge++)
                    {
                        var precursorIon = new Ion(precursorComposition + Composition.H2O, charge);
                        var imsScorer = imsScorerFactory.GetImsScorer(imsData, precursorIon);
                        var precursorMz = precursorIon.GetIsotopeMz(precursorIon.Composition.GetMostAbundantIsotopeZeroBasedIndex());
                        if (precursorMz > imsData.MaxPrecursorMz || precursorMz < imsData.MinPrecursorMz) continue;
                        var precursorFeatures = imsData.GetPrecursorFeatures(precursorMz);
                        //Console.WriteLine("Precursor: {0}, Charge: {1}\n", precursorMz, charge + "\t" + precursorComposition);
                        foreach (var precursorFeature in precursorFeatures)
                        {
                            //  Console.WriteLine("Precursor Feature: " + precursorFeature + "\n");
                            var score = imsScorer.GetPrecursorScore(precursorFeature);
                            //   Console.WriteLine("Feature: " + precursorFeature);
                            // Console.WriteLine("Precursor score: " + score);
                            // if (score < -0.5) continue; 
                            var portionExplainedFrags = 0.0;
                            for (var cutNumber = 1; cutNumber < pep.Length; cutNumber++)
                            {
                                //Console.WriteLine("Cut " + cutNumber);
                                var cutScore = imsScorer.GetCutScore(pep[cutNumber - 1], pep[cutNumber], sequence.GetComposition(0, cutNumber), precursorFeature);
                                //Console.WriteLine("{0} {1} {2} {3}", pep[cutNumber-1], pep[cutNumber], sequence.GetComposition(0, cutNumber), cutScore);
                                //   var cutNodeScore = imsScorer.GetCutNodeScore(pep[cutNumber - 1], pep[cutNumber], sequence.GetComposition(0, cutNumber), precursorFeature);
                                //   var cutRatioScore = imsScorer.GetCutRatioScore(pep[cutNumber - 1], pep[cutNumber], sequence.GetComposition(0, cutNumber), precursorFeature);
                                //   var cutLCScore = imsScorer.GetCutLcScore(pep[cutNumber - 1], pep[cutNumber], sequence.GetComposition(0, cutNumber), precursorFeature);
                                //   var cutIMSScore = imsScorer.GetCutImsScore(pep[cutNumber - 1], pep[cutNumber], sequence.GetComposition(0, cutNumber), precursorFeature);
                                //Console.Write(cutNumber + "\t" + cutNodeScore + "\t" + cutRatioScore + "\t" + cutLCScore + "\t" + cutIMSScore + "\t" +  cutScore);
                                //  Console.Write(cutNumber + "\t" + cutScore + "\t");
                                // foreach(var ion in imsScorer.supportingIonTypes) Console.Write("\t" + ion.Name+", ");
                                // Console.WriteLine();
                                score += cutScore;
                                portionExplainedFrags += imsScorer.supportingIonTypes.Count == 0 ? 0 : 1;
                            }
                            if (!(maxScore < score)) continue;
                            maxScore = score;
                            maxFeature = precursorFeature;
                            maxPortionOfExplainedFrag = portionExplainedFrags/sequence.Count;
                            //Console.WriteLine(i + " Score = " + score + "\n");
                        }
                    }

                    if (maxFeature != null)
                    {
                        if (!highestScorePerFeature.ContainsKey(maxFeature))
                            highestScorePerFeature[maxFeature] = new Tuple<double, bool>(maxScore, q == 0);
                        else
                        {
                            var prevScore = highestScorePerFeature[maxFeature].Item1;
                            if(maxScore > prevScore)
                                highestScorePerFeature[maxFeature] = new Tuple<double, bool>(maxScore, q == 0);
                        }  
                    }


                    if (!(maxScore > 0) || maxFeature == null || (highestScorePerFeature[maxFeature].Item1 > maxScore))
                        continue;

                    Console.WriteLine((q == 0 ? "T" : "D") + " " + num++ + " Max Score of the peptide " + pep + " is " + maxScore + " Portion of expalined fragmentation is " + maxPortionOfExplainedFrag);
                    Console.WriteLine("Corresponding max feature is " + maxFeature);
                    if (q == 0)
                    {
                        targetMatches.Add(new Tuple<double, string, Feature>(maxScore, pep, maxFeature));
                    }
                    else
                    {
                        decoyMatches.Add(new Tuple<double, string, Feature>(maxScore, pep, maxFeature));
                    }
                    // break;
                   // var scoreIndex = (int)Math.Min(targetDist.Length - 1, Math.Max(0, maxScore + 50));
                   // dist[scoreIndex] = dist[scoreIndex] + 1;

                }
            }


            var threshold = double.NegativeInfinity;

            foreach (var score in highestScorePerFeature.Values)
            {
                if (score.Item2) continue;
                threshold = Math.Max(threshold, score.Item1);
            }

            var numTarget = 0;
            foreach (var score in highestScorePerFeature.Values)
            {
                if (!score.Item2) continue;
                if (score.Item1 > threshold)
                {
                    numTarget++;
                }
            }

            // Print out target matches
            Console.WriteLine("Target matches");
            foreach (var match in targetMatches)
            {
                if (match.Item1 > threshold)
                {
                    Console.WriteLine("{0}\t{1}\t{2}", match.Item1, match.Item2, match.Item3);
                }
            }
            // Print out target matches
            Console.WriteLine("Decoy matches");
            foreach (var match in decoyMatches)
            {
                if (match.Item1 > threshold)
                {
                    Console.WriteLine("{0}\t{1}\t{2}", match.Item1, match.Item2, match.Item3);
                }
            }


            Console.WriteLine("Threshold is "+threshold + "\nNumber of target is "+numTarget);

            //var twriter = new StreamWriter(@"..\..\..\TestFiles\target.m");
            //var dwriter = new StreamWriter(@"..\..\..\TestFiles\decoy.m");
            
            //twriter.Write("t=[");
            //for (var j = 0; j < targetDist.Length;j++ )
            //{
            //   twriter.WriteLine(j+"\t"+targetDist[j]);
            //}
            //twriter.WriteLine("];");
            //twriter.Close();

            //dwriter.Write("d=[");
            //for (var j = 0; j < decoyDist.Length; j++)
            //{
            //    dwriter.WriteLine(j + "\t" + decoyDist[j]);
            //}
            //dwriter.WriteLine("];");
            //dwriter.Close();
        }
    }
}
