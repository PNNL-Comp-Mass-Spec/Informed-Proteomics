using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using DeconTools.Backend;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.FeatureFinding.MassFeature;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;

namespace InformedProteomics.Test.TopDownAnalysis
{
    class TestLikelihoodScorer
    {
        public void TestLoadTariningParam()
        {
            //const string paramPath = @"D:\MassSpecFiles\training\IdScoring\likelihoodTable";
            //var model = new LikelihoodScoringModel(paramPath);
        }

        [Test]
        public void TestGenerateFrequencyData()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string idFileFolder = @"D:\MassSpecFiles\training\IdScoring\MSPF_trainset";
            const string outFileFolder = @"D:\MassSpecFiles\training\IdScoring";

            if (!Directory.Exists(idFileFolder))
            {
                Assert.Ignore(@"Skipping test {0} since folder not found: {1}", methodName, idFileFolder);
            }

            Modification.RegisterAndGetModification(Modification.Cysteinyl.Name, Modification.Cysteinyl.Composition);
            Modification.RegisterAndGetModification(Modification.Phosphorylation.Name, Modification.Phosphorylation.Composition);
            Modification.RegisterAndGetModification(Modification.Methylation.Name, Modification.Methylation.Composition);
            Modification.RegisterAndGetModification(Modification.DiMethylation.Name, Modification.DiMethylation.Composition);
            Modification.RegisterAndGetModification(Modification.TriMethylation.Name, Modification.TriMethylation.Composition);
            Modification.RegisterAndGetModification("Trioxidation", new Composition(0, 0, 0, 3, 0));
            var aaSet = new AminoAcidSet(@"D:\MassSpecFiles\training\Mods.txt");

            var n = 0;

            for (var d = 0; d < TrainSetFileLists.Length; d++)
            {
                var dataset = TrainSetFileLists[d];
                var dataname = Path.GetFileNameWithoutExtension(dataset);
                var idFile = string.Format(@"{0}\{1}_IcTda.tsv", idFileFolder, dataname);
                var decoyFile = string.Format(@"{0}\{1}_IcDecoy.tsv", idFileFolder, dataname);
                var targetFile = string.Format(@"{0}\{1}_IcTarget.tsv", idFileFolder, dataname);

                if (!File.Exists(idFile)) continue;

                var prsmReader = new ProteinSpectrumMatchReader(0.01);
                var prsmList = prsmReader.LoadIdentificationResult(idFile);

                var minScore = prsmList.Last().Score;
                var decoyMatches = prsmReader.ReadMsPathFinderResult(decoyFile, int.MaxValue, 1, Math.Max(minScore - 5, 10));
                var run = PbfLcMsRun.GetLcMsRun(dataset);

                var spectrumMatchSet = LcMsFeatureTrain.CollectTrainSet(dataset, idFile);
                Console.WriteLine(spectrumMatchSet.Count);
                var writer = new StreamWriter(string.Format(@"{0}\{1}_target.tsv", outFileFolder, dataname));

                foreach (var matches in spectrumMatchSet)
                {
                    foreach (var match in matches)
                    {
                        var spec = run.GetSpectrum(match.ScanNum) as ProductSpectrum;
                        GetMatchStatistics(spec, match.GetSequence(), match.Charge, writer);
                    }
                }
                writer.Close();

                writer = new StreamWriter(string.Format(@"{0}\{1}_decoy.tsv", outFileFolder, dataname));
                foreach (var match in decoyMatches)
                {
                    var sequence = match.GetSequence();
                    var spec = run.GetSpectrum(match.ScanNum) as ProductSpectrum;
                    GetMatchStatistics(spec, sequence, match.Charge, writer);
                }
                writer.Close();
                n++;
            }
        }

        internal class FragmentStat
        {
            internal double Intensity;
            internal int Count;
            internal double MassError;
            internal double Corr;
            internal double Dist;
        }

        private void GetMatchStatistics(ProductSpectrum ms2Spec, Sequence sequence, int parentIonCharge, StreamWriter writer)
        {
            if (ms2Spec == null) return;
            if (sequence == null) return;

            var BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            var BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
            var tolerance = new Tolerance(12);
            var MinProductCharge = 1;
            var MaxProductCharge = Math.Min(parentIonCharge+2, 20);

            var baseIonTypes = ms2Spec.ActivationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;
            var refIntensity = CompositeScorer.GetRefIntensity(ms2Spec.Peaks);

            var activationMethodFlag = ms2Spec.ActivationMethod == ActivationMethod.ETD ? 2 : 1;
            var cleavages = sequence.GetInternalCleavages();
            var nComplementaryFrags = 0;

            var prefixStat = new FragmentStat();
            var suffixStat = new FragmentStat();

            var minMz = ms2Spec.Peaks.First().Mz;
            var maxMz = ms2Spec.Peaks.Last().Mz;

            var cleavageIndex = 0;

            var preFixIonCheck = new bool[sequence.Count + 1];
            var sufFixIonCheck = new bool[sequence.Count + 1];

            foreach (var c in cleavages)
            {
                var prefixHit = false;
                var suffixHit = false;

                foreach (var baseIonType in baseIonTypes)
                {
                    var stat = baseIonType.IsPrefix ? prefixStat : suffixStat;

                    var fragmentComposition = baseIonType.IsPrefix
                                  ? c.PrefixComposition + baseIonType.OffsetComposition
                                  : c.SuffixComposition + baseIonType.OffsetComposition;

                    if (fragmentComposition.Mass < ms2Spec.Peaks[0].Mz) continue;

                    var curFragMass = fragmentComposition.Mass;
                    /*var curObsIonCharge = 0;
                    var curObsIonDist = 1.0d;
                    var curObsIonCorr = 0d;
                    var curObsIonIntensity = 0d;
                    var curObsIonMassError = 0d;*/

                    var mostAbundantIsotopeIndex = fragmentComposition.GetMostAbundantIsotopeZeroBasedIndex();
                    var fragmentIonMostAbuMass = fragmentComposition.Mass + Constants.C13MinusC12 * mostAbundantIsotopeIndex;

                    var maxCharge = (int)Math.Floor(fragmentIonMostAbuMass / (minMz - Constants.Proton));
                    var minCharge = (int)Math.Ceiling(fragmentIonMostAbuMass / (maxMz - Constants.Proton));
                    if (maxCharge < 1 || maxCharge > MaxProductCharge) maxCharge = MaxProductCharge;
                    if (minCharge < 1 || minCharge < MinProductCharge) minCharge = MinProductCharge;

                    //var ionMatch = false;
                    for (var charge = minCharge; charge <= maxCharge; charge++)
                    {
                        var ion = new Ion(fragmentComposition, charge);

                        var isotopePeaks = ms2Spec.GetAllIsotopePeaks(ion, tolerance, 0.1);

                        if (isotopePeaks == null) continue;

                        var distCorr = CompositeScorer.GetDistCorr(ion, isotopePeaks);
                        if (distCorr.Item2 < 0.7 && distCorr.Item1 > 0.03) continue;
                        var mostAbuPeak = isotopePeaks[mostAbundantIsotopeIndex];
                        var intScore = mostAbuPeak.Intensity / refIntensity;
                        /*
                        if (ionMatch == false || curObsIonIntensity < intScore)
                        {
                            curObsIonCharge = charge;
                            curObsIonCorr = distCorr.Item2;
                            curObsIonDist = distCorr.Item1;
                            curObsIonIntensity = intScore;

                            var mostAbuPeakMz = Ion.GetIsotopeMz(curFragMass, charge, mostAbundantIsotopeIndex);
                            curObsIonMassError = (Math.Abs(mostAbuPeak.Mz - mostAbuPeakMz) / mostAbuPeakMz) * 1e6;

                            //var curObsIonMass = Ion.GetMonoIsotopicMass(mostAbuPeak.Mz, charge, mostAbundantIsotopeIndex);
                            //curObsIonMassError = (Math.Abs(curFragMass - curObsIonMass) / curFragMass) * 1e6;
                        }
                        ionMatch = true;
                        */
                        var mostAbuPeakMz = Ion.GetIsotopeMz(curFragMass, charge, mostAbundantIsotopeIndex);
                        var curObsIonMassError = (Math.Abs(mostAbuPeak.Mz - mostAbuPeakMz) / mostAbuPeakMz) * 1e6;

                        stat.Count++;
                        stat.Intensity += Math.Min(intScore, 1.0);
                        stat.Corr += distCorr.Item2;
                        stat.Dist += distCorr.Item1;
                        stat.MassError += curObsIonMassError;

                        if (baseIonType.IsPrefix) prefixHit = true;
                        else suffixHit = true;
                    }

                    //if (!ionMatch) continue;
                }

                if (prefixHit) preFixIonCheck[cleavageIndex] = true;
                if (suffixHit) sufFixIonCheck[cleavageIndex] = true;

                if (prefixHit && suffixHit)
                {
                    nComplementaryFrags++;
                }
                cleavageIndex++;
            }

            var preContCount = 0;
            var sufContCount = 0;
            for (var i = 0; i < preFixIonCheck.Length - 1; i++)
            {
                if (preFixIonCheck[i] && preFixIonCheck[i + 1]) preContCount++;
                if (sufFixIonCheck[i] && sufFixIonCheck[i + 1]) sufContCount++;
            }

            writer.Write(activationMethodFlag);
            writer.Write("\t");
            writer.Write(sequence.Composition.Mass);
            writer.Write("\t");
            writer.Write(sequence.Count);
            writer.Write("\t");
            writer.Write(nComplementaryFrags);
            writer.Write("\t");

            writer.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", prefixStat.Count, preContCount, prefixStat.Intensity, prefixStat.Corr, prefixStat.Dist, prefixStat.MassError);
            writer.Write("\t");
            writer.Write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", suffixStat.Count, sufContCount, suffixStat.Intensity, suffixStat.Corr, suffixStat.Dist, suffixStat.MassError);
            writer.Write("\n");
        }

        private void GetNodeStatistics(bool isDecoy, ProductSpectrum ms2Spec, Sequence sequence, StreamWriter writer) //, StreamWriter mzErrorWriter)
        {
            if (ms2Spec == null) return;
            if (sequence == null) return;
            //var refIntensity = ms2Spec.Peaks.Max(p => p.Intensity) * 0.01;
            //refIntensity = Math.Min(ms2Spec.Peaks.Select(p => p.Intensity).Median(), refIntensity);

            var BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            var BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
            var tolerance = new Tolerance(15);
            var minCharge = 1;
            var maxCharge = 20;
            var baseIonTypes = ms2Spec.ActivationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;

            var refIntensity = ms2Spec.Peaks.Max(p => p.Intensity);

            var activationMethodFlag = ms2Spec.ActivationMethod == ActivationMethod.ETD ? 1 : 2;
            var cleavages = sequence.GetInternalCleavages();

            var prevPrefixFragMass = 0d;
            var prevPrefixObsIonMass = 0d;
            var prevPrefixObsIonCharge = 0;
            var prevPrefixObsIonIntensity = 0d;

            var prevSuffixFragMass = 0d;
            var prevSuffixObsIonMass = 0d;
            var prevSuffixObsIonCharge = 0;
            var prevSuffixObsIonIntensity = 0d;

            var nComplementaryFrags = 0;

            var cleavageIndex = 0;

            foreach (var c in cleavages)
            {
                var bothObs = true;
                foreach (var baseIonType in baseIonTypes)
                {
                    var peakType = baseIonType.IsPrefix ? 1 : 2; // unexplained
                    var fragmentComposition = baseIonType.IsPrefix
                                  ? c.PrefixComposition + baseIonType.OffsetComposition
                                  : c.SuffixComposition + baseIonType.OffsetComposition;

                    var curFragMass = fragmentComposition.Mass;
                    var curObsIonMass = 0d;
                    var curObsIonCharge = 0;
                    var curObsIonDist = 1.0d;
                    var curObsIonCorr = 0d;
                    var curObsIonIntensity = 0d;

                    var ionMatch = false;
                    for (var charge = minCharge; charge <= maxCharge; charge++)
                    {
                        var ion = new Ion(fragmentComposition, charge);

                        var isotopePeaks = ms2Spec.GetAllIsotopePeaks(ion, tolerance, 0.1);

                        if (isotopePeaks == null) continue;

                        var distCorr = AbstractFragmentScorer.GetDistCorr(ion, isotopePeaks);

                        if (distCorr.Item2 < 0.7 && distCorr.Item1 > 0.07) continue;
                        var mostAbundantIsotopeIndex = ion.Composition.GetMostAbundantIsotopeZeroBasedIndex();
                        var mostAbuPeak = isotopePeaks[mostAbundantIsotopeIndex];

                        var summedIntensity = isotopePeaks.Where(p => p != null).Sum(p => p.Intensity);
                        var intScore = summedIntensity/refIntensity;
                        //var intScore = mostAbuPeak.Intensity / medIntensity;
                        //var intScore = summedIntensity / refIntensity;

                        if (ionMatch == false || curObsIonIntensity < intScore)
                        {
                            curObsIonMass = Ion.GetMonoIsotopicMass(mostAbuPeak.Mz, charge, mostAbundantIsotopeIndex);
                            curObsIonCharge = charge;
                            curObsIonCorr = distCorr.Item2;
                            curObsIonDist = distCorr.Item1;
                            curObsIonIntensity = intScore;
                        }
                        ionMatch = true;
                    }

                    if (!ionMatch)
                    {
                        bothObs = false;
                        continue;
                    }

                    writer.Write(activationMethodFlag);
                    writer.Write("\t");

                    writer.Write(peakType);
                    writer.Write("\t");

                    writer.Write("{0:0.000}", curFragMass);
                    writer.Write("\t");

                    writer.Write("{0}", curObsIonCharge);
                    writer.Write("\t");

                    writer.Write("{0:0.000}", curObsIonDist);
                    writer.Write("\t");

                    writer.Write("{0:0.000}", curObsIonCorr);
                    writer.Write("\t");

                    writer.Write("{0:0.000}", curObsIonIntensity);
                    writer.Write("\t");

                    writer.Write("{0:0.000}", (Math.Abs(curFragMass - curObsIonMass)/curFragMass)*1e6);
                    writer.Write("\n");

                    // mz error output
                    /*
                    if (baseIonType.IsPrefix && prevPrefixFragMass > 0 & prevPrefixObsIonMass > 0)
                    {
                        var aaMass = Math.Abs(prevPrefixFragMass - curFragMass);
                        var massError = Math.Abs(Math.Abs(prevPrefixObsIonMass - curObsIonMass) - aaMass);
                        var massErrorPpm = (massError / curObsIonMass) * 1e6;
                         mzErrorWriter.WriteLine("{0}\t{1:0.000}\t{2}", activationMethodFlag, massErrorPpm, Math.Abs(prevPrefixObsIonCharge - curObsIonCharge));
                    }
                    else if (prevSuffixFragMass > 0 & prevSuffixObsIonMass > 0)
                    {
                        var aaMass = Math.Abs(prevSuffixFragMass - curFragMass);
                        var massError = Math.Abs(Math.Abs(prevSuffixObsIonMass - curObsIonMass) - aaMass);
                        var massErrorPpm = (massError / curObsIonMass) * 1e6;
                        mzErrorWriter.WriteLine("{0}\t{1:0.000}\t{2}", activationMethodFlag, massErrorPpm, Math.Abs(prevSuffixObsIonCharge - curObsIonCharge));
                    }
                    */
                    if (baseIonType.IsPrefix)
                    {
                        prevPrefixFragMass = curFragMass;
                        prevPrefixObsIonMass = curObsIonMass;
                        prevPrefixObsIonCharge = curObsIonCharge;
                        prevPrefixObsIonIntensity = curObsIonIntensity;
                        //Array.Copy(curObsIonMass, prevPrefixObsIonMass, curObsIonMass.Length);
                    }
                    else
                    {
                        prevSuffixFragMass = curFragMass;
                        prevSuffixObsIonMass = curObsIonMass;
                        prevSuffixObsIonCharge = curObsIonCharge;
                        prevSuffixObsIonIntensity = curObsIonIntensity;
                        //Array.Copy(curObsIonMass, prevSuffixObsIonMass, curObsIonMass.Length);
                    }
                }

                if (bothObs)
                {
                    //pairWriter.Write("{0}\t{1}\t", prevPrefixObsIonIntensity, prevSuffixObsIonIntensity);
                    nComplementaryFrags++;
                }
                cleavageIndex++;
            }

            Console.WriteLine("{0}\t{1}", nComplementaryFrags, sequence.Count);

            //if (!isDecoy) Console.WriteLine("{0}", totalExplainedAbundanceRatio);
        }

        public Peak[] GetAllIsotopePeaks(Spectrum spec, Ion ion, Tolerance tolerance, double relativeIntensityThreshold, out int[] peakIndexList)
        {
            var mostAbundantIsotopeIndex = ion.Composition.GetMostAbundantIsotopeZeroBasedIndex();
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelopeRelativeIntensities();

            peakIndexList = new int[isotopomerEnvelope.Length];

            var mostAbundantIsotopeMz = ion.GetIsotopeMz(mostAbundantIsotopeIndex);
            var mostAbundantIsotopeMatchedPeakIndex = spec.FindPeakIndex(mostAbundantIsotopeMz, tolerance);
            if (mostAbundantIsotopeMatchedPeakIndex < 0) return null;

            var observedPeaks = new Peak[isotopomerEnvelope.Length];
            observedPeaks[mostAbundantIsotopeIndex] = spec.Peaks[mostAbundantIsotopeMatchedPeakIndex];
            peakIndexList[mostAbundantIsotopeIndex] = mostAbundantIsotopeMatchedPeakIndex;

            // go down
            var peakIndex = mostAbundantIsotopeMatchedPeakIndex - 1;
            for (var isotopeIndex = mostAbundantIsotopeIndex - 1; isotopeIndex >= 0; isotopeIndex--)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex; i >= 0; i--)
                {
                    var peakMz = spec.Peaks[i].Mz;
                    if (peakMz < minMz)
                    {
                        peakIndex = i;
                        break;
                    }
                    if (peakMz <= maxMz)    // find match, move to prev isotope
                    {
                        var peak = spec.Peaks[i];
                        if (observedPeaks[isotopeIndex] == null ||
                            peak.Intensity > observedPeaks[isotopeIndex].Intensity)
                        {
                            observedPeaks[isotopeIndex] = peak;
                            peakIndexList[isotopeIndex] = i;
                        }
                    }
                }
            }

            // go up
            peakIndex = mostAbundantIsotopeMatchedPeakIndex + 1;
            for (var isotopeIndex = mostAbundantIsotopeIndex + 1; isotopeIndex < isotopomerEnvelope.Length; isotopeIndex++)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;
                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex; i < spec.Peaks.Length; i++)
                {
                    var peakMz = spec.Peaks[i].Mz;
                    if (peakMz > maxMz)
                    {
                        peakIndex = i;
                        break;
                    }
                    if (peakMz >= minMz)    // find match, move to prev isotope
                    {
                        var peak = spec.Peaks[i];
                        if (observedPeaks[isotopeIndex] == null ||
                            peak.Intensity > observedPeaks[isotopeIndex].Intensity)
                        {
                            observedPeaks[isotopeIndex] = peak;
                            peakIndexList[isotopeIndex] = i;
                        }
                    }
                }
            }

            return observedPeaks;
        }

        public static string[] TrainSetFileLists = new string[]
            {
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_1_11Feb15_Bane_C2Column5.pbf",// top-down datasets
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_2_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_BE100_PO4_3_11Feb15_Bane_C2Column5.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep3_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep4_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep5_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep7_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep8_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_01_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_02_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_03_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_FA_01_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_FA_02_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_1\CPTAC_Intact_SDS_T_FA_03_15Jan15_Bane_C2-14-08-02RZ.pbf",
                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_3_2Feb15_Bane_C2Column4.pbf",
                @"D:\MassSpecFiles\training\raw\QC_Shew_Intact_26Sep14_Bane_C2Column3.pbf",
                @"D:\MassSpecFiles\test\NewQC_LongSep_29Sep14_141001104925.pbf",
                @"D:\MassSpecFiles\training\raw\YS_Shew_testHCD_CID.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2014_2\yufeng_column_test2.pbf",

                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\MZ20150529DS_histone41.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\MZ20150529DS_histone43.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\MZ20150529DS_histone44.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\MZ20150529DS_histone55.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\MZ20150529DS_histone46.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\MZ20150529DS_histone47.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\MZ20150529DS_histone49.pbf",
                @"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2015_2\MZ20150529DS_histone50.pbf",

                //@"\\proto-11\MSXML_Cache\PBF_Gen_1_193\2014_3\Lewy_intact_01.pbf"
            };
    }
}
