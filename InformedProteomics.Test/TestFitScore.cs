using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
//using System.Text;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    class TestFitScore
    {
        private const double FdrThreshold = 0.01;
        private const double EValueThreshold = 0.1;
        private const int TotalCharges = 10;
        private const int NumMutations = 3;

        private static string Trim(string prot)
        {
            int start = prot.IndexOf('.') + 1;
            int length = prot.LastIndexOf('.') - start;
            return prot.Substring(start, length);
        }
        private static List<Tuple<string, Spectrum>> CleanScans(LcMsRun lcms, IEnumerable<string> scans,
                                IEnumerable<string> peptides, IEnumerable<string> minValues,
                                double thresh)
        {
            var clean = new List<Tuple<string, Spectrum>>();
            using (var scani = scans.GetEnumerator())
            using (var peptidei = peptides.GetEnumerator())
            using (var min = minValues.GetEnumerator())
            {
                while (scani.MoveNext() && peptidei.MoveNext() && min.MoveNext())
                {
                    var spec = lcms.GetSpectrum(Convert.ToInt32(scani.Current));
                    double minValue = Convert.ToDouble(min.Current);
                    if (peptidei.Current.Contains('[') || peptidei.Current.Contains('U') || minValue > thresh) continue;
                    clean.Add(new Tuple<string, Spectrum>(Trim(peptidei.Current), spec));
                }
            }
            return clean;
        }

        public Peak GetMostIntensePeak(IEnumerable<Peak> peaks)
        {
            Peak mostIntense = null;
            double intensity = 0.0;
            if (peaks != null)
            {
                foreach (var peak in peaks)
                {
                    if (peak != null && peak.Intensity > intensity)
                    {
                        intensity = peak.Intensity;
                        mostIntense = peak;
                    }
                }
            }
            return mostIntense;
        }

        private Tuple<double, double> getBestScore(IEnumerable<Tuple<double,double>> scores)
        {
            Tuple<double, double> score = null;
            double bestScore = 0.0;

            foreach (var sc in scores)
            {
                if (sc.Item2 > bestScore)
                {
                    score = sc;
                    bestScore = sc.Item2;
                }
            }
            return score;

        }

        private string ChargeToString(int charge)
        {
            string chargeStr = "";
            if (charge > 1)
            {
                chargeStr = charge.ToString(CultureInfo.InvariantCulture);
            }
            return chargeStr;
        }

        private void WriteFitScore(string txtFileName, string rawFileName, string preOutFile, string suffOutFile, bool useDecoy)
        {
            Console.WriteLine(rawFileName);
            var tsvParser = new TsvFileParser(txtFileName);

            var scans = tsvParser.GetData("Scan(s)");
            var peptides = tsvParser.GetData("Peptide");
            var fdrs = tsvParser.GetData("FDR");
            var evalues = tsvParser.GetData("E-value");
            IEnumerable<string> comp = fdrs;
            double thresh = FdrThreshold;
            if (fdrs == null)
            {
                comp = evalues;
                thresh = EValueThreshold;
            }
            var lcms = LcMsRun.GetLcMsRun(rawFileName, MassSpecDataType.XCaliburRun, 0, 0);

            var cleanScans = CleanScans(lcms, scans, peptides, comp, thresh);

            var aset = new AminoAcidSet();
            var ionTypeFactory =
            new IonTypeFactory(new[] {BaseIonType.B, BaseIonType.Y, BaseIonType.C, BaseIonType.Z},
                                new[] {NeutralLoss.NoLoss}, TotalCharges);

            StreamWriter preFile = File.AppendText(preOutFile);
            StreamWriter suffFile = File.AppendText(suffOutFile);

            foreach (var node in cleanScans)
            {
                var protein = node.Item1;
                var spectrum = node.Item2;
                if (useDecoy)
                {
                    var shuffled = SimpleStringProcessing.Shuffle(protein);
                    protein = SimpleStringProcessing.Mutate(shuffled, NumMutations);
                }
                var sequence = new Sequence(protein, aset);
                var spec = spectrum as ProductSpectrum;
                if (spec == null || spec.ActivationMethod != ActivationMethod.CID) continue;
                for (int i = 1; i < protein.Length; i++)
                {
                    var suffix = sequence.GetComposition(protein.Length - i, protein.Length);
                    var bFitScores = new List<Tuple<double, double>>();
                    Peak[] bPeaks;
                    var yFitScores = new List<Tuple<double,double>>();
                    Peak[] yPeaks;
                    for (int charge = 1; charge <= TotalCharges; charge++)
                    {
                            // b
                            var bIon = ionTypeFactory.GetIonType("b" + ChargeToString(charge));
                            var b = bIon.GetIon(suffix);
                            b.Composition.ComputeApproximateIsotopomerEnvelop();
                            var bFitScore = spectrum.GetConsineScore(b, new Tolerance(15, ToleranceUnit.Ppm), 0.1);
                            bPeaks = spectrum.GetAllIsotopePeaks(b, new Tolerance(15, ToleranceUnit.Ppm), 0.1);
                            var highestBPeak = GetMostIntensePeak(bPeaks);
                            double highestBIntensity = -1.0;
                            if (highestBPeak != null)
                                highestBIntensity = highestBPeak.Intensity;
                            bFitScores.Add(new Tuple<double, double>(highestBIntensity, bFitScore));
                            // y
                            var yIon = ionTypeFactory.GetIonType("y" + ChargeToString(charge));
                            var y = yIon.GetIon(suffix);
                            y.Composition.ComputeApproximateIsotopomerEnvelop();
                            var yFitScore = spectrum.GetConsineScore(y, new Tolerance(15, ToleranceUnit.Ppm), 0.1);
                            yPeaks = spectrum.GetAllIsotopePeaks(y, new Tolerance(15, ToleranceUnit.Ppm), 0.1);
                            var highestYPeak = GetMostIntensePeak(yPeaks);
                            double highestYIntensity = -1.0;
                            if (highestYPeak != null)
                                highestYIntensity = highestYPeak.Intensity;
                            yFitScores.Add(new Tuple<double, double>(highestYIntensity, yFitScore));
                    }

                    var bBestScore = getBestScore(bFitScores);
                    preFile.WriteLine("{0}\t{1}", bBestScore.Item1, bBestScore.Item2);

                    var yBestScore = getBestScore(yFitScores);
                    suffFile.WriteLine("{0}\t{1}", yBestScore.Item1, yBestScore.Item2);
                }
            }
            preFile.Close();
            suffFile.Close();
        }

        [Test]
        public void FitScore()
        {
            const string fileList = @"C:\Users\wilk011\Documents\DataFiles\oldfileList.txt";
            const string preRes = @"C:\Users\wilk011\Documents\DataFiles\ForChris\";
            const string preRaw = @"C:\Users\wilk011\Documents\DataFiles\ForChris\";
            const string outFile1 = @"C:\Users\wilk011\Documents\DataFiles\Cosine\bCID.txt";
            const string outFile2 = @"C:\Users\wilk011\Documents\DataFiles\Cosine\yCID.txt";
            const bool decoy = false;

            var fileNameParser = new TsvFileParser(fileList);

            var txtFiles = fileNameParser.GetData("text");
            var rawFiles = fileNameParser.GetData("raw");

            using (var txtFileIt = txtFiles.GetEnumerator())
            using (var rawFileIt = rawFiles.GetEnumerator())
            {
                while (txtFileIt.MoveNext() && rawFileIt.MoveNext())
                {
                    WriteFitScore(preRes + txtFileIt.Current, preRaw + rawFileIt.Current, outFile1, outFile2, decoy);
                }
            }
        }

    }
}
