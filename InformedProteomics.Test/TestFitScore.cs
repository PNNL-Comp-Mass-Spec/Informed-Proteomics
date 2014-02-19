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
                    if (peak.Intensity > intensity)
                    {
                        intensity = peak.Intensity;
                        mostIntense = peak;
                    }
                }
            }
            return mostIntense;
        }

        private Tuple<double, double> getBestScore(IEnumerable<Tuple<double, double>> scores)
        {
            Tuple<double, double> bestScore = null;
            double score = 0.0;
            if (scores != null)
            {
                foreach (var sc in scores)
                {
                    if (sc.Item2 > score)
                    {
                        score = sc.Item2;
                        bestScore = sc;
                    }
                }
            }
            return bestScore;

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

        private void WriteFitScore(string txtFileName, string rawFileName, string outFile, bool useDecoy)
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
            var lcms = LcMsRun.GetLcMsRun(rawFileName, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

            var cleanScans = CleanScans(lcms, scans, peptides, comp, thresh);

            var aset = new AminoAcidSet();
            var ionTypeFactory =
            new IonTypeFactory(new[] {BaseIonType.B, BaseIonType.Y, BaseIonType.C, BaseIonType.Z},
                                new[] {NeutralLoss.NoLoss}, TotalCharges);

            StreamWriter suffHCDFile = File.AppendText(outFile);

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
                if (spec == null) continue;
                for (int i = 0; i < protein.Length - 1; i++)
                {
                    var suffix = sequence.GetComposition(protein.Length - i, protein.Length);
                    var fitScores = new List<Tuple<double, double>>();
                    for (int charge = 1; charge <= TotalCharges; charge++)
                    {
                        var yIon = ionTypeFactory.GetIonType("y" + ChargeToString(charge));
                        var y = yIon.GetIon(suffix);
                        y.Composition.ComputeApproximateIsotopomerEnvelop();

                        if (spec.ActivationMethod == ActivationMethod.HCD)
                        {
                            // HCD Y
                            var yPeaks = spectrum.GetIonPeaks(y, new Tolerance(15, ToleranceUnit.Ppm), 0.8);
                            Peak tallesty = null;
                            var yFitScore = spectrum.GetFitScore(y, new Tolerance(15, ToleranceUnit.Ppm), 0.1);
                            if (yPeaks != null)
                                tallesty = GetMostIntensePeak(yPeaks);
                            if (tallesty != null)
                            {
                                fitScores.Add(new Tuple<double, double>(tallesty.Intensity, yFitScore));
                            }
                        }

                    }
                    var bestScore = getBestScore(fitScores);
                    if (bestScore != null)
                        suffHCDFile.WriteLine(bestScore.Item1 + "\t" + bestScore.Item2);
                }
            }
            suffHCDFile.Close();
        }

        [Test]
        public void FitScore()
        {
            const string fileList = @"C:\Users\wilk011\Documents\DataFiles\fileList.txt";
            const string pre = @"\\protoapps\UserData\Sangtae\ForChris\";
            const string outFile = @"C:\Users\wilk011\Documents\DataFiles\IntensityHistograms\FitScore\suffHCD2.txt";
            const bool decoy = false;

            var fileNameParser = new TsvFileParser(fileList);

            var txtFiles = fileNameParser.GetData("text");
            var rawFiles = fileNameParser.GetData("raw");

            using (var txtFileIt = txtFiles.GetEnumerator())
            using (var rawFileIt = rawFiles.GetEnumerator())
            {
                while (txtFileIt.MoveNext() && rawFileIt.MoveNext())
                {
                    WriteFitScore(pre + txtFileIt.Current, pre + rawFileIt.Current, outFile, decoy);
                }
            }
        }

    }
}
