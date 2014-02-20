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

        private double getBestScore(IEnumerable<double> scores)
        {
            double bestScore = 0.0;
            double score = 0.0;
            if (scores != null)
            {
                foreach (var sc in scores)
                {
                    if (sc > score)
                    {
                        score = sc;
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
            var lcms = LcMsRun.GetLcMsRun(rawFileName, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);

            var cleanScans = CleanScans(lcms, scans, peptides, comp, thresh);

            var aset = new AminoAcidSet();
            var ionTypeFactory =
            new IonTypeFactory(new[] {BaseIonType.B, BaseIonType.Y, BaseIonType.C, BaseIonType.Z},
                                new[] {NeutralLoss.NoLoss}, TotalCharges);

            StreamWriter preHCDFile = File.AppendText(preOutFile);
            StreamWriter suffHCDFile = File.AppendText(suffOutFile);

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
                    var bFitScores = new List<double>();
                    var yFitScores = new List<double>();
                    for (int charge = 1; charge <= TotalCharges; charge++)
                    {
                        if (spec.ActivationMethod == ActivationMethod.HCD)
                        {
                            // HCD b
                            var bIon = ionTypeFactory.GetIonType("b" + ChargeToString(charge));
                            var b = bIon.GetIon(suffix);
                            b.Composition.ComputeApproximateIsotopomerEnvelop();
                            var bFitScore = spectrum.GetFitScore(b, new Tolerance(15, ToleranceUnit.Ppm), 0.1);
                            bFitScores.Add(bFitScore);
                            // HCD y
                            var yIon = ionTypeFactory.GetIonType("y" + ChargeToString(charge));
                            var y = yIon.GetIon(suffix);
                            y.Composition.ComputeApproximateIsotopomerEnvelop();
                            var yFitScore = spectrum.GetFitScore(y, new Tolerance(15, ToleranceUnit.Ppm), 0.1);
                            yFitScores.Add(yFitScore);
                        }

                    }
                    var bBestScore = getBestScore(bFitScores);
                    preHCDFile.WriteLine(bBestScore);
                }
            }
            preHCDFile.Close();
            suffHCDFile.Close();
        }

        [Test]
        public void FitScore()
        {
            const string fileList = @"\\protoapps\UserData\Wilkins\fileList.txt";
            const string preRes = @"\\protoapps\UserData\Sangtae\TopDownQCShew\msalign\";
            const string preRaw = @"\\protoapps\UserData\Sangtae\TopDownQCShew\raw\";
            const string outFile1 = @"\\protoapps\UserData\Wilkins\preHCD.txt";
            const string outFile2 = @"\\protoapps\UserData\Wilkins\suffHCD.txt";
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
