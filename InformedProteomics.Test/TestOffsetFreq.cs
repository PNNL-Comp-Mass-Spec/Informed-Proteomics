using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class OffsetTable
    {
        private const double FdrThreshold = 0.01;
        private const int TotalCharges = 20;
        private static string Trim(string prot)
        {
            int start = prot.IndexOf('.') + 1;
            int length = prot.LastIndexOf('.') - start;
            return prot.Substring(start, length);
        }
        private static List<int> CleanScans(ref IList<string> scans, ref IList<string> peptides)
        {
            var cleanScans = new List<int>();
            var cleanPeptides = new List<string>();
            using (var scani = scans.GetEnumerator())
            using (var peptidei = peptides.GetEnumerator())
            //            using (var fdri = fdrs.GetEnumerator())
            {
                while (scani.MoveNext() && peptidei.MoveNext())
                {
                    // double fdr = Convert.ToDouble(fdri.Current);
                    if  (peptidei.Current.Contains('[') || peptidei.Current.Contains('U')) continue;
                    cleanScans.Add(Convert.ToInt32(scani.Current));
                    cleanPeptides.Add(Trim(peptidei.Current));
                }
            }
            peptides = cleanPeptides;
            return cleanScans;
        }
        private static string ChargeToString(int charge)
        {
            string chargeStr = "";
            if (charge > 1)
            {
                chargeStr = charge.ToString();
            }
            return chargeStr;
        }

        private static void WriteOffSet(string txtFileName, string rawFileName)
        {
            Console.WriteLine(rawFileName);
            string outFile = rawFileName + ".out";
            var tsvParser = new TsvFileParser(txtFileName);

            var scans = tsvParser.GetData("Scan(s)");
            var peptides = tsvParser.GetData("Peptide");

            var cleanScans = CleanScans(ref scans, ref peptides);
            var lcms = LcMsRun.GetLcMsRun(rawFileName, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
//            var lcms = new LcMsRun(new XCaliburReader(rawFileName));
            var spectra = cleanScans.Select(lcms.GetSpectrum).ToList();

            var aset = new AminoAcidSet();
            var ionTypeFactory =
                new IonTypeFactory(new[] { BaseIonType.B, BaseIonType.Y, BaseIonType.C, BaseIonType.Z },
                                    new[] { NeutralLoss.NoLoss }, TotalCharges);


            using (var file = new StreamWriter(outFile))
            {
                file.WriteLine("charge\tb\ty\tc\tz");
                for (int charge = 1; charge <= TotalCharges; charge++)
                {
                    file.Write("{0}\t", charge);
                    int btotal = 0, bfound = 0;
                    int ytotal = 0, yfound = 0;
                    int ctotal = 0, cfound = 0;
                    int ztotal = 0, zfound = 0;
                    using (var spectrum = spectra.GetEnumerator())
                    using (var peptide = peptides.GetEnumerator())
                    {
                        while (spectrum.MoveNext() && peptide.MoveNext())
                        {
                            var sequence = new Sequence(peptide.Current, aset);
                            var bIon = ionTypeFactory.GetIonType("b" + ChargeToString(charge));
                            var yIon = ionTypeFactory.GetIonType("y" + ChargeToString(charge));
                            var cIon = ionTypeFactory.GetIonType("c" + ChargeToString(charge));
                            var zIon = ionTypeFactory.GetIonType("z" + ChargeToString(charge));
                            for (int i = 0; i < peptide.Current.Length - 1; i++)
                            {
                                btotal++;
                                ytotal++;
                                ctotal++;
                                ztotal++;
                                var prefix = sequence.GetComposition(0, i);
                                var suffix = sequence.GetComposition(peptide.Current.Length - i,
                                                                        peptide.Current.Length);
                                var b = bIon.GetIon(prefix);
                                var y = yIon.GetIon(suffix);
                                var c = cIon.GetIon(prefix);
                                var z = zIon.GetIon(suffix);
                                var bmz = b.GetMonoIsotopicMz();
                                var ymz = y.GetMonoIsotopicMz();
                                var cmz = c.GetMonoIsotopicMz();
                                var zmz = z.GetMonoIsotopicMz();
                                var bpeak = spectrum.Current.FindPeak(bmz, new Tolerance(15, ToleranceUnit.Ppm));
                                var ypeak = spectrum.Current.FindPeak(ymz, new Tolerance(15, ToleranceUnit.Ppm));
                                var cpeak = spectrum.Current.FindPeak(cmz, new Tolerance(15, ToleranceUnit.Ppm));
                                var zpeak = spectrum.Current.FindPeak(zmz, new Tolerance(15, ToleranceUnit.Ppm));
                                if (bpeak != null) bfound++;
                                if (ypeak != null) yfound++;
                                if (cpeak != null) cfound++;
                                if (zpeak != null) zfound++;
                            }
                        }
                    }
                    file.Write("{0}/{1}\t", bfound, btotal);
                    file.Write("{0}/{1}\t", yfound, ytotal);
                    file.Write("{0}/{1}\t", cfound, ctotal);
                    file.Write("{0}/{1}", zfound, ztotal);
                    file.WriteLine();
                }
            }
        }

        private static void InitList(int size, ref int[] foundList, ref int[] totalList)
        {
            for (int i = 0; i < size; i++)
            {
                foundList[i] = 0;
                totalList[i] = 0;
            }
        }

        private static void CalcData(IList<string> data, ref int[] foundList, ref int[] totalList)
        {
            int count = 0;
            foreach (var datum in data)
            {
                var segments = datum.Split('/');
                int found = Convert.ToInt32(segments[0].Trim());
                int total = Convert.ToInt32(segments[1].Trim());
                foundList[count] += found;
                totalList[count] += total;
                count++;
            }
        }

        [Test]
        public static void OffsetFreq(string[] args)
        {

            const string pre = @"\\protoapps\UserData\Sangtae\ForChris\";
            const string offsetFile = @"\\protoapps\UserData\Sangtae\ForChris\offsetCounts.txt";

            var fileNameParser = new TsvFileParser(@"C:\Users\wilk011\Documents\DataFiles\fileList.txt");

            var txtFiles = fileNameParser.GetData("text");
            var rawFiles = fileNameParser.GetData("raw");

            using (var txtFileIt = txtFiles.GetEnumerator())
            using (var rawFileIt = rawFiles.GetEnumerator())
            {
                while (txtFileIt.MoveNext() && rawFileIt.MoveNext())
                {
                    WriteOffSet(pre+txtFileIt.Current, pre+rawFileIt.Current);
                }
            }

            int[] bFound = new int[TotalCharges];
            int[] bTotal = new int[TotalCharges];
            InitList(TotalCharges, ref bFound, ref bTotal);
            int[] cFound = new int[TotalCharges];
            int[] cTotal = new int[TotalCharges];
            InitList(TotalCharges, ref cFound, ref cTotal);
            int[] yFound = new int[TotalCharges];
            int[] yTotal = new int[TotalCharges]; 
            InitList(TotalCharges, ref yFound, ref yTotal);
            int[] zFound = new int[TotalCharges];
            int[] zTotal = new int[TotalCharges];
            InitList(TotalCharges, ref zFound, ref zTotal);

            foreach (var rawFile in rawFiles)
            {
                string outFile = pre + rawFile + ".out";
                var outParser = new TsvFileParser(outFile);
                var bData = outParser.GetData("b");
                var cData = outParser.GetData("c");
                var yData = outParser.GetData("y");
                var zData = outParser.GetData("z");
                CalcData(bData, ref bFound, ref bTotal);
                CalcData(cData, ref cFound, ref cTotal);
                CalcData(yData, ref yFound, ref yTotal);
                CalcData(zData, ref zFound, ref zTotal);
            }

            using (var finalOutFile = new StreamWriter(offsetFile))
            {
                finalOutFile.WriteLine("charge\tb\tc\ty\tz");
                for (int i = 0; i < TotalCharges; i++)
                {
                    finalOutFile.Write((i+1)+"\t");
                    finalOutFile.Write("{0}\t", Math.Round((double)(bFound[i])/ (bTotal[i]), 5));
                    finalOutFile.Write("{0}\t", Math.Round((double)(cFound[i]) / (cTotal[i]), 5));
                    finalOutFile.Write("{0}\t", Math.Round((double)(yFound[i]) / (yTotal[i]), 5));
                    finalOutFile.Write("{0}", Math.Round((double)(zFound[i]) / (zTotal[i]), 5));
                    finalOutFile.WriteLine();
                }
            }
        }
    }
}
