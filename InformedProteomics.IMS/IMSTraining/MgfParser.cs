using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.IMS.IMSTraining
{
    public class MgfParser
    {
        private static readonly AminoAcidSet AminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);

        public static List<MSMSSpectrum> Parse(string fileName)
        {
            return Parse(fileName, int.MaxValue);
        }
        
        public static List<MSMSSpectrum> Parse(string fileName, int maxNumberOfSpectra)
        {
            var  spectra = new List<MSMSSpectrum>();
            var reader = new StreamReader(fileName);
            String s;
            var t = "";
            while ((s=reader.ReadLine())!=null)
            {
                if (s.StartsWith("BEGIN ION"))
                {
                    t = "";
                    continue;
                }
                t = t + s + "\n";
                if (s.StartsWith("END ION"))
                {
                    spectra.Add(ReadSpectrum(t));
                    if (spectra.Count >= maxNumberOfSpectra) break;
                }
            }
            reader.Close();
            return spectra;
        }

        private static MSMSSpectrum ReadSpectrum(string s)
        {
            var precursorMz = 0.0;
            var charge = 0;
            Sequence annotation = null;
            var peaks = new List<MSMSSpectrumPeak>();
            var token = s.Split('\n');
            foreach (var t in token)
            {
                if (t.Length == 0) continue;
                if (char.IsDigit(t[0]))
                {
                    var p = t.Split(new [] {"\t", " "}, StringSplitOptions.None);
                    peaks.Add(new MSMSSpectrumPeak(double.Parse(p[0]), double.Parse(p[1])));
                }
                else if (t.StartsWith("CHARGE"))
                {
                    var chargeStr = t.Substring(t.IndexOf('=') + 1).Trim();
                    if (chargeStr.StartsWith("+")) chargeStr = chargeStr.Substring(1);
                    if (chargeStr.EndsWith("+")) chargeStr = chargeStr.Substring(0, chargeStr.Length - 1);
                    charge = int.Parse(chargeStr);
                }
                else if (t.StartsWith("SEQ"))
                {
                    var annotationStr = t.Substring(t.LastIndexOf('=') + 1);
                    var precursorComposition = AminoAcidSet.GetComposition(annotationStr);
                    var peptideComposition = precursorComposition + Composition.H2O;
                    annotation = new Sequence(peptideComposition, annotationStr, AminoAcidSet);
                }
                else if (t.StartsWith("PEPMASS"))
                {
                    var p = t.Substring(t.IndexOf('=') + 1).Split(new[] {"\\s+"}, StringSplitOptions.None);
                    precursorMz = double.Parse(p[0]);
                }
            }
            return new MSMSSpectrum(charge, precursorMz, annotation, peaks);   
        }

    }
}
