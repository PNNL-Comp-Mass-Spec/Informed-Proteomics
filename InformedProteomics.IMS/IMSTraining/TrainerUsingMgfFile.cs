using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.IMS.IMSTraining
{
    public class TrainerUsingMgfFile
    {
        public static void Train(string outFileName, List<MSMSSpectrum> spectra, Tolerance tolerance, int maxCharge)
        {
            var spectraTmp = new List<MSMSSpectrum>();
            foreach (var spectrum in spectra)
            {
                if (spectrum.Charge <= maxCharge) spectraTmp.Add(spectrum);
            }
            spectra = spectraTmp;
            Console.WriteLine("Number of spectra: " + spectra.Count);
            var ionTypeTrainer = new IonTypeTrainerUsingMgfFile(spectra, tolerance, maxCharge);
            ionTypeTrainer.Train();

            /*Console.WriteLine("Number of ion types: " + ionTypeTrainer.IonTypes.Count);
            foreach (var groupParameter in ionTypeTrainer.IonTypes.Keys)
            {
                Console.WriteLine(groupParameter.ToFileString());
                foreach (var ionType in ionTypeTrainer.IonTypes[groupParameter])
                {
                    Console.WriteLine("\t" + ionType);
                }
            }*/
            Console.WriteLine("Ion Type training Done");
           // return;
            var isotopeIntensityCorrelationScoreTrainerT = new IsotopeIntensityCorrelationScoreTrainerUsingMgfFile(spectra, ionTypeTrainer.IonTypes, tolerance, maxCharge);
            isotopeIntensityCorrelationScoreTrainerT.Train();
            WriteIonTypes(outFileName, ionTypeTrainer);
            Console.WriteLine("Isotope Training Done (Target)");

            var ratioScoreTrainerT = new RatioScoreTrainerUsingMgfFile(spectra, ionTypeTrainer.IonTypes, tolerance, maxCharge);
            ratioScoreTrainerT.Train();
            Write(outFileName, ionTypeTrainer, isotopeIntensityCorrelationScoreTrainerT, ratioScoreTrainerT, false);
            Console.WriteLine("Ion Ratio Training Done (Target)");

            foreach (var spectrum in spectra)
            {
                spectrum.Annotation = GetReversedSequence(spectrum.Annotation);
            }

            var isotopeIntensityCorrelationScoreTrainerD = new IsotopeIntensityCorrelationScoreTrainerUsingMgfFile(spectra, ionTypeTrainer.IonTypes, tolerance, maxCharge);
            isotopeIntensityCorrelationScoreTrainerD.Train();
            Console.WriteLine("Isotope Training Done (Decoy)");

            var ratioScoreTrainerD = new RatioScoreTrainerUsingMgfFile(spectra, ionTypeTrainer.IonTypes, tolerance, maxCharge);
            ratioScoreTrainerD.Train();
            Write(outFileName, ionTypeTrainer, isotopeIntensityCorrelationScoreTrainerD, ratioScoreTrainerD, true);
            Console.WriteLine("Ion Ratio Training Done (Decoy)");

        }

        private static void WriteIonTypes(string outFileName, IonTypeTrainerUsingMgfFile ionTypeTrainer)
        {
            var writer = new StreamWriter(outFileName);
            var ionTypes = ionTypeTrainer.IonTypes;
            writer.Write("##IONTYPES\n");
            foreach (var groupParameter in ionTypes.Keys)
            {
                writer.Write("#G\t" + groupParameter+"\n");
                foreach (var ionType in ionTypes[groupParameter])
                {
                    writer.Write("#I\t" + ionType+"\n");
                }
            }
            writer.Close();
        }

        private static void Write(string outFileName, IonTypeTrainerUsingMgfFile ionTypeTrainer, IsotopeIntensityCorrelationScoreTrainerUsingMgfFile isotopeIntensityCorrelationScoreTrainer, RatioScoreTrainerUsingMgfFile ratioScoreTrainer, bool isDecoy)
        {
            var writer = new StreamWriter(outFileName, true);
            if(isDecoy) writer.WriteLine("###DECOY");
            var isotope = isotopeIntensityCorrelationScoreTrainer.IsotopeIntensityCorrProbDictionary;
            writer.Write("##ISOTOPE\n");
            foreach (var groupParameter in isotope.Keys)
            {
                writer.Write("#G\t" + groupParameter + "\n");
                var s = isotope[groupParameter];
                var si = ionTypeTrainer.IonTypes[groupParameter];
                foreach (var ionType in s.Keys)
                {
                    writer.Write("#I\t" + si.IndexOf(ionType) + "\n");
                    var t = s[ionType];
                    foreach (var k in t.Keys)
                    {
                        writer.Write(k+","+t[k]+"\t");
                    }
                    writer.Write("\n");
                }
            }
            var ratio = ratioScoreTrainer.RatioProbDictionary;
            writer.Write("##RATIO\n");
            foreach (var groupParameter in ratio.Keys)
            {
                writer.Write("#G\t" + groupParameter + "\n");
                var s = ratio[groupParameter];
                var si = ionTypeTrainer.IonTypes[groupParameter];
                foreach (var ionTypes in s.Keys)
                {
                    writer.Write("#I");
                    writer.Write("\t" + si.IndexOf(ionTypes.Item1) + "\t" + si.IndexOf(ionTypes.Item2));

                    writer.Write("\n");
                    var t = s[ionTypes];
                    foreach (var k in t.Keys)
                    {
                        writer.Write(k + "," + t[k] + "\t");
                    }
                    writer.Write("\n");
                }
            }
            writer.Write("##NOION\n");
            var noIon = ratioScoreTrainer.NoIonProbDictionary;
            foreach (var groupParameter in noIon.Keys)
            {
                writer.Write("#G\t" + groupParameter + "\n");
                writer.Write(noIon[groupParameter]+"\n");
            }
            writer.Close();
        }

        public static Sequence GetReversedSequence(Sequence sequence)
        {
            var aaList = new List<AminoAcid>();
            for(var i=sequence.Count-2;i>=0;i--)
                aaList.Add(sequence[i]);
            aaList.Add(sequence[sequence.Count-1]); // keep C-term aa
            return new Sequence(aaList);
        }
    }
}
