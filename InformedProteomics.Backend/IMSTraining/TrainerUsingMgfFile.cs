using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.IMSTraining
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
            isotopeIntensityCorrelationScoreTrainerT.Train(false);
            Console.WriteLine("Isotope Training Done (Target)");

            var ratioScoreTrainerT = new RatioScoreTrainerUsingMgfFile(spectra, ionTypeTrainer.IonTypes, tolerance, maxCharge);
            ratioScoreTrainerT.Train(false);
            Write(outFileName, ionTypeTrainer, isotopeIntensityCorrelationScoreTrainerT, ratioScoreTrainerT, false);
            Console.WriteLine("Ion Ratio Training Done (Target)");

            var isotopeIntensityCorrelationScoreTrainerD = new IsotopeIntensityCorrelationScoreTrainerUsingMgfFile(spectra, ionTypeTrainer.IonTypes, tolerance, maxCharge);
            isotopeIntensityCorrelationScoreTrainerD.Train(true);
            Console.WriteLine("Isotope Training Done (Decoy)");

            var ratioScoreTrainerD = new RatioScoreTrainerUsingMgfFile(spectra, ionTypeTrainer.IonTypes, tolerance, maxCharge);
            ratioScoreTrainerD.Train(true);
            Write(outFileName, isotopeIntensityCorrelationScoreTrainerD, ratioScoreTrainerD, true);
            Console.WriteLine("Ion Ratio Training Done (Decoy)");
            
        }

        private static void Write(string outFileName, IonTypeTrainerUsingMgfFile ionTypeTrainer, IsotopeIntensityCorrelationScoreTrainerUsingMgfFile isotopeIntensityCorrelationScoreTrainer, RatioScoreTrainerUsingMgfFile ratioScoreTrainer, bool isDecoy)
        {
            var writer = new StreamWriter(outFileName);
            var ionTypes = ionTypeTrainer.IonTypes;
            writer.Write("#IONTYPES\n");
            foreach (var groupParameter in ionTypes.Keys)
            {
                writer.Write("#GROUP\t" + groupParameter.ToFileString()+"\n");
                foreach (var ionType in ionTypes[groupParameter])
                {
                    writer.Write("#IONTYPE\t" + ionType+"\n");
                }
            }
            writer.Close();
            Write(outFileName, isotopeIntensityCorrelationScoreTrainer, ratioScoreTrainer, isDecoy);
        }

        private static void Write(string outFileName, IsotopeIntensityCorrelationScoreTrainerUsingMgfFile isotopeIntensityCorrelationScoreTrainer, RatioScoreTrainerUsingMgfFile ratioScoreTrainer, bool isDecoy)
        {
            var writer = new StreamWriter(outFileName, true);
            if(isDecoy) writer.WriteLine("###DECOY");
            var isotope = isotopeIntensityCorrelationScoreTrainer.IsotopeIntensityCorrProbDictionary; 
            writer.Write("#ISOTOPE\n");
            foreach (var groupParameter in isotope.Keys)
            {
                writer.Write("#GROUP\t" + groupParameter.ToFileString() + "\n");
                var s = isotope[groupParameter];
                foreach (var ionType in s.Keys)
                {
                    writer.Write("IONTYPE\t" + ionType + "\n");
                    var t = s[ionType];
                    foreach (var k in t.Keys)
                    {
                        writer.Write(k+","+t[k]+"\t");
                    }
                    writer.Write("\n");
                }
            }
            var ratio = ratioScoreTrainer.RatioProbDictionary;
            writer.Write("#RATIO\n");
            foreach (var groupParameter in ratio.Keys)
            {
                writer.Write("#GROUP\t" + groupParameter.ToFileString() + "\n");
                var s = ratio[groupParameter];
                foreach (var ionTypes in s.Keys)
                {
                    writer.Write("IONTYPES");
                    foreach (var ionType in ionTypes)
                    {
                        writer.Write("\t" + ionType);
                    }    
                    writer.Write("\n");
                    var t = s[ionTypes];
                    foreach (var k in t.Keys)
                    {
                        writer.Write(k + "," + t[k] + "\t");
                    }
                    writer.Write("\n");
                }
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
