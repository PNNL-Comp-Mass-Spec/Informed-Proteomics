using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.IMSTraining
{
    public class TrainerUsingMgfFile
    {
        public TrainerUsingMgfFile(string mgfFileName, string outFileName, Tolerance tolerance)
        {
            Train(outFileName, new MgfParser(mgfFileName).Spectra, tolerance);
        }

        public static void Train(string outFileName, List<MSMSSpectrum> spectra, Tolerance tolerance)
        {
            var ionTypeTrainer = new IonTypeTrainerUsingMgfFile(spectra, tolerance);
            ionTypeTrainer.Train();
            var isotopeIntensityCorrelationScoreTrainer = new IsotopeIntensityCorrelationScoreTrainerUsingMgfFile(spectra, ionTypeTrainer.IonTypes, tolerance);
            isotopeIntensityCorrelationScoreTrainer.Train(false);
            var ratioScoreTrainer = new RatioScoreTrainerUsingMgfFile(spectra, ionTypeTrainer.IonTypes, tolerance);
            ratioScoreTrainer.Train(false);
            Write(outFileName, ionTypeTrainer, isotopeIntensityCorrelationScoreTrainer, ratioScoreTrainer);

            isotopeIntensityCorrelationScoreTrainer = new IsotopeIntensityCorrelationScoreTrainerUsingMgfFile(spectra, ionTypeTrainer.IonTypes, tolerance);
            isotopeIntensityCorrelationScoreTrainer.Train(true);
            ratioScoreTrainer = new RatioScoreTrainerUsingMgfFile(spectra, ionTypeTrainer.IonTypes, tolerance);
            ratioScoreTrainer.Train(true);
            Write(outFileName, isotopeIntensityCorrelationScoreTrainer, ratioScoreTrainer);
        }

        private static void Write(string outfileName, IonTypeTrainerUsingMgfFile ionTypeTrainer, IsotopeIntensityCorrelationScoreTrainerUsingMgfFile isotopeIntensityCorrelationScoreTrainer, RatioScoreTrainerUsingMgfFile ratioScoreTrainer)
        {
            
        }

        private static void Write(string outfileName, IsotopeIntensityCorrelationScoreTrainerUsingMgfFile isotopeIntensityCorrelationScoreTrainer, RatioScoreTrainerUsingMgfFile ratioScoreTrainer)
        {

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
