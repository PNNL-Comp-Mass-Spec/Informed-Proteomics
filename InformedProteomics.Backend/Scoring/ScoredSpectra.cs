using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Science;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Scoring
{
    public class ScoredSpectra
    {
        public float Ms1ScoreThreshold { get; set; }

        public ScoredSpectra(Spectra spectra)
        {
            Spectra = spectra;
        }

        public Spectra Spectra
        {
            get;
            private set;
        }

        public float GetProductIonScore(Composition productComp, Composition precursorComp)
        {
            throw new System.NotImplementedException();
        }

        public float GetPrecursorIonScore(Composition precursorComp)
        {
            throw new System.NotImplementedException();
        }

        public Feature[] GetMatchedPrecursors(Composition precursorComp)
        {
            throw new System.NotImplementedException();
        }

        public IList<Feature> GetMatchedPrecursors(Ion ion)
        {
            IList<Feature> matchedPrecursorList = new List<Feature>();

            // enumerate MS1 spectra
            foreach(Spectrum spec in Spectra.GetSpectra(1))
            {
                IList<Peak> matchedPeaks = spec.GetMatchedPeaks(ion);
                String specId = spec.GetId();
//                Spectra.GetProductIonSpectra()
                float ms1Score = IsotopomerEnvScorer.GetScore(ion, matchedPeaks);
                if(ms1Score > Ms1ScoreThreshold)
                {
                    matchedPrecursorList.Add(null);
                }
            }
            return null;
        }

    }
}
