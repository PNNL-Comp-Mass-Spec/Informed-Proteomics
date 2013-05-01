using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.IMSScoring;

namespace InformedProteomics.Backend.IMSTraining
{
    public class RatioScoreTrainerUsingMgfFile
    {
        private readonly List<MSMSSpectrum> _spectra;
        private readonly Dictionary<GroupParameter, List<IonType>> _ionTypes;
        private readonly Tolerance _tolerance;
        private readonly int _maxCharge;
        public Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>> RatioProbDictionary { get; private set; }

        public RatioScoreTrainerUsingMgfFile(List<MSMSSpectrum> spectra, Dictionary<GroupParameter, List<IonType>> ionTypes, Tolerance tolerance, int maxCharge)
        {
            _spectra = spectra;
            _ionTypes = ionTypes;
            _tolerance = tolerance;
            _maxCharge = maxCharge;
            RatioProbDictionary = new Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>>();
        }

        public void Train(bool isDecoy)
        {
            foreach (var spectrum in _spectra)
            {
                var annotation = isDecoy ? TrainerUsingMgfFile.GetReversedSequence(spectrum.Annotation) : spectrum.Annotation;
                for (var cutNumber = 1; cutNumber < annotation.Count; cutNumber++)
                {
                    var groupParameter = spectrum.GetGroupParameter(annotation, cutNumber);
                    var ionTypes = _ionTypes[groupParameter];
                    var explainedPeaks = spectrum.GetExplainedPeaks(annotation, cutNumber, ionTypes, _tolerance);
                    if(!RatioProbDictionary.ContainsKey(groupParameter))
                        RatioProbDictionary[groupParameter] = new Dictionary<List<IonType>, Dictionary<int, double>>();
                    var s = RatioProbDictionary[groupParameter];
                    for (var i = 0; i < ionTypes.Count; i++)
                    {
                        for (var j = 0; j < ionTypes.Count; j++)
                        {
                            //if (i == j) continue; //TODO THis is for precursor but there is not precursor ion in mgf file..
                            var ionPair = new List<IonType> {ionTypes[i], ionTypes[j]};
                            if(!s.ContainsKey(ionPair)) s[ionPair] = new Dictionary<int, double>();
                            var t = s[ionPair];
                            var ratio = FeatureEdge.GetRatio(explainedPeaks[i].Intensity, explainedPeaks[j].Intensity);
                            if (!t.ContainsKey(ratio)) t[ratio] = 0;
                            t[ratio] = t[ratio] + 1;
                        }   
                    }
                }
            }
            Normalize();
        }

        private void Normalize()
        {
            foreach (var groupParameter in GroupParameter.GetAllFragmentParameters(_maxCharge))
            {
                if(!RatioProbDictionary.ContainsKey(groupParameter))
                    RatioProbDictionary[groupParameter] = new Dictionary<List<IonType>, Dictionary<int, double>>();

                foreach (var t in RatioProbDictionary[groupParameter].Values)
                {
                    var sum = 0.0;
                    foreach (var m in t.Values)
                    {
                        sum = sum + m;
                    }
                    var keys = new List<int>(t.Keys);
                    foreach (var ratio in keys)
                    {
                        t[ratio] = Math.Max(double.MinValue, t[ratio] / sum);
                    }
                }
            }
        }
    }
}
