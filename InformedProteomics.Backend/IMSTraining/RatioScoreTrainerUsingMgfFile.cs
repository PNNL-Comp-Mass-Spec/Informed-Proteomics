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
        public Dictionary<GroupParameter, Dictionary<Tuple<IonType, IonType>, Dictionary<int, double>>> RatioProbDictionary { get; private set; }
        public Dictionary<GroupParameter, double> NoIonProbDictionary { get; private set; }
        private Dictionary<GroupParameter, double> AllNumberDictionary;

        public RatioScoreTrainerUsingMgfFile(List<MSMSSpectrum> spectra, Dictionary<GroupParameter, List<IonType>> ionTypes, Tolerance tolerance, int maxCharge)
        {
            _spectra = spectra;
            _ionTypes = ionTypes;
            _tolerance = tolerance;
            _maxCharge = maxCharge;
            RatioProbDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<IonType, IonType>, Dictionary<int, double>>>();
            NoIonProbDictionary = new Dictionary<GroupParameter, double>();
            AllNumberDictionary = new Dictionary<GroupParameter, double>();
        }

        public void Train()
        {
            foreach (var spectrum in _spectra)
            {
                var annotation = spectrum.Annotation;
                for (var cutNumber = 1; cutNumber < annotation.Count; cutNumber++)
                {
                    var groupParameter = spectrum.GetGroupParameter(annotation, cutNumber);
                    var ionTypes = _ionTypes[groupParameter];
                    var explainedPeaks = spectrum.GetExplainedPeaks(annotation, cutNumber, ionTypes, _tolerance);
                    if(!RatioProbDictionary.ContainsKey(groupParameter))
                        RatioProbDictionary[groupParameter] = new Dictionary<Tuple<IonType, IonType>, Dictionary<int, double>>();
                    if (!AllNumberDictionary.ContainsKey(groupParameter)) AllNumberDictionary[groupParameter] = 0.0;
                    if (!NoIonProbDictionary.ContainsKey(groupParameter)) NoIonProbDictionary[groupParameter] = 0.0;
                    AllNumberDictionary[groupParameter] = AllNumberDictionary[groupParameter] + 1;
                    var s = RatioProbDictionary[groupParameter];
                    /*var writeNoIon = true;
                    foreach (var peak in explainedPeaks)
                    {
                        if (peak.Intensity > 0)
                        {
                            writeNoIon = false;
                            break;
                        }
                    }*/
                    //if (writeNoIon)
                    NoIonProbDictionary[groupParameter] = NoIonProbDictionary[groupParameter] + 1;
                    for (var i = 0; i < ionTypes.Count; i++)
                    {
                        for (var j = 0; j < ionTypes.Count; j++)
                        {
                            //if (i == j) continue; //TODO THis is for precursor but there is not precursor ion in mgf file..
                            var ionPair = new Tuple<IonType, IonType>(ionTypes[i], ionTypes[j]);
                            if(!s.ContainsKey(ionPair)) s[ionPair] = new Dictionary<int, double>();
                            var t = s[ionPair];
                            //if (explainedPeaks[i].Intensity <= 0 || explainedPeaks[j].Intensity <= 0)
                            //{
                             //   continue;
                            //}
                            var ratio = FeatureEdge.GetRatioScore(explainedPeaks[i].Intensity, explainedPeaks[j].Intensity);
                            if (!t.ContainsKey(ratio)) t[ratio] = 0;
                            t[ratio] = t[ratio] + 1;
                        }   
                    }
                }
            }
            Normalize();  
            Console.WriteLine(RatioProbDictionary.Count);
        }

        private void Normalize()
        {
            foreach (var groupParameter in GroupParameter.GetAllFragmentGroupParameters(_maxCharge))
            {
                if(!RatioProbDictionary.ContainsKey(groupParameter))
                    RatioProbDictionary[groupParameter] = new Dictionary<Tuple<IonType, IonType>, Dictionary<int, double>>();
                if (!AllNumberDictionary.ContainsKey(groupParameter))
                    AllNumberDictionary[groupParameter] = 1.0;
                if (!NoIonProbDictionary.ContainsKey(groupParameter))
                    NoIonProbDictionary[groupParameter] = 1.0;

                NoIonProbDictionary[groupParameter] = NoIonProbDictionary[groupParameter]/
                                                      Math.Max(AllNumberDictionary[groupParameter], 1.0);

                var sr = RatioProbDictionary[groupParameter];

                foreach (var ionType in _ionTypes[groupParameter])
                {
                    foreach (var ionType2 in _ionTypes[groupParameter])
                    {
                        var sum = 0.0;
                        var ionTypes = new Tuple<IonType, IonType> (ionType, ionType2);
                        if(!sr.ContainsKey(ionTypes)) sr[ionTypes] = new Dictionary<int, double>();
                        var ssr = sr[ionTypes];
                        foreach (var ratio in FeatureEdge.GetAllRatioScores())
                        {
                            if (!ssr.ContainsKey(ratio)) ssr[ratio] = 1.0;
                            sum = sum + ssr[ratio];
                        }
                        var keys = new List<int> (ssr.Keys);
                        foreach (var ratio in keys)
                        {
                            ssr[ratio] = Math.Max(Double.MinValue, ssr[ratio] / sum);
                        }
                    }
                }
            }
        }
    }
}
