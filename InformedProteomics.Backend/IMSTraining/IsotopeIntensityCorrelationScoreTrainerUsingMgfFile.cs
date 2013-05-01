using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.IMSScoring;

namespace InformedProteomics.Backend.IMSTraining
{
    public class IsotopeIntensityCorrelationScoreTrainerUsingMgfFile
    {
        private readonly List<MSMSSpectrum> _spectra;
        private readonly Dictionary<GroupParameter, List<IonType>> _ionTypes;
        private readonly Tolerance _tolerance;
        private readonly int _maxCharge;
        public Dictionary<GroupParameter, Dictionary<IonType, Dictionary<int, double>>> IsotopeIntensityCorrProbDictionary { get; private set; }

        public IsotopeIntensityCorrelationScoreTrainerUsingMgfFile(List<MSMSSpectrum> spectra, Dictionary<GroupParameter, List<IonType>> ionTypes, Tolerance tolerance, int maxCharge)
        {
            _spectra = spectra;
            _ionTypes = ionTypes;
            _tolerance = tolerance;
            _maxCharge = maxCharge;
            IsotopeIntensityCorrProbDictionary = new Dictionary<GroupParameter, Dictionary<IonType, Dictionary<int, double>>>();
        }

        public void Train(bool isDecoy)
        {
            foreach (var spectrum in _spectra)
            {
                var annotation = isDecoy ? TrainerUsingMgfFile.GetReversedSequence(spectrum.Annotation) : spectrum.Annotation;
                for (var cutNumber = 1; cutNumber < annotation.Count; cutNumber++)
                {
                    var groupParameter = spectrum.GetGroupParameter(annotation, cutNumber);
                    if (!IsotopeIntensityCorrProbDictionary.ContainsKey(groupParameter))
                        IsotopeIntensityCorrProbDictionary[groupParameter] = new Dictionary<IonType, Dictionary<int, double>>();
                    var s = IsotopeIntensityCorrProbDictionary[groupParameter];
                    foreach (var ionType in _ionTypes[groupParameter])
                    {
                        if(!s.ContainsKey(ionType)) s[ionType] = new Dictionary<int, double>();
                        var t = s[ionType];
                        var isotopomerTuple = spectrum.GetIsotopomerEnvelop(annotation, cutNumber, ionType, _tolerance);
                        var score = GetCorrelationScore(isotopomerTuple);
                        if(!t.ContainsKey(score))
                            t[score] = 0;
                        t[score] = t[score] + 1;
                    }
                }
            }
            Normalize();
        }

        private void Normalize()
        {
            foreach (var groupParameter in GroupParameter.GetAllFragmentParameters(_maxCharge))
            {
                if (!IsotopeIntensityCorrProbDictionary.ContainsKey(groupParameter))
                    IsotopeIntensityCorrProbDictionary[groupParameter] = new Dictionary<IonType, Dictionary<int, double>>();

                foreach (var t in IsotopeIntensityCorrProbDictionary[groupParameter].Values)
                {
                    var sum = 0.0;
                    foreach (var m in t.Values)
                    {
                        sum = sum + m;
                    }
                    var keys = new List<int>(t.Keys);
                    foreach (var score in keys)
                    {
                        t[score] = Math.Max(double.MinValue, t[score]/sum);
                    }
                }
            }
        }

        private int GetCorrelationScore(Tuple<List<MSMSSpectrumPeak>, float[]> isotopomerTuple)
        {
            var observed = isotopomerTuple.Item1;
            var theoretical = isotopomerTuple.Item2;
            var maxIndex = 0;
            var maxAbundance = float.NegativeInfinity;
            for (var i = 0; i < theoretical.Length; i++)
            {
                if (!(maxAbundance < theoretical[i])) continue;
                maxAbundance = theoretical[i];
                maxIndex = i;
            }
            var c1 = new double[FeatureNode.NumSupport];
            var c2 = new double[FeatureNode.NumSupport];

            for (var i = 0; i < c1.Length; i++)
            {
                var i1 = maxIndex + i; 
                var i2 = -FeatureNode.NumMinusIsotope + maxIndex + i;
                c1[i] = observed[Math.Min(i1, observed.Count -1)].Intensity;
                c2[i] = i2 < 0 ? 0 : theoretical[Math.Min(i2, theoretical.Length-1)];
            }
            return SubScoreFactory.CorrToInt(StatisticsTools.GetCorrelation(c1, c2));
        }
    }
}
