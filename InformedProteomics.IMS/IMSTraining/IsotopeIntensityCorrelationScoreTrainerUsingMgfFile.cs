using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;
using InformedProteomics.IMS.IMSScoring;

namespace InformedProteomics.IMS.IMSTraining
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

        public void Train()
        {
            foreach (var spectrum in _spectra)
            {
                var annotation = spectrum.Annotation;
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
                        if (isotopomerTuple == null) continue;
                        var score = GetCorrelationScore(isotopomerTuple);
                        if(!t.ContainsKey(score))
                            t[score] = 0;
                        t[score] = t[score] + 1;
                    }
                }
            }
            Normalize();
            Console.WriteLine(IsotopeIntensityCorrProbDictionary.Count);
        }

        private void Normalize()
        {
            foreach (var groupParameter in GroupParameter.GetAllFragmentGroupParameters(_maxCharge))
            {
                if (!IsotopeIntensityCorrProbDictionary.ContainsKey(groupParameter))
                    IsotopeIntensityCorrProbDictionary[groupParameter] = new Dictionary<IonType, Dictionary<int, double>>();
                var si = IsotopeIntensityCorrProbDictionary[groupParameter];
                foreach (var ionType in _ionTypes[groupParameter])
                {
                    if(!si.ContainsKey(ionType)) si[ionType] = new Dictionary<int, double>();
                    var ssi = si[ionType];
                    var sum = 0.0;
                    foreach (var score in SubScoreFactory.GetAllCorrIntScore())
                    {
                        if (!ssi.ContainsKey(score)) ssi[score] = 1.0;
                        sum = sum + ssi[score];
                    }
                    var keys = new List<int>(ssi.Keys);
                    foreach (var score in keys)
                    {
                        ssi[score] = Math.Max(double.MinValue,ssi[score] / sum);
                    }
                }
            }
        }

        private int GetCorrelationScore(Tuple<List<MSMSSpectrumPeak>, double[]> isotopomerTuple)
        {
            var observed = isotopomerTuple.Item1;
            var theoretical = isotopomerTuple.Item2;
            var maxIndex = 0;
            var maxAbundance = double.NegativeInfinity;
            for (var i = 0; i < theoretical.Length; i++)
            {
                if (!(maxAbundance < theoretical[i])) continue;
                maxAbundance = theoretical[i];
                maxIndex = i;
            }
            var c1 = new double[theoretical.Length + FeatureNode.OffsetFromMonoIsotope];
            var c2 = new double[theoretical.Length + FeatureNode.OffsetFromMonoIsotope];

            for (var i = 0; i < c1.Length; i++)
            {
                var i1 = maxIndex + i; 
                var i2 = -FeatureNode.OffsetFromMonoIsotope + maxIndex + i;
                c1[i] = observed[Math.Min(i1, observed.Count -1)].Intensity;
                c2[i] = i2 < 0 ? 0 : theoretical[Math.Min(i2, theoretical.Length-1)];
            }
        //    Console.WriteLine(c1.Length);
            return SubScoreFactory.CorrToIntForIsotopeScore(SimpleMath.GetCorrelation(c1, c2));
        }
    }
}
