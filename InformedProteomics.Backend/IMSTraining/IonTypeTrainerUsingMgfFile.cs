using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.IMSScoring;

namespace InformedProteomics.Backend.IMSTraining
{
    public class IonTypeTrainerUsingMgfFile
    {
        private const int NumberPerGroup = 6;
        private readonly List<MSMSSpectrum> _spectra;
        private readonly Tolerance _tolerance;
        private readonly Dictionary<GroupParameter, Dictionary<IonType, double>> _ionFrequencyFunction;
        public Dictionary<GroupParameter, List<IonType>> IonTypes { get; private set; } 

        public IonTypeTrainerUsingMgfFile(List<MSMSSpectrum> spectra, Tolerance tolerance)
        {
            _spectra = spectra;
            _tolerance = tolerance;
            _ionFrequencyFunction = new Dictionary<GroupParameter, Dictionary<IonType, double>>();
            IonTypes = new Dictionary<GroupParameter, List<IonType>>();

        }

        public void Train()
        {
            var allKnownIonTypes = new List<IonType>(new IonTypeFactory().GetAllKnownIonTypes());
            foreach (var spectrum in _spectra)
            {
                var annotation = spectrum.Annotation;
                for (var cutNumber = 1; cutNumber < annotation.Count; cutNumber++)
                {
                    var groupParameter = spectrum.GetGroupParameter(annotation, cutNumber);
                    
                    if (!_ionFrequencyFunction.ContainsKey(groupParameter))
                        _ionFrequencyFunction[groupParameter] = new Dictionary<IonType, double>();
                    var subIonFrequencyFunction = _ionFrequencyFunction[groupParameter];
                    var ionTypes = spectrum.GetExplainingIonTypes(annotation, cutNumber, allKnownIonTypes, _tolerance);
                    foreach (var ionType in ionTypes)
                    {
                        if (!subIonFrequencyFunction.ContainsKey(ionType)) subIonFrequencyFunction[ionType] = 0;
                        subIonFrequencyFunction[ionType] = subIonFrequencyFunction[ionType] + 1;
                    }
                }
            }
            GetIonTypes();
        }

        private void GetIonTypes()
        {
            foreach (var groupParameter in _ionFrequencyFunction.Keys)
            {
                IonTypes[groupParameter] = new List<IonType>();
                var subIonTypes = IonTypes[groupParameter];
                var subIonFrequencyFunction = _ionFrequencyFunction[groupParameter];
                var offsets = new List<double>(subIonFrequencyFunction.Values);
                offsets.Sort();//ascending
                foreach (var ionType in subIonFrequencyFunction.Keys)
                {
                    if (!(subIonFrequencyFunction[ionType] >= offsets[offsets.Count - NumberPerGroup])) continue;
                    if (subIonTypes.Count < NumberPerGroup) 
                        subIonTypes.Add(ionType);
                }
            }
        }
    }
}
