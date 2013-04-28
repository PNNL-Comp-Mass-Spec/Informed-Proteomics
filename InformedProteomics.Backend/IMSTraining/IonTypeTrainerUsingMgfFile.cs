using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
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
            var allKnownIonTypes = new IonTypeFactory().GetAllKnownIonTypes();
            foreach (var spectrum in _spectra)
            {
                var annotation = spectrum.Annotation;
                if (annotation == null) continue;

                var precursorIon = new Ion(annotation.Composition + Composition.H2O, spectrum.Charge);
                for (var i = 1; i < annotation.Count; i++)
                {
                    var prefixComposition = annotation.GetComposition(0, i);
                    var suffixCompostion = annotation.Composition - prefixComposition;
                    var groupParameter = new GroupParameter(prefixComposition, annotation[i - 1].Residue, annotation[i].Residue, precursorIon);
                    if (!_ionFrequencyFunction.ContainsKey(groupParameter))
                        _ionFrequencyFunction[groupParameter] = new Dictionary<IonType, double>();
                    var subIonFrequencyFunction = _ionFrequencyFunction[groupParameter];
                    // ReSharper disable PossibleMultipleEnumeration
                    foreach (var ionType in allKnownIonTypes)
                    // ReSharper restore PossibleMultipleEnumeration
                    {
                        var composition = ionType.IsPrefixIon ? prefixComposition : suffixCompostion;
                        var mz = ionType.GetMz(composition);
                        if (spectrum.GetPeaks(mz, _tolerance).Count <= 0) continue;
                        if (!subIonFrequencyFunction.ContainsKey(ionType)) subIonFrequencyFunction[ionType] = 0;
                        subIonFrequencyFunction[ionType]++;
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
