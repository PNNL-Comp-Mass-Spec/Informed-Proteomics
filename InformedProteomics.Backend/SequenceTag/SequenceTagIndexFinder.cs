using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.SequenceTag
{
    public class SequenceTagIndexFinder
    {

        public SequenceTagIndexFinder(Tolerance tolerance, int minCharge, int maxCharge)
        {
            _tolerance = tolerance;
            _minCharge = minCharge;
            _maxCharge = maxCharge;

        }

        public Tuple<int, int, int, int> GetLongestSequence(ProductSpectrum spectrum, Sequence sequence)
        {
            _spectrum = spectrum;
            _sequence = sequence;
            _baseIonTypes = _spectrum.ActivationMethod != ActivationMethod.ETD ? BaseIonTypesCid : BaseIonTypesEtd;

            var cleavages = _sequence.GetInternalCleavages().ToArray();
            var prefixValueArr = new int[cleavages.Length];
            var suffixValueArr = new int[cleavages.Length];

            int cleavageIndex = 0;

            foreach (var c in cleavages)
            {
                foreach (var baseIonType in _baseIonTypes)
                {
                    var fragmentComposition = baseIonType.IsPrefix
                                 ? c.PrefixComposition + baseIonType.OffsetComposition
                                 : c.SuffixComposition + baseIonType.OffsetComposition;
                    for (int charge = _minCharge; charge <= _maxCharge; charge++)
                    {

                        var ion = new Ion(fragmentComposition, charge);
                        int baseIsotopePeakIndex;
                        int nIsotopes;
                        int nMatchedIsotopes;

                        if (FindIon(ion, _tolerance, RelativeIsotopeIntensityThreshold, out baseIsotopePeakIndex, out nIsotopes, out nMatchedIsotopes))
                        {
                            if (baseIonType.IsPrefix) prefixValueArr[cleavageIndex] = 1;
                            else suffixValueArr[cleavageIndex] = 1;

                        }
                    }
                }
                cleavageIndex++;
            }

            var sequenceArr = new int[_sequence.Count];
            sequenceArr[0] = prefixValueArr[0];
            sequenceArr[sequenceArr.Length - 1] = suffixValueArr[suffixValueArr.Length - 1];

            /*for (int i = 1; i < prefixValueArr.Length ; i++)
            {
                if(prefixValueArr[i] == 1 && prefixValueArr[i - 1] == 1) sequenceArr[i] =  1;
            }

            for (int i = suffixValueArr.Length - 2; i >= 0; i--)
            {
                if (suffixValueArr[i] == 1 && suffixValueArr[i + 1] == 1) sequenceArr[i + 1] = 1;
            }*/


            var prefixString = FindLongestSubstring(prefixValueArr);
            var suffixString = FindLongestSubstring(suffixValueArr);
            var prefixIndex = string.Concat(prefixValueArr);
            var prefixStartIndex = prefixIndex.IndexOf(prefixString) + 1;
            var suffixIndex = string.Concat(suffixValueArr);
            var suffixStartIndex = suffixIndex.IndexOf(suffixString) + 1;
            return new Tuple<int, int, int, int>(prefixStartIndex, prefixStartIndex + prefixString.Length - 1, suffixStartIndex, suffixStartIndex + suffixString.Length - 1);
        }


        private readonly Tolerance _tolerance;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private BaseIonType[] _baseIonTypes;

        private const double RelativeIsotopeIntensityThreshold = 0.8;
        public static readonly BaseIonType[] BaseIonTypesCid, BaseIonTypesEtd;

        private ProductSpectrum _spectrum;
        private Sequence _sequence;


        private string FindLongestSubstring(int[] booleanSequence)
        {
            var sequenceString = "";
            foreach (var s in booleanSequence)
            {
                sequenceString += s;
            }

            Regex regex = new Regex(@"[1]+");
            var matches = regex.Matches(sequenceString);
            var largestSubstring = "";
            foreach (Match m in matches)
            {
                if (m.Length > largestSubstring.Length) largestSubstring = m.Value;
            }


            return largestSubstring;
        }

        private bool FindIon(Ion ion, Tolerance tolerance, double relativeIntensityThreshold, out int baseIsotopePeakIndex, out int nIsotopes, out int nMatchedIsotopes)
        {

            var corrScore = _spectrum.GetCorrScore(ion, tolerance, relativeIntensityThreshold);

            ///isotopePeaks
            /// 
            //FitScoreCalculator.GetPearsonCorrelation()
            

            
            //matchedPeakIndex = new List<int>();
            var baseIsotopeIndex = ion.Composition.GetMostAbundantIsotopeZeroBasedIndex();
            var isotopomerEnvelope = ion.Composition.GetIsotopomerEnvelopeRelativeIntensities();
            var baseIsotopMz = ion.GetIsotopeMz(baseIsotopeIndex);
            baseIsotopePeakIndex = _spectrum.FindPeakIndex(baseIsotopMz, tolerance);

            nIsotopes = isotopomerEnvelope.Select(x => x >= relativeIntensityThreshold).Count();
            nMatchedIsotopes = 0;

            if (baseIsotopePeakIndex < 0) return false;
            //if (baseIsotopePeakIndex < 0) baseIsotopePeakIndex = ~baseIsotopePeakIndex;
            nMatchedIsotopes++;

            // go down
            var peakIndex = baseIsotopePeakIndex;
            //matchedPeakIndex.Add(peakIndex);
            for (var isotopeIndex = baseIsotopeIndex - 1; isotopeIndex >= 0; isotopeIndex--)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;

                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex - 1; i >= 0; i--)
                {
                    var peakMz = _spectrum.Peaks[i].Mz;
                    if (peakMz < minMz)
                    {
                        //peakIndex = i;
                        //break;
                        return false;
                    }
                    if (peakMz <= maxMz)    // find match, move to prev isotope
                    {
                        peakIndex = i;
                        //matchedPeakIndex.Add(peakIndex);
                        nMatchedIsotopes++;
                        break;
                    }
                }
            }

            // go up
            peakIndex = baseIsotopePeakIndex;
            for (var isotopeIndex = baseIsotopeIndex + 1; isotopeIndex < isotopomerEnvelope.Length; isotopeIndex++)
            {
                if (isotopomerEnvelope[isotopeIndex] < relativeIntensityThreshold) break;

                var isotopeMz = ion.GetIsotopeMz(isotopeIndex);
                var tolTh = tolerance.GetToleranceAsTh(isotopeMz);
                var minMz = isotopeMz - tolTh;
                var maxMz = isotopeMz + tolTh;
                for (var i = peakIndex + 1; i < _spectrum.Peaks.Length; i++)
                {
                    var peakMz = _spectrum.Peaks[i].Mz;
                    if (peakMz > maxMz)
                    {
                        //peakIndex = i;
                        //break;
                        return false;
                    }
                    if (peakMz >= minMz)    // find match, move to prev isotope
                    {
                        peakIndex = i;
                        //matchedPeakIndex.Add(peakIndex);
                        nMatchedIsotopes++;
                        break;
                    }
                }
            }
            return true;
        }

        static SequenceTagIndexFinder()
        {
            BaseIonTypesCid = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesEtd = new[] { BaseIonType.C, BaseIonType.Z };
        }
    }
}
