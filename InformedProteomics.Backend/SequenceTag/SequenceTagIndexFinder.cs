using System;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using System.Security.Cryptography;
using System.Text;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
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

        public Tuple<int, int, int, int, string, string> GetLongestSequence(ProductSpectrum spectrum, Sequence sequence)
        {
            _spectrum = spectrum;
            _sequence = sequence;
            _baseIonTypes = _spectrum.ActivationMethod != ActivationMethod.ETD ? BaseIonTypesCid : BaseIonTypesEtd;

            var cleavages = _sequence.GetInternalCleavages().ToArray();
            var prefixValueArr = new int[cleavages.Length];
            var suffixValueArr = new int[cleavages.Length];
            var prefixPeakArr = new Peak[cleavages.Length];
            var suffixPeakArr = new Peak[cleavages.Length];
            //var peakList = new double[_spectrum.Peaks.Length];

            int cleavageIndex = 0;

            /*
            for (int i = 0; i < peakList.Length; i++)
            {
                peakList[i] = _spectrum.Peaks[i].Intensity;
            }*/

            //var rankings = ArrayUtil.GetRankings(peakList);

            foreach (var c in cleavages)
            {
                foreach (var baseIonType in _baseIonTypes)
                {
                    var fragmentComposition = baseIonType.IsPrefix
                                 ? c.PrefixComposition + baseIonType.OffsetComposition
                                 : c.SuffixComposition + baseIonType.OffsetComposition;
                    for (var charge = _minCharge; charge <= _maxCharge; charge++)
                    {

                        var ion = new Ion(fragmentComposition, charge);
                        if (_spectrum.GetCorrScore(ion, _tolerance, RelativeIsotopeIntensityThreshold) < .7) continue;
                        if (baseIonType.IsPrefix) prefixValueArr[cleavageIndex] = 1;
                        else suffixValueArr[cleavageIndex] = 1;


                    }
                }
                cleavageIndex++;
            }

            var prefixSequenceArr = new int[_sequence.Count];
            var suffixSequenceArr = new int[_sequence.Count];
            prefixSequenceArr[0] = prefixValueArr[0];
            suffixSequenceArr[suffixSequenceArr.Length - 1] = suffixValueArr[suffixValueArr.Length - 1];

            for (int i = 1; i < prefixValueArr.Length; i++)
            {
                if (prefixValueArr[i] == 1 && prefixValueArr[i - 1] == 1)
                {
                    if (_sequence[i] is ModifiedAminoAcid) continue;
                    prefixSequenceArr[i] = 1;

                }
            }

            for (int i = suffixValueArr.Length - 2; i >= 0; i--)
            {
                if (suffixValueArr[i] == 1 && suffixValueArr[i + 1] == 1)
                {
                    if (_sequence[i + 1] is ModifiedAminoAcid) continue;
                    suffixSequenceArr[i + 1] = 1;
                }
            }



            var prefixSubString = FindLongestSubstring(prefixSequenceArr);
            var prefixStartIndex = -1;
            var prefixEndIndex = -1;
            //var prefixSequencePeaks = new List<Peak>();
            //var prefixPval = -1.0;
            var prefixSequence = "";
            if (prefixSubString != "")
            {
                var prefixIndex = string.Concat(prefixSequenceArr);
                prefixStartIndex = prefixIndex.IndexOf(prefixSubString) + 1;
                prefixEndIndex = (prefixStartIndex == 1) ? 1 : prefixStartIndex + prefixSubString.Length - 1;
                //prefixSequencePeaks = GetPrefixSequencePeaks(prefixPeakArr, prefixStartIndex, prefixEndIndex);
                //var prefixRankSum = GetSequenceRankSum(prefixSequencePeaks, rankings, peakList);
                //prefixPval = FitScoreCalculator.GetRankSumPvalue(peakList.Length, prefixSequencePeaks.Count, prefixRankSum);
                prefixSequence = GetStringSubSequence(_sequence, prefixStartIndex, prefixEndIndex);
            }

            var suffixSubString = FindLongestSubstring(suffixSequenceArr);
            var suffixStartIndex = -1;
            var suffixEndIndex = -1;
            //var suffixSequencePeaks = new List<Peak>();
            //var suffixPval = -1.0;
            var suffixSequence = "";

            if (suffixSubString != "")
            {
                var suffixIndex = string.Concat(suffixSequenceArr);
                suffixStartIndex = suffixIndex.IndexOf(suffixSubString) + 1;
                suffixEndIndex = (suffixStartIndex == 1) ? 1 : suffixStartIndex + suffixSubString.Length - 1;
                //suffixSequencePeaks = GetSuffixSequencePeaks(suffixPeakArr, suffixStartIndex, suffixEndIndex);
                //var suffixRankSum = GetSequenceRankSum(suffixSequencePeaks, rankings, peakList);
                //suffixPval = FitScoreCalculator.GetRankSumPvalue(peakList.Length, suffixSequencePeaks.Count, suffixRankSum);
                suffixSequence = GetStringSubSequence(_sequence, suffixStartIndex, suffixEndIndex);
            }

            return new Tuple<int, int, int, int, string, string>(prefixStartIndex, prefixEndIndex, suffixStartIndex, suffixEndIndex, prefixSequence, suffixSequence);
        }


        private readonly Tolerance _tolerance;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private BaseIonType[] _baseIonTypes;

        private const double RelativeIsotopeIntensityThreshold = 0.8;
        public static readonly BaseIonType[] BaseIonTypesCid, BaseIonTypesEtd;

        private ProductSpectrum _spectrum;
        private Sequence _sequence;


        /*private double GetSequenceRankSum(List<Peak> seqPeaks, int[] rankArr, double[] peakInt)
        {
            var rankSum = 0.0;
            foreach (var p in seqPeaks)
            {
                var i = Array.IndexOf(peakInt, p.Intensity);
                rankSum += rankArr[i];
            }
            return rankSum;
        }*/

        private string GetStringSubSequence(Sequence seq, int start, int end)
        {
            var startIndex = start - 1;
            var endIndex = end - 1;
            var subString = "";
            if (startIndex == endIndex)
            {
                subString += seq[start].Residue;
                return subString;
            }
            for (int i = startIndex; i < endIndex + 1; i++)
            {
                subString += seq[i].Residue;
            }
            return subString;
        }

        /*private List<Peak> GetPrefixSequencePeaks(Peak[] prefixPeaks, int start, int end)
        {
            var startIndex = start - 1;
            var endIndex = end - 1;
            var peaks = new List<Peak>();
           
            for (var i = startIndex; i <= endIndex; i++)
            {
                if (i > prefixPeaks.Length) continue;
                if (i == startIndex && startIndex != 0)
                {
                    peaks.Add(prefixPeaks[i - 1]);
                    peaks.Add(prefixPeaks[i]);
                }
                else peaks.Add(prefixPeaks[i]);
            }
            return peaks;
        }*/

        /*private List<Peak> GetSuffixSequencePeaks(Peak[] suffixPeaks, int start, int end)
        {
            var startIndex = start - 1;
            var endIndex = end - 1;
            var peaks = new List<Peak>();
            
            for (var i = startIndex; i <= endIndex; i++)
            {
                if (i > suffixPeaks.Length-1) continue;
                if (i == startIndex && startIndex != 0)
                {
                    peaks.Add(suffixPeaks[i - 1]);
                    peaks.Add(suffixPeaks[i]);
                }
                else peaks.Add(suffixPeaks[i]);
            }
            return peaks;
        }*/

        private string FindLongestSubstring(int[] booleanSequence)
        {
            var sequenceString = String.Concat(booleanSequence);
            Regex regex = new Regex(@"[1]+");
            var matches = regex.Matches(sequenceString);
            var largestSubstring = "";
            foreach (Match m in matches)
            {
                if (m.Length > largestSubstring.Length) largestSubstring = m.Value;
            }

            return largestSubstring;
        }

        static SequenceTagIndexFinder()
        {
            BaseIonTypesCid = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesEtd = new[] { BaseIonType.C, BaseIonType.Z };
        }
    }
}
