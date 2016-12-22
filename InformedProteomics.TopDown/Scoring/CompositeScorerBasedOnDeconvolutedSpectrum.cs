using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.TopDown;

namespace InformedProteomics.TopDown.Scoring
{
    using InformedProteomics.Backend.Data.Sequence;

    public class CompositeScorerBasedOnDeconvolutedSpectrum : CompositeScorer
    {
        public CompositeScorerBasedOnDeconvolutedSpectrum(DeconvolutedSpectrum deconvolutedSpectrum, ProductSpectrum spec, Tolerance productTolerance, IMassBinning comparer, ActivationMethod actMethod = ActivationMethod.Unknown)
            : base(deconvolutedSpectrum, productTolerance, activationMethod: actMethod)
        {
            ReferencePeakIntensity = GetRefIntensity(spec.Peaks);
            _comparer = comparer;
            _massBinToPeakMap = new Dictionary<int, DeconvolutedPeak>();

            foreach (var p in deconvolutedSpectrum.Peaks)
            {
                var mass = p.Mz;
                var deltaMass = productTolerance.GetToleranceAsDa(mass, 1);
                var minMass = mass - deltaMass;
                var maxMass = mass + deltaMass;

                var binNum = comparer.GetBinNumber(mass);

                if (binNum < 0)
                {
                    binNum = comparer.GetBinNumber(minMass);
                    if (binNum < 0) binNum = comparer.GetBinNumber(maxMass);
                }

                // filter out
                if (binNum < 0) continue;

                UpdateDeconvPeak(binNum, p as DeconvolutedPeak);
                // going up
                for (var nextBinNum = binNum + 1; nextBinNum < comparer.NumberOfBins; nextBinNum++)
                {
                    var nextBinMass = comparer.GetMassStart(nextBinNum);
                    if (minMass < nextBinMass && nextBinMass < maxMass) UpdateDeconvPeak(nextBinNum, p as DeconvolutedPeak); //_ionMassChkBins[nextBinNum] = true;
                    else break;
                }

                // going down
                for (var prevBinNum = binNum - 1; prevBinNum < comparer.NumberOfBins; prevBinNum--)
                {
                    var prevBinMass = comparer.GetMassEnd(prevBinNum);
                    if (minMass < prevBinMass && prevBinMass < maxMass) UpdateDeconvPeak(prevBinNum, p as DeconvolutedPeak); //_ionMassChkBins[prevBinNum] = true;
                    else break;
                }
            }
        }

        private void UpdateDeconvPeak(int binNum, DeconvolutedPeak newPeak)
        {
            if (newPeak == null) return;
            DeconvolutedPeak existingPeak;
            if (_massBinToPeakMap.TryGetValue(binNum, out existingPeak))
            {
                if (existingPeak.Intensity < newPeak.Intensity) _massBinToPeakMap[binNum] = newPeak;
            }
            else
            {
                _massBinToPeakMap[binNum] = newPeak;
            }
        }

        public override double GetFragmentScore(Composition prefixFragmentComposition, Composition suffixFragmentComposition,
            AminoAcid nTerminalResidue = null,
            AminoAcid cTerminalResidue = null)
        {
            var score = 0.0;
            var prefixHit = false;
            var suffixHit = false;

            var ionsFound = new Dictionary<BaseIonType, double>();
            foreach (var baseIonType in BaseIonTypes)
            {
                var fragmentComposition = baseIonType.IsPrefix
                    ? prefixFragmentComposition + baseIonType.OffsetComposition
                    : suffixFragmentComposition + baseIonType.OffsetComposition;

                var param = baseIonType.IsPrefix ? ScoreParam.Prefix : ScoreParam.Suffix;

                var massBinNum = _comparer.GetBinNumber(fragmentComposition.Mass);
                if (massBinNum >= 0 && massBinNum < _comparer.NumberOfBins)
                {
                    DeconvolutedPeak existingPeak;
                    if (_massBinToPeakMap.TryGetValue(massBinNum, out existingPeak))
                    {
                        var massErrorPpm = 1e6 * (Math.Abs(existingPeak.Mass - fragmentComposition.Mass)/fragmentComposition.Mass);

                        double ionscore = 0;
                        ionscore += param.Count;
                        ionscore += param.Intensity * Math.Min(existingPeak.Intensity / ReferencePeakIntensity, 1.0); // intensity-based scoring
                        ionscore += param.Dist * existingPeak.Dist; // Envelope distance-based scoring
                        ionscore += param.Corr * existingPeak.Corr; // Envelope correlation-based scoring
                        ionscore += param.MassError * massErrorPpm; // Envelope correlation-based scoring

                        if (ionsFound.ContainsKey(baseIonType)) ionsFound.Add(baseIonType, ionscore);
                        if (baseIonType == BaseIonType.Ar && ionsFound.ContainsKey(BaseIonType.A) && ionscore < ionsFound[BaseIonType.A]) continue;
                        if (baseIonType == BaseIonType.YM1 && ionsFound.ContainsKey(BaseIonType.Y) && ionscore < ionsFound[BaseIonType.Y]) continue;

                        score += ionscore;

                        if (baseIonType.IsPrefix) prefixHit = true;
                        else suffixHit = true;
                    }
                }
            }

            if (prefixHit && suffixHit) score += ScoreParam.ComplementaryIonCount;
            return score;
        }

        public double?[][] GetNodeScores(double proteinMass)
        {


            var numNodes = _comparer.GetBinNumber(proteinMass) + 1;
            var nodeScores = new double?[2][];
            var prefixFragScores = nodeScores[0] = new double?[numNodes];
            var suffixFragScores = nodeScores[1] = new double?[numNodes];
            var deconvPeaks = (DeconvolutedPeak[])Ms2Spectrum.Peaks;

            // assume that peaks are prefixFragment ions
            foreach (var peak in deconvPeaks)
            {
                var prefixIonMass = peak.Mass;
                var prefixFragmentMass = prefixIonMass - PrefixOffsetMass;
                var binIndex = _comparer.GetBinNumber(prefixFragmentMass);

                if (binIndex >= 0 && binIndex < numNodes)
                    prefixFragScores[binIndex] = GetMatchedIonPeakScore(true, peak);

                var suffixIonMass = peak.Mass;
                var suffixFragmentMass = suffixIonMass - SuffixOffsetMass;
                prefixFragmentMass = proteinMass - suffixFragmentMass;
                binIndex = _comparer.GetBinNumber(prefixFragmentMass);

                if (binIndex >= 0 && binIndex < numNodes)
                {
                    suffixFragScores[binIndex] = GetMatchedIonPeakScore(false, peak);

                    if (prefixFragScores[binIndex] != null) // there are complementary pair ions
                        suffixFragScores[binIndex] += ScoreParam.ComplementaryIonCount;
                }
            }

            return nodeScores;
        }

        private double GetMatchedIonPeakScore(bool isPrefixIon, DeconvolutedPeak peak)
        {
            var param = (isPrefixIon) ? ScoreParam.Prefix : ScoreParam.Suffix;
            var score = param.Count;
            score += param.Intensity * Math.Min(peak.Intensity / ReferencePeakIntensity, 1.0d);
            score += param.Corr * peak.Corr;
            score += param.Dist * peak.Dist;
            return score;
        }

        private readonly Dictionary<int, DeconvolutedPeak> _massBinToPeakMap;
        private readonly IMassBinning _comparer;
    }
}
