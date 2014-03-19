using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public enum ScoreMethod { Cosine, FitScore, Pearson }
    public class ScoreTable
    {
        private readonly double[] _binEdges;
        private readonly ScoreMethod _method;
        private readonly Histogram<FitScore> _intensityHistogram;
        private readonly int _intensityBinCount;
        public double[] IntensityBins { get; private set; }
        public Probability WorstScore { get; private set; }

        public List<Histogram<FitScore>> Histograms
        {
            get
            {
                if (IntensityBins == null)
                {
                    _intensityHistogram.Equalize(_intensityBinCount);
                    var edgeList = new FitScoreList(_intensityHistogram.BinEdges);
                    IntensityBins = edgeList.Intensities;
                }
                else
                {
                    _intensityHistogram.BinEdges = new FitScoreList(IntensityBins, null).ToArray();
                }
                var bins = _intensityHistogram.Bins;
                return bins.Select(bin => new Histogram<FitScore>(bin, (new FitScoreList(null, _binEdges)).ToArray(), new CompareFitScoreByScore())).ToList();
            }
        }

        public ScoreTable(ScoreMethod method, int intensityBins, double[] binEdges)
        {
            _method = method;
            _binEdges = binEdges;
            IntensityBins = null;
            _intensityBinCount = intensityBins;
            WorstScore = new Probability();

            _intensityHistogram = new Histogram<FitScore>(new CompareFitScoreByIntensity());
        }

        public ScoreTable(ScoreMethod method, double[] intensityBins, double[] binEdges)
        {
            _method = method;
            _binEdges = binEdges;
            IntensityBins = intensityBins;
            _intensityBinCount = intensityBins.Length;
            WorstScore = new Probability();

            _intensityHistogram = new Histogram<FitScore>(new CompareFitScoreByIntensity());

        }

        private string ReducedChargeName(IonType ionType)
        {
            string name = ionType.Name;
            if (ionType.Charge > 1 && ionType.Charge < 10)
                name = name.Remove(1, 1);
            if (ionType.Charge >= 10)
                name = name.Remove(1, 2);
            return name;
        }

        private List<FitScore> SelectBestScores(Dictionary<string, FitScoreList> ionTypeScores)
        {
            return ionTypeScores.Values.Select(scoreList => scoreList.MaxScore(_method)).ToList();
        }

        private double GetScore(Ion ion, Spectrum spectrum, Tolerance tolerance, double relativeIntenistyThreshold)
        {
            double score;
            switch (_method)
            {
                case ScoreMethod.Cosine:
                    score = spectrum.GetConsineScore(ion, tolerance, relativeIntenistyThreshold);
                    break;
                case ScoreMethod.FitScore:
                    score = spectrum.GetFitScore(ion, tolerance, relativeIntenistyThreshold);
                    break;
                case ScoreMethod.Pearson:
                    score = spectrum.GetCorrScore(ion, tolerance, relativeIntenistyThreshold);
                    break;
                default:
                    score = spectrum.GetConsineScore(ion, tolerance, relativeIntenistyThreshold);
                    break;
            }
            return score;
        }

        private Peak GetHighestPeak(Ion ion, Spectrum spectrum, Tolerance tolerance, double relativeIntensityThreshold)
        {
            var peaks = spectrum.GetAllIsotopePeaks(ion, tolerance, relativeIntensityThreshold);

            Peak highestPeak = null;
            double highestIntensity = 0.0;

            if (peaks == null) return null;

            foreach (var peak in peaks)
            {
                if (peak != null && peak.Intensity >= highestIntensity)
                {
                    highestPeak = peak;
                    highestIntensity = peak.Intensity;
                }
            }
            return highestPeak;
        }

        public void AddMatches(List<SpectrumMatch> matches, IonType[] ionTypes, Tolerance tolerance, double relativeIntensityThreshold, bool reduceCharges=true)
        {
            var scores = new List<FitScore>();

            foreach (var match in matches)
            {
                var spectrum = match.Spectrum;

                var prefixes = match.Prefixes;
                var suffixes = match.Suffixes;

                for (int i = 0; i < prefixes.Count; i++)
                {
                    var ionTypeScores = new Dictionary<string, FitScoreList>();
                    foreach (var ionType in ionTypes)
                    {
                        var cleavagePoints = new List<Composition>();
                        if (ionType.BaseIonType == BaseIonType.A || ionType.BaseIonType == BaseIonType.B || ionType.BaseIonType == BaseIonType.C)
                            cleavagePoints = prefixes;
                        else if (ionType.BaseIonType == BaseIonType.X || ionType.BaseIonType == BaseIonType.Y || ionType.BaseIonType == BaseIonType.Z)
                            cleavagePoints = suffixes;

                        var ion = ionType.GetIon(cleavagePoints[i]);

                        var mostIntensePeak = GetHighestPeak(ion, spectrum, tolerance, relativeIntensityThreshold);
                        double score = GetScore(ion, spectrum, tolerance, relativeIntensityThreshold);

                        double intensity = -1.0;

                        if (mostIntensePeak != null)
                            intensity = mostIntensePeak.Intensity;

                        var name = ionType.Name;
                        if (reduceCharges)
                            name = ReducedChargeName(ionType);

                        if (intensity >= 0)
                        {
                            if (!ionTypeScores.ContainsKey(name))
                                ionTypeScores.Add(name, new FitScoreList());

                            ionTypeScores[name].Add(new FitScore(intensity, score));
                        }
                        else
                            WorstScore.Found++;

                        WorstScore.Total++;
                    }
                    var bestScores = SelectBestScores(ionTypeScores);
                    scores.AddRange(bestScores);
                }
            }

            _intensityHistogram.AddData(scores);
        }
    }
}
