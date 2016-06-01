using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables
{
    public enum ScoreMethod { Cosine, FitScore, Pearson }
    public class CorrelationScoreTable
    {
        public double[] IntensityBins { get; private set; }
        public Probability<int> WorstScore { get; private set; }

        public CorrelationScoreTable(ScoreMethod method, int intensityBins, double[] binEdges)
        {
            _method = method;
            _binEdges = binEdges;
            IntensityBins = null;
            _intensityBinCount = intensityBins;
            WorstScore = new Probability<int>(0);

            _intensityHistogram = new Histogram<FitScore>(new CompareFitScoreByIntensity());
        }

        public CorrelationScoreTable(ScoreMethod method, double[] intensityBins, double[] binEdges)
        {
            _method = method;
            _binEdges = binEdges;
            IntensityBins = intensityBins;
            _intensityBinCount = intensityBins.Length;
            WorstScore = new Probability<int>(0);

            _intensityHistogram = new Histogram<FitScore>(new CompareFitScoreByIntensity());

        }

        public List<Histogram<FitScore>> Histograms
        {
            get
            {
                if (IntensityBins == null)
                {
                    _intensityHistogram.Equalize(_intensityBinCount, new FitScore(0, 0));
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

        public void AddMatches(List<SpectrumMatch> matches, IonType[] ionTypes, Tolerance tolerance, double relativeIntensityThreshold, bool reduceCharges=true)
        {
            foreach (var match in matches)
            {
                var spectrum = match.Spectrum;

                var prefixes = match.Prefixes;
                var suffixes = match.Suffixes;

                StreamWriter debugFile = File.AppendText(@"C:\Users\wilk011\Documents\DataFiles\TestFolder\bCID.txt");

                for (int i = 0; i < prefixes.Count; i++)
                {
                    var ionTypeScores = new Dictionary<string, FitScoreList>();
                    foreach (var ionType in ionTypes)
                    {
                        var cleavagePoints = ionType.BaseIonType.IsPrefix ? prefixes : suffixes;

                        var ion = ionType.GetIon(cleavagePoints[i]);

                        var mostIntensePeak = GetHighestPeak(ion, spectrum, tolerance, relativeIntensityThreshold);
                        double score = GetScore(ion, spectrum, tolerance, relativeIntensityThreshold);

                        double intensity = -1.0;

                        if (mostIntensePeak != null)
                            intensity = mostIntensePeak.Intensity;

                        var name = ionType.Name;
                        if (reduceCharges)
                            name = ReducedChargeName(ionType);

                        if (!ionTypeScores.ContainsKey(name))
                            ionTypeScores.Add(name, new FitScoreList());

                        ionTypeScores[name].Add(new FitScore(intensity, score));
                    }
                    var bestScores = SelectBestScores(ionTypeScores);
                    foreach (var bestscore in bestScores)
                    {
                        var score = bestscore.Score;
                        if (_method == ScoreMethod.FitScore)
                            score = 1 - score;
                        WorstScore.Total++;
                        if (score.Equals(0))
                            WorstScore.Found++;

                        debugFile.WriteLine("{0}\t{1}", bestscore.Intensity, score);
                    }
                    _intensityHistogram.AddData(bestScores);
                }
                debugFile.Close();
            }
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

        private static string ReducedChargeName(IonType ionType)
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

        private readonly double[] _binEdges;
        private readonly ScoreMethod _method;
        private readonly Histogram<FitScore> _intensityHistogram;
        private readonly int _intensityBinCount;
    }
}
