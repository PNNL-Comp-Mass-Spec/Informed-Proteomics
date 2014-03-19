using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.Scoring
{
    public class LikelihoodScoringModel
    {
        public LikelihoodScoringModel(string scoringParamPath)
        {
            _scoringParamPath = scoringParamPath;
            ParseScoringParametersOldFormat();
        }
        public double GetScore(BaseIonType baseIonType, double score, double intensity)
        {
            var baseIonTypeIndex = baseIonType.IsPrefix ? 0 : 1;
            if (intensity <= 0) return _missingIonScore[baseIonTypeIndex];

            var intensityIndex = Array.BinarySearch(_intensityBoundaries, intensity);
            if (intensityIndex < 0) intensityIndex = ~intensityIndex;

            var scoreIndex = (int)(score*NumScoreBins);
            if (scoreIndex >= NumScoreBins) scoreIndex = NumScoreBins - 1;

            return _ionScore[intensityIndex, scoreIndex, baseIonTypeIndex];
        }

        private readonly string _scoringParamPath;
        private double[] _intensityBoundaries;
        private double[] _missingIonScore;  // 0: b, 1: y
        private double[,,] _ionScore;   // IntensityBin, ScoreBin, IonBin
        private const int NumScoreBins = 20;
        private const int NumBaseIons = 2;

        private void ParseScoringParameters()
        {
            var lines = File.ReadAllLines(_scoringParamPath);

            var intensityBoundList = new List<double>();

            foreach (var line in lines)
            {
                if (!line.StartsWith("Intensity")) continue;
                var token = line.Split('\t');
                if (token.Length != 2) throw new Exception("Illegal scoring parameter file!");
                var intensity = Convert.ToDouble(token[1]);
                if (Math.Abs(intensity) <= 0) continue; // ignore zero intensity == 0
                intensityBoundList.Add(intensity);
            }
            _intensityBoundaries = intensityBoundList.ToArray();

            _ionScore = new double[_intensityBoundaries.Length + 1, NumScoreBins, NumBaseIons];
            _missingIonScore = new double[NumBaseIons];

            var intensityIndex = -1;
            var lineNum = 0;
            while (lineNum < lines.Length)
            {
                var line = lines[lineNum++];
                if (line.StartsWith("Range"))
                {
                    ++intensityIndex;
                    for (var scoreBin = 0; scoreBin < NumScoreBins; scoreBin++)  // this needs to be fixed
                    {
                        var data = lines[lineNum++];
                        var token = data.Split('\t');
                        if (token.Length != NumBaseIons * 2 + 1) throw new Exception("Illegal scoring parameter file!");

                        var numbers = token.Select(Convert.ToDouble).ToArray();

                        for (int i = 1; i < 5; i++)
                        {
                            if (numbers[i] <= 0) numbers[i] = 1;
                        }

                        _ionScore[intensityIndex, scoreBin, 0] =
                            Math.Log10(numbers[1] / numbers[2]);
                        _ionScore[intensityIndex, scoreBin, 1] =
                            Math.Log10(numbers[3] / numbers[4]);
                    }
                }
                else if (line.StartsWith("b"))
                {
                    var token = line.Split('\t');
                    if (token.Length != 3) throw new Exception("Illegal scoring parameter file!");
                    _missingIonScore[0] = Math.Log10(Convert.ToDouble(token[1]) / Convert.ToDouble(token[2]));
                }
                else if (line.StartsWith("y"))
                {
                    var token = line.Split('\t');
                    if (token.Length != 3) throw new Exception("Illegal scoring parameter file!");
                    _missingIonScore[1] = Math.Log10(Convert.ToDouble(token[1]) / Convert.ToDouble(token[2]));
                }
            }
        }

        private void ParseScoringParametersOldFormat()
        {
            var lines = File.ReadAllLines(_scoringParamPath);

            var intensityBoundList = new List<double>();

            foreach (var line in lines)
            {
                if (!line.StartsWith("Intensity")) continue;
                var token = line.Split('\t');
                if (token.Length != 2) throw new Exception("Illegal scoring parameter file!");
                var intensity = Convert.ToDouble(token[1]);
                if (Math.Abs(intensity) <= 0) continue; // ignore zero intensity == 0
                intensityBoundList.Add(intensity);
            }
            _intensityBoundaries = intensityBoundList.ToArray();

            _ionScore = new double[_intensityBoundaries.Length+1, NumScoreBins, NumBaseIons];
            _missingIonScore = new double[NumBaseIons];

            var intensityIndex = -1;
            var lineNum = 0;
            while(lineNum < lines.Length)
            {
                var line = lines[lineNum++];
                if (line.StartsWith("Range"))
                {
                    ++intensityIndex;
                    for (var scoreBin = 0; scoreBin < NumScoreBins; scoreBin++)  // this needs to be fixed
                    {
                        var data = lines[lineNum++];
                        var token = data.Split('\t');
                        if (token.Length != NumBaseIons*2 + 1) throw new Exception("Illegal scoring parameter file!");

                        var numbers = token.Select(Convert.ToDouble).ToArray();

                        for (int i = 1; i < 5; i++)
                        {
                            if (numbers[i] <= 0) numbers[i] = 1;
                        }

                        _ionScore[intensityIndex, scoreBin, 0] =
                            Math.Log10(numbers[1] / numbers[2]);
                        _ionScore[intensityIndex, scoreBin, 1] =
                            Math.Log10(numbers[3] / numbers[4]);
                    }
                }
                else if (line.StartsWith("b"))
                {
                    var token = line.Split('\t');
                    if (token.Length != 3) throw new Exception("Illegal scoring parameter file!");
                    _missingIonScore[0] = Math.Log10(Convert.ToDouble(token[1])/Convert.ToDouble(token[2]));
                }
                else if (line.StartsWith("y"))
                {
                    var token = line.Split('\t');
                    if (token.Length != 3) throw new Exception("Illegal scoring parameter file!");
                    _missingIonScore[1] = Math.Log10(Convert.ToDouble(token[1]) / Convert.ToDouble(token[2]));
                }
            }
        }
    }
}
