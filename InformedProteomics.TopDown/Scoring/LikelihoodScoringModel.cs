using System;
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
            ParseScoringParameters();
        }

        public double GetScore(BaseIonType baseIonType, double score, double intensity)
        {
            var ionTypeIndex = baseIonType.IsPrefix ? 0 : 1;
            if (intensity <= 0) return _missingIonScore[ionTypeIndex];

            var intensityIndex = Array.BinarySearch(_intensityBoundaries[ionTypeIndex], intensity);
            if (intensityIndex < 0) intensityIndex = ~intensityIndex;
            if(intensityIndex > 0) intensityIndex -= 1;

            var scoreIndex = (int)(score*_numScoreBins);
            if (scoreIndex >= _numScoreBins) scoreIndex = _numScoreBins - 1;

            //if (intensityIndex < 0 || intensityIndex >= _numIntensityBins)
            //{
            //    Console.Write("IntensityIndex: {0}", intensityIndex);
            //}
            //if (scoreIndex < 0 || scoreIndex >= _numScoreBins)
            //{
            //    Console.Write("ScoreIndex: {0}", scoreIndex);
            //}
            return _ionScore[ionTypeIndex, intensityIndex, scoreIndex];
        }

        private readonly string _scoringParamPath;
        private int _numIons;
        private string[] _ions;
        private int _numScoreBins;
        private int _numIntensityBins;

        private double[] _missingIonScore;  // ion
        private double[][] _intensityBoundaries; // ion, intensity
        private double[][] _scoreBoundaries; // ion, score
        private double[, ,] _ionScore;   // ion, intensity, score

        private void ParseScoringParameters()
        {
            var lines = File.ReadAllLines(_scoringParamPath);

            if (lines.Length == 0) throw new Exception("Illegal scoring parameter file!");

            var lineNum = 0;

            // First line: Ions
            if (!lines[lineNum++].StartsWith("Ions")) throw new Exception("Illegal scoring parameter file!");
            var ionsToken = lines[0].Split('\t');
            _numIons = Convert.ToInt32(ionsToken[1]);
            _ions = new string[_numIons];
            for (var i = 0; i < _numIons; i++) _ions[i] = ionsToken[2 + i];
            _missingIonScore = new double[_numIons];

            // Second line: NumScoreBins
            if (!lines[lineNum++].StartsWith("NumScoreBins")) throw new Exception("Illegal scoring parameter file!");
            var numScoreBinsToken = lines[1].Split('\t');
            _numScoreBins = Convert.ToInt32(numScoreBinsToken[1]);
            _scoreBoundaries = new double[_numIons][];

            // Third line: NumIntensityBins
            if (!lines[lineNum++].StartsWith("NumIntensityBins")) throw new Exception("Illegal scoring parameter file!");
            var numIntensityBinsToken = lines[2].Split('\t');
            _numIntensityBins = Convert.ToInt32(numIntensityBinsToken[1]);
            _intensityBoundaries = new double[_numIons][];

            _ionScore = new double[_numIons, _numIntensityBins, _numScoreBins];

            for (var ionIndex = 0; ionIndex < _numIons; ionIndex++)
            {
                var line = lines[lineNum++];
                var ionToken = line.Split('\t');
                if (ionToken.Length != 2 || !ionToken[1].Equals(_ions[ionIndex])) throw new Exception("Illegal scoring parameter file!");

                line = lines[lineNum++];
                var intensityToken = line.Split('\t');
                if (intensityToken.Length != _numIntensityBins + 1) throw new Exception("Illegal scoring parameter file!");
                _intensityBoundaries[ionIndex] = new double[_numIntensityBins];
                for (var i = 0; i < _numIntensityBins; i++) _intensityBoundaries[ionIndex][i] = Convert.ToDouble(intensityToken[i + 1]);

                line = lines[lineNum++];
                var scoreToken = line.Split('\t');
                if (scoreToken.Length != _numScoreBins + 1) throw new Exception("Illegal scoring parameter file!");
                _scoreBoundaries[ionIndex] = new double[_numScoreBins];
                for (var i = 0; i < _numScoreBins; i++) _scoreBoundaries[ionIndex][i] = Convert.ToDouble(scoreToken[i + 1]);

                line = lines[lineNum++];
                var missingIonToken = line.Split('\t');
                if (missingIonToken.Length != 3) throw new Exception("Illegal scoring parameter file!");
                var missingIonTarget = Convert.ToDouble(missingIonToken[1]);
                var missingIonDecoy = Convert.ToDouble(missingIonToken[2]);
                _missingIonScore[ionIndex] = Math.Log10(missingIonTarget/missingIonDecoy);

                for (var intensityIndex = 0; intensityIndex < _numIntensityBins; intensityIndex++)
                {
                    // Target
                    line = lines[lineNum++];
                    var targetToken = line.Split('\t');
                    if (targetToken.Length != _numScoreBins) throw new Exception("Illegal scoring parameter file!");
                    var targetNumbers = targetToken.Select(t => Convert.ToInt32(t)).ToArray();

                    // Decoy
                    line = lines[lineNum++];
                    var decoyToken = line.Split('\t');
                    if (decoyToken.Length != _numScoreBins) throw new Exception("Illegal scoring parameter file!");
                    var decoyNumbers = decoyToken.Select(t => Convert.ToInt32(t)).ToArray();

                    for (var scoreIndex = 0; scoreIndex < _numScoreBins; scoreIndex++)
                    {
                        double score;
                        if (targetNumbers[scoreIndex] <= 0 || decoyNumbers[scoreIndex] <= 0)
                        {
                            score = _missingIonScore[ionIndex];
                        }
                        else
                        {
                            score = Math.Log10(targetNumbers[scoreIndex]/(double) decoyNumbers[scoreIndex]);
                        }
                        _ionScore[ionIndex, intensityIndex, scoreIndex] = score;
                    }
                }
            }

            if (lineNum != lines.Length) throw new Exception("Illegal scoring parameter file!");
        }

        //private void ParseScoringParametersOldFormat()
        //{
        //    var lines = File.ReadAllLines(_scoringParamPath);

        //    var intensityBoundList = new List<double>();

        //    foreach (var line in lines)
        //    {
        //        if (!line.StartsWith("Intensity")) continue;
        //        var token = line.Split('\t');
        //        if (token.Length != 2) throw new Exception("Illegal scoring parameter file!");
        //        var intensity = Convert.ToDouble(token[1]);
        //        if (Math.Abs(intensity) <= 0) continue; // ignore zero intensity == 0
        //        intensityBoundList.Add(intensity);
        //    }
        //    _intensityBoundaries = intensityBoundList.ToArray();

        //    _ionScore = new double[_intensityBoundaries.Length+1, _numScoreBins, _numIons];
        //    _missingIonScore = new double[_numIons];

        //    var intensityIndex = -1;
        //    var lineNum = 0;
        //    while(lineNum < lines.Length)
        //    {
        //        var line = lines[lineNum++];
        //        if (line.StartsWith("Range"))
        //        {
        //            ++intensityIndex;
        //            for (var scoreBin = 0; scoreBin < _numScoreBins; scoreBin++)  // this needs to be fixed
        //            {
        //                var data = lines[lineNum++];
        //                var token = data.Split('\t');
        //                if (token.Length != _numIons*2 + 1) throw new Exception("Illegal scoring parameter file!");

        //                var numbers = token.Select(Convert.ToDouble).ToArray();

        //                for (int i = 1; i < 5; i++)
        //                {
        //                    if (numbers[i] <= 0) numbers[i] = 1;
        //                }

        //                _ionScore[intensityIndex, scoreBin, 0] =
        //                    Math.Log10(numbers[1] / numbers[2]);
        //                _ionScore[intensityIndex, scoreBin, 1] =
        //                    Math.Log10(numbers[3] / numbers[4]);
        //            }
        //        }
        //        else if (line.StartsWith("b"))
        //        {
        //            var token = line.Split('\t');
        //            if (token.Length != 3) throw new Exception("Illegal scoring parameter file!");
        //            _missingIonScore[0] = Math.Log10(Convert.ToDouble(token[1])/Convert.ToDouble(token[2]));
        //        }
        //        else if (line.StartsWith("y"))
        //        {
        //            var token = line.Split('\t');
        //            if (token.Length != 3) throw new Exception("Illegal scoring parameter file!");
        //            _missingIonScore[1] = Math.Log10(Convert.ToDouble(token[1]) / Convert.ToDouble(token[2]));
        //        }
        //    }
        //}
    }
}
