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

            return _ionScore[ionTypeIndex, intensityIndex, scoreIndex];
        }

        public void PrintAllScores()
        {
            for (var ionBin = 0; ionBin < _numIons; ionBin++)
            {
                Console.WriteLine("{0} ion", ionBin == 0 ? "b" : "y");
                Console.WriteLine("\t" + string.Join("\t", _scoreBoundaries[ionBin]));
                for (var intBin = 0; intBin < _numIntensityBins; intBin++)
                {
                    Console.Write("{0:F2}-{1}",
                        _intensityBoundaries[ionBin][intBin], 
                        (intBin == _numIntensityBins-1 ? "" : string.Format("{0:F2}",_intensityBoundaries[ionBin][intBin + 1])));
                    for (var scoreBin = 0; scoreBin < _numScoreBins; scoreBin++)
                    {
                        Console.Write("\t"+_ionScore[ionBin, intBin, scoreBin]);
                    }
                    Console.WriteLine();
                }
            }
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

        private double[] _existingIonScore; // Dancik

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

            _existingIonScore = new double[_numIons];

        }
    }
}
