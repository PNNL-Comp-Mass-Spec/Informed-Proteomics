using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables
{
    public class MassErrorTable
    {
        public MassErrorTable(IonType[] ionTypes, Tolerance tolerance, double width=0.2, double binWidth=0.01, double offset=0.0)
        {
            _ionTypes = ionTypes;
            _totalPairs = 0;
            _tolerance = tolerance;
            _width = width;
            _binWidth = binWidth;
            _offset = offset;
            _massError = new Dictionary<IonType, Histogram<double>>();
            _ionPairFrequency = new Dictionary<IonType, Histogram<IonPairFound>>();
            var binEdges = GenerateEdges();
            foreach (var ionType in _ionTypes)
            {
                _massError.Add(ionType, new Histogram<double>(binEdges));
                _ionPairFrequency.Add(ionType, new Histogram<IonPairFound>((IonPairFound[])Enum.GetValues(typeof(IonPairFound))));
            }
        }

        public void AddMatch(SpectrumMatch match)
        {
            foreach (var ionType in _ionTypes)
            {
                var charge = ionType.Charge;
                var sequence = match.Sequence;
                var pepSeq = ionType.IsPrefixIon
                    ? sequence.GetRange(0, sequence.Count - 1)
                    : sequence.GetRange(1, sequence.Count - 1);
                var ions = match.GetCleavageIons(ionType);

                var nextIonIndex = 1;
                while (nextIonIndex < ions.Count)
                {
                    // look for peaks for current ion and next ion
                    _totalPairs++;
                    var currentIonIndex = nextIonIndex - 1;
                    var currentMz = ions[currentIonIndex].GetMonoIsotopicMz();
                    var currentPeak = match.Spectrum.FindPeak(currentMz, _tolerance);
                    var nextMz = ions[nextIonIndex].GetMonoIsotopicMz();
                    var nextPeak = match.Spectrum.FindPeak(nextMz, _tolerance);

                    if (currentPeak == null && nextPeak == null)
                        _ionPairFrequency[ionType].AddDatum(IonPairFound.Neither);
                    else if (nextPeak == null)
                        _ionPairFrequency[ionType].AddDatum(IonPairFound.First);
                    else if (currentPeak == null)
                        _ionPairFrequency[ionType].AddDatum(IonPairFound.Second);
                    else
                    {
                        // found both peaks, compute mass error
                        _ionPairFrequency[ionType].AddDatum(IonPairFound.Both);
                        var aaIndex = (ionType.IsPrefixIon ? nextIonIndex : currentIonIndex);
                        var aaMz = pepSeq[aaIndex].Mass/charge;
                        var massError = Math.Abs(nextPeak.Mz - currentPeak.Mz) - aaMz;
                        _massError[ionType].AddDatum(massError);
                    }
                    nextIonIndex++;
                }
            }
        }

        public void AddMatches(List<SpectrumMatch> matchList)
        {
            foreach (var match in matchList)
            {
                AddMatch(match);
            }
        }

        public void WriteToFile(StreamWriter file, IonType[] selectedIonTypes)
        {
        }

        public Probability<double> GetMassErrorProbability(double massError, IonType ionType)
        {
            var index = _massError[ionType].GetBinIndex(massError);
            if (index < 0) return new Probability<double>(massError, 0, _massError[ionType].Total);
            var found = _massError[ionType].Bins[index].Count;
            var total = _massError[ionType].Total;
            return new Probability<double>(massError, found, total);
        }

        public Probability<double>[] GetProbabilities(IonType ionType)
        {
            var bins = _massError[ionType].Bins;
            var binEdges = _massError[ionType].GetAlignedBinEdges(BinEdgeAlignment.Center, _binWidth);
            return binEdges.Select((t, i) => new Probability<double>(t, bins[i].Count, _totalPairs)).ToArray();
        }

        public double[] GetBinEdges(IonType ionType)
        {
            return _massError[ionType].BinEdges;
        }

        private double[] GenerateEdges()
        {
            var binEdges = new List<double>();
            for (double width = -1 * _offset; width >= -1 * _width; width -= _binWidth)
            {
                binEdges.Add(width);
            }
            for (double width = _offset; width < _width; width += _binWidth)
            {
                binEdges.Add(width);
            }

            binEdges = binEdges.Distinct().ToList();
            binEdges.Sort();
            return binEdges.ToArray();
        }

        private readonly IonType[] _ionTypes;
        private readonly Tolerance _tolerance;
        private readonly double _width;
        private readonly double _binWidth;
        private readonly double _offset;

        private readonly Dictionary<IonType, Histogram<double>> _massError;
        private readonly Dictionary<IonType, Histogram<IonPairFound>> _ionPairFrequency;
        private int _totalPairs;
    }

    public enum IonPairFound
    {
        Neither = 0,
        Second = 1,
        First = 2,
        Both = 3
    };
}
