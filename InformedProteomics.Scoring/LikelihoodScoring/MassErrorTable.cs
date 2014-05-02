using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public enum IonPairFound
    {
        Neither = 0,
        Second = 1,
        First = 2,
        Both = 3
    };

    public class MassErrorTable: I1DProbabilityTable<double>
    {
        private readonly IonType[] _ionTypes;
        private readonly Tolerance _tolerance;
        private readonly double _width;
        private readonly double _binWidth;
        private readonly double _offset;

        private readonly Histogram<double> _massError;
        private readonly Histogram<IonPairFound> _ionPairFrequency;
        private int _totalPairs;

        public List<Probability<IonPairFound>> IonPairFrequency
        {
            get
            {
                var bins = _ionPairFrequency.Bins;
                var binEdges = _ionPairFrequency.BinEdges;
                return binEdges.Select((t, i) => new Probability<IonPairFound>(t, bins[i].Count, _ionPairFrequency.Total)).ToList();
            }
        }

        public MassErrorTable(IonType[] ionTypes, Tolerance tolerance, double width=0.2, double binWidth=0.01, double offset=0.0)
        {
            _ionTypes = ionTypes;
            _totalPairs = 0;
            _tolerance = tolerance;
            _width = width;
            _binWidth = binWidth;
            _offset = offset;
            _massError = new Histogram<double>();
            _ionPairFrequency = new Histogram<IonPairFound>((IonPairFound[])Enum.GetValues(typeof(IonPairFound)));
            GenerateEdges();
        }

        public MassErrorTable(IonType[] ionTypes, Tolerance tolerance, double[] bins, List<double> data)
        {
            _ionTypes = ionTypes;
            _totalPairs = 0;
            _width = bins[bins.Length - 1] - bins[0];
            _binWidth = bins[1] - bins[0];
            double offset = 0.0;
            _massError = new Histogram<double>(bins);
            for (int i = 0; i < data.Count; i++)
            {
                var dataList = new List<double>();
                for (int j = 0; j < data[i]; j++)
                {
                    dataList.Add(bins[i]);
                }
                _massError.AddData(dataList);
            }
            _ionPairFrequency = new Histogram<IonPairFound>((IonPairFound[]) Enum.GetValues(typeof (IonPairFound)));
        }

        private void GenerateEdges()
        {
            var binEdges = new List<double>();
            for (double width = -1*_offset; width >= -1 * _width; width -= _binWidth)
            {
                binEdges.Add(width);
            }
            for (double width = _offset; width < _width; width += _binWidth)
            {
                binEdges.Add(width);
            }

            binEdges = binEdges.Distinct().ToList();
            binEdges.Sort();
            _massError.BinEdges = binEdges.ToArray();
        }

        public void AddMatches(List<SpectrumMatch> matchList)
        {
            foreach (var match in matchList)
            {
                foreach (var ionType in _ionTypes)
                {
                    var charge = ionType.Charge;
                    var sequence = match.Sequence;
                    var pepSeq = ionType.IsPrefixIon ? sequence.GetRange(0, sequence.Count - 1) : 
                                                        sequence.GetRange(1, sequence.Count - 1);
                    var ions = match.GetCleavageIons(ionType);

                    var nextIonIndex = 1;
                    while(nextIonIndex < ions.Count)
                    {
                        // look for peaks for current ion and next ion
                        _totalPairs++;
                        var currIonIndex = nextIonIndex - 1;
                        var currMz = ions[currIonIndex].GetMonoIsotopicMz();
                        var currPeak = match.Spectrum.FindPeak(currMz, _tolerance);
                        var nextMz = ions[nextIonIndex].GetMonoIsotopicMz();
                        var nextPeak = match.Spectrum.FindPeak(nextMz, _tolerance);

                        if (currPeak == null && nextPeak == null)
                            _ionPairFrequency.AddDatum(IonPairFound.Neither);
                        else if (nextPeak == null)
                            _ionPairFrequency.AddDatum(IonPairFound.First);
                        else if (currPeak == null)
                            _ionPairFrequency.AddDatum(IonPairFound.Second);
                        else
                        {
                            // found both peaks, compute mass error
                            _ionPairFrequency.AddDatum(IonPairFound.Both);
                            var aaIndex = (ionType.IsPrefixIon ? nextIonIndex : currIonIndex);
                            var aaMz = pepSeq[aaIndex].GetMass() / charge;
                            var massError = Math.Abs(nextPeak.Mz - currPeak.Mz) - aaMz;
                            _massError.AddDatum(massError);
                        }
                        nextIonIndex++;
                    }
                }
            }
        }

        public Probability<double> GetMassErrorProbability(double massError)
        {
            var index = _massError.GetBinIndex(massError);
            var found = _massError.Bins[index].Count;
            var total = _massError.Total;
            return new Probability<double>(massError, found, total);
        }

        public Probability<double>[] GetProbabilities()
        {
            var bins = _massError.Bins;
            var binEdges = _massError.GetAlignedBinEdges(BinEdgeAlignment.Center, _binWidth);
            return binEdges.Select((t, i) => new Probability<double>(t, bins[i].Count, _totalPairs)).ToArray();
        }

        public double[] GetBinEdges()
        {
            return _massError.BinEdges;
        }
    }
}
