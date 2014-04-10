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

    public class MassErrorTable
    {
        private readonly IonType[] _ionTypes;
        private readonly Tolerance _tolerance;

        private readonly Histogram<double> _massError;
        private int _totalPairs;
        private readonly Histogram<IonPairFound> _ionPairFrequency;

        public List<Probability<double>> MassError
        {
            get
            {
                var bins = _massError.Bins;
                var binEdges = _massError.BinEdges;
                return binEdges.Select((t, i) => new Probability<double>(t, bins[i].Count, _totalPairs)).ToList();
            }
        }

        public List<Probability<IonPairFound>> IonPairFrequency
        {
            get
            {
                var bins = _ionPairFrequency.Bins;
                var binEdges = _ionPairFrequency.BinEdges;
                return binEdges.Select((t, i) => new Probability<IonPairFound>(t, bins[i].Count, _ionPairFrequency.Total)).ToList();
            }
        }

        public MassErrorTable(IonType[] ionTypes, Tolerance tolerance, double searchWidth=0.2, double binWidth=0.01)
        {
            _ionTypes = ionTypes;
            _totalPairs = 0;
            _tolerance = tolerance;
            _massError = new Histogram<double>();
            _ionPairFrequency = new Histogram<IonPairFound>((IonPairFound[])Enum.GetValues(typeof(IonPairFound)));
            GenerateEdges(searchWidth, binWidth);
        }

        private void GenerateEdges(double searchWidth, double binWidth)
        {
            var binEdges = new List<double>();
            for (double width = 0; width >= -1 * searchWidth; width -= binWidth)
            {
                binEdges.Add(width);
            }
            for (double width = 0; width < searchWidth; width += binWidth)
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
    }
}
