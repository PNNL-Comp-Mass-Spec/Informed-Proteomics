using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public enum IonPairFound
    {
        Neither = 0,
        SecondIon = 1,
        FirstIon = 2,
        Both = 3
    };

    public class MassErrorTable
    {
        private readonly IonType[] _ionTypes;
        private readonly Tolerance _tolerance;

        public Histogram<double> MassError { get; private set; }
        public Histogram<IonPairFound> IonPairFrequency { get; private set; } 

        public MassErrorTable(IonType[] ionTypes, Tolerance tolerance, double binWidth=0.01)
        {
            _ionTypes = ionTypes;
            _tolerance = tolerance;
            MassError = new Histogram<double>();
            IonPairFrequency = new Histogram<IonPairFound>((IonPairFound[])Enum.GetValues(typeof(IonPairFound)));
            GenerateEdges(binWidth);
        }

        private void GenerateEdges(double binWidth)
        {
            const double searchWidth = 100;
            var binEdges = new List<double>();
            for (double width = 0; width < searchWidth; width += binWidth)
            {
                binEdges.Add(width);
            }
            MassError.BinEdges = binEdges.ToArray();
        }

        public void AddMatches(List<SpectrumMatch> matchList)
        {
            var aminoAcidSet = new AminoAcidSet();
            foreach (var match in matchList)
            {
                foreach (var ionType in _ionTypes)
                {
                    var charge = ionType.Charge;
                    var peptide = ionType.IsPrefixIon ? match.GetPeptidePrefix() : match.GetPeptideSuffix();
                    var ions = match.GetCleavageIons(ionType);
                    int nextIonIndex = 1;
                    while(nextIonIndex < ions.Count)
                    {
                        // look for peak for current ion and next ion
                        var currIonIndex = nextIonIndex - 1;
                        var currMz = ions[currIonIndex].GetMonoIsotopicMz();
                        var nextMz = ions[nextIonIndex].GetMonoIsotopicMz();
                        var currPeak = match.Spectrum.FindPeak(currMz, _tolerance);
                        var nextPeak = match.Spectrum.FindPeak(nextMz, _tolerance);

                        if (currPeak == null && nextPeak == null)
                            IonPairFrequency.AddDatum(IonPairFound.Neither);
                        else if (nextPeak == null)
                            IonPairFrequency.AddDatum(IonPairFound.FirstIon);
                        else if (currPeak == null)
                            IonPairFrequency.AddDatum(IonPairFound.SecondIon);
                        else
                        {
                            // found both peaks, compute mass error
                            IonPairFrequency.AddDatum(IonPairFound.Both);
                            var aminoAcid = aminoAcidSet.GetAminoAcid(peptide[currIonIndex]);
                            var aaMass = aminoAcid.GetMass() / charge;
                            var massError = Math.Abs(Math.Abs(currMz - nextMz) - aaMass);
                            MassError.AddDatum(massError);
                        }
                        nextIonIndex++;
                    }
                }
            }
        }
    }
}
