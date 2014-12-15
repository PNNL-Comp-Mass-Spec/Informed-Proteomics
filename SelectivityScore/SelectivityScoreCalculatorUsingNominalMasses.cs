using System;
using System.Collections.Generic;

namespace SelectivityScore
{
    public class SelectivityScoreCalculatorUsingNominalMasses
    {
        public SelectivityScoreCalculatorUsingNominalMasses(
            IList<int> integerAminoAcidMasses,
            IList<double> aminoAcidProbabilities,
            int minPrecursorCharge = 2,
            int maxPrecursorCharge = 3,
            int maxMass = 3000
            )
        {
            _integerAminoAcidMasses = integerAminoAcidMasses;
            _aminoAcidProbabilities = aminoAcidProbabilities;
            _minPrecursorCharge = minPrecursorCharge;
            _maxPrecursorCharge = maxPrecursorCharge;
            _maxMass = maxMass;
            PrecomputeProbabilities();
        }

        public double[] GetPValues(
            double precursorMzTarget,
            double[] productMzTargets
            )
        {
            var numTransitions = productMzTargets.Length;
            var combinations = CreateNtoTheKCombinations(3, numTransitions);

            var pValues = new double[numTransitions + 1];
            
            for (var precursorCharge = _minPrecursorCharge; precursorCharge <= _maxPrecursorCharge; precursorCharge++)
            {
                var intPeptideMass = GetIntMass((precursorMzTarget - Proton) * precursorCharge - Water);
                foreach (var combination in combinations)
                {
                    var suffixList = new List<int>();
                    var score = 0;
                    for (var transitionIndex = 0; transitionIndex < combination.Length; transitionIndex++)
                    {
                        // Considering ith transition
                        var type = combination[transitionIndex];
                        if (type == 1) // Use as b
                        {
                            var intSuffixMass = intPeptideMass - GetIntMass(productMzTargets[transitionIndex] - OffsetB);
                            if(intSuffixMass > 0 && intSuffixMass < intPeptideMass) suffixList.Add(intSuffixMass);
                            ++score;
                        }
                        else if (type == 2) // Use as y
                        {
                            var intSuffixMass = GetIntMass(productMzTargets[transitionIndex] - OffsetY);
                            if (intSuffixMass > 0 && intSuffixMass < intPeptideMass) suffixList.Add(intSuffixMass);
                            ++score;
                        }
                        // Not used if type == 0
                    }
                    suffixList.Sort();
                    suffixList.Add(intPeptideMass);
                    var pValue = 1.0;
                    for (var i = 1; i < suffixList.Count; i++)
                    {
                        pValue *= _probability[suffixList[i] - suffixList[i - 1]];
                    }
                    pValues[score] += pValue * GetPriorProbability(intPeptideMass, precursorCharge);
                }
            }
            return pValues;
        }

        public double[] GetPValues(
            double precursorMzTarget,
            double precursorMzWindow,
            double[] productMzTargets,
            double[] productMzWindows
            )
        {
            var numTransitions = productMzTargets.Length;
            var combinations = CreateNtoTheKCombinations(3, numTransitions);

            var pValues = new double[numTransitions + 1];

            for (var precursorCharge = _minPrecursorCharge; precursorCharge <= _maxPrecursorCharge; precursorCharge++)
            {
                var minIntPeptideMass = GetIntMass((precursorMzTarget - precursorMzWindow / 2 - Proton) * precursorCharge - Water);
                var maxIntPeptideMass = GetIntMass((precursorMzTarget + precursorMzWindow / 2 - Proton) * precursorCharge - Water);

                for (var intPeptideMass = minIntPeptideMass; intPeptideMass <= maxIntPeptideMass; intPeptideMass++)
                {
                    foreach (var combination in combinations)
                    {
                        var suffixList = new List<int>();
                        var score = 0;
                        for (var transitionIndex = 0; transitionIndex < combination.Length; transitionIndex++)
                        {
                            // Considering ith transition
                            var type = combination[transitionIndex];
                            if (type == 1) // Use as b
                            {
                                var intSuffixMass = intPeptideMass - GetIntMass(productMzTargets[transitionIndex] - OffsetB);
                                if (intSuffixMass > 0 && intSuffixMass < intPeptideMass) suffixList.Add(intSuffixMass);
                                ++score;
                            }
                            else if (type == 2) // Use as y
                            {
                                var intSuffixMass = GetIntMass(productMzTargets[transitionIndex] - OffsetY);
                                if (intSuffixMass > 0 && intSuffixMass < intPeptideMass) suffixList.Add(intSuffixMass);
                                ++score;
                            }
                            // Not used if type == 0
                        }
                        suffixList.Sort();
                        suffixList.Add(intPeptideMass);
                        var pValue = 1.0;
                        for (var i = 1; i < suffixList.Count; i++)
                        {
                            pValue *= _probability[suffixList[i] - suffixList[i - 1]];
                        }
                        pValues[score] += pValue * GetPriorProbability(intPeptideMass, precursorCharge);
                    }                    
                }
            }
            return pValues;
        }

        public const double Proton = 1.00727649;
        public const double Water = 18.0105647;
        public const double OffsetB = Proton;
        public const double OffsetY = Water + Proton;

        private readonly IList<int> _integerAminoAcidMasses;
        private readonly IList<double> _aminoAcidProbabilities;
        private const double RescalingConstant = 0.9995;

        private readonly int _maxMass;
        private readonly int _minPrecursorCharge;
        private readonly int _maxPrecursorCharge;

        private double[] _probability;

        // TODO: this must be computed from real data
        private double GetPriorProbability(int peptideMass, int charge)
        {
            if (charge == 2) return 0.6;
            return 0.4;
        }

        private void PrecomputeProbabilities()
        {
            _probability = new double[_maxMass + 1];
            _probability[0] = 1.0;

            for (var i = 1; i < _probability.Length; i++)
            {
                for (var j = 0; j < _integerAminoAcidMasses.Count; j++)
                {
                    var prevMass = i - _integerAminoAcidMasses[j];
                    if (prevMass < 0) continue;
                    _probability[i] += _probability[prevMass] * _aminoAcidProbabilities[j];
                }
            }
        }

        private int GetIntMass(double m)
        {
            return (int)Math.Round(m * RescalingConstant);
        }

        private static int[][] CreateNtoTheKCombinations(int n, int length)
        {
            if (n <= 0)
                return null;

            if (length == 0)
            {
                return new[] { new int[0] };
            }
            if (length == 1)
            {
                var combinations = new int[n][];
                for (var i = 0; i < n; i++)
                {
                    combinations[i] = new[] { i };
                }
                return combinations;
            }
            else
            {
                var prevCombinations = CreateNtoTheKCombinations(n, length - 1);
                var combinations = new List<int[]>();
                foreach (var combination in prevCombinations)
                {
                    for (var j = 0; j < n; j++)
                    {
                        var newCombination = new int[combination.Length + 1];
                        Array.Copy(combination, newCombination, combination.Length);
                        newCombination[newCombination.Length - 1] = j;
                        combinations.Add(newCombination);
                    }
                }
                return combinations.ToArray();
            }
        }
    }
}
