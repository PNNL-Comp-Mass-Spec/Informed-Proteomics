using System;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public class ProteinMassComparerWithBinning
    {
        public ProteinMassComparerWithBinning(AminoAcid[] aaSet, double maxProteinMass = 50000, int numBits = 27)
        {
            _maxProteinMass = maxProteinMass;
            _mzComparer = new MzComparerWithBinning(numBits);
            var maxBinIndex = _mzComparer.GetBinNumber(maxProteinMass);
            var nBins = maxBinIndex + 1;
            var binMap = new int[nBins];
            binMap[0] = 1;
            var binCounter = 1;

            for (var i = 0; i < nBins; i++)
            {
                var minMass = _mzComparer.GetMzStart(i);
                var maxMass = _mzComparer.GetMzEnd(i);
                if (binMap[i] == 1)
                {
                    foreach (var aa in aaSet)
                    {
                        for (var j = _mzComparer.GetBinNumber(minMass + aa.Mass);
                            j <= _mzComparer.GetBinNumber(maxMass + aa.Mass);
                            j++)
                        {
                            if (j > maxBinIndex) continue;
                            if (binMap[j] == 0)
                            {
                                binMap[j] = 1;
                                binCounter++;
                            }
                        }
                    }
                }
            }
            
            _binNumberToIndexMap = new int[binCounter];
            var idx = 0;
            for (var i = 0; i < nBins; i++)
            {
                if (binMap[i] > 0)
                {
                    binMap[i] = idx;
                    _binNumberToIndexMap[idx] = i;
                    idx++;
                }
            }
            _binMap = binMap;
        }

        public int GetBinNumber(double mass)
        {
            var originalBinIndex = _mzComparer.GetBinNumber(mass);

            if (originalBinIndex >= _binMap.Length) return -1;
            

            if (originalBinIndex == 0 || _binMap[originalBinIndex] > 0) return _binMap[originalBinIndex];

            var originalBinMass = _mzComparer.GetMzAverage(originalBinIndex);

            if (mass < originalBinMass && _binMap[originalBinIndex - 1] > 0)
            {
                return _binMap[originalBinIndex - 1];
            }
            if (mass >= originalBinMass && _binMap[originalBinIndex + 1] > 0)
            {
                return _binMap[originalBinIndex + 1];
            }

            return -1;
        }

        public double GetMass(int binNumber)
        {
            if(binNumber >= _binNumberToIndexMap.Length || binNumber < 0) 
                throw new IndexOutOfRangeException(string.Format("Out of bin range(max bin number={0})", _binNumberToIndexMap.Length-1));
            var originalBinIndex = _binNumberToIndexMap[binNumber];
            return _mzComparer.GetMzAverage(originalBinIndex);
        }
        
        private readonly int[] _binMap;
        private readonly int[] _binNumberToIndexMap;
        private readonly MzComparerWithBinning _mzComparer;
        private double _maxProteinMass;
    }
}
