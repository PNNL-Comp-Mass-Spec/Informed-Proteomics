using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.TopDown
{
    public class ProteinMassBinning : IMassBinning
    {
        public ProteinMassBinning(double minMass = 50, double maxMass = 50000, bool padZeroMassBin = false, int numBits = 28)
        {
            _mzComparer = new MzComparerWithBinning(numBits);
            MinMass = minMass;
            MaxMass = maxMass;
            _minBinIndex = _mzComparer.GetBinNumber(MinMass);
            _maxBinIndex = _mzComparer.GetBinNumber(MaxMass);
            _padZeroMassBin = padZeroMassBin;
            NumberOfBins = (_padZeroMassBin) ? _maxBinIndex - _minBinIndex + 2 : _maxBinIndex - _minBinIndex + 1;
        }

        public bool Filtered => false;

        public bool CheckMassRange(double mass)
        {
            return (mass >= MinMass && mass <= MaxMass);
        }

        public int GetBinNumber(double mass)
        {
            var binNumber = _mzComparer.GetBinNumber(mass);
            if (_padZeroMassBin)
            {
                if (binNumber == 0)
                {
                    return 0;
                }

                if (binNumber < _minBinIndex)
                {
                    return 1;
                }

                if (binNumber >= _maxBinIndex)
                {
                    return NumberOfBins - 1;
                }

                return binNumber - _minBinIndex + 1;
            }

            if (binNumber < _minBinIndex)
            {
                return 0;
            }

            if (binNumber >= _maxBinIndex)
            {
                return NumberOfBins - 1;
            }

            return binNumber - _minBinIndex;
        }

        public double GetMass(int binNumber)
        {
            var oriBinNum = GetOriginalBinNumber(binNumber);
            return _mzComparer.GetMzAverage(oriBinNum);
        }

        public double GetMassStart(int binNumber)
        {
            var oriBinNum = GetOriginalBinNumber(binNumber);
            return _mzComparer.GetMzStart(oriBinNum);
        }

        public double GetMassEnd(int binNumber)
        {
            var oriBinNum = GetOriginalBinNumber(binNumber);
            return _mzComparer.GetMzEnd(oriBinNum);
        }

        private int GetOriginalBinNumber(int effectiveBinNumber)
        {
            if (_padZeroMassBin)
            {
                if (effectiveBinNumber == 0)
                {
                    return 0;
                }

                return effectiveBinNumber + _minBinIndex - 1;
            }
            return effectiveBinNumber + _minBinIndex;
        }

        private readonly MzComparerWithBinning _mzComparer;
        private readonly int _minBinIndex;
        private readonly int _maxBinIndex;
        private readonly bool _padZeroMassBin;

        public int NumberOfBins { get; }
        public double MaxMass { get; }
        public double MinMass { get; }
    }
}
