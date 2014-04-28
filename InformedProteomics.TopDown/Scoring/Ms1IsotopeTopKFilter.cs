using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Scoring
{
    public class Ms1IsotopeTopKFilter : ISequenceFilter
    {
        public Ms1IsotopeTopKFilter(
            LcMsRun run,
            int minCharge = 3, int maxCharge = 30,
            double ppmTolerance = 15)
        {
            _run = run;
            _minCharge = minCharge;
            _maxCharge = maxCharge;
            _ppmTolerance = ppmTolerance;
            _tolerance = new Tolerance(_ppmTolerance);
            _comparer = new MzComparerWithPpmTolerance(ppmTolerance);
            _sequenceMassBinToScanNumsMap = new Dictionary<int, IList<int>>();
            //PrecomputePossibleSequenceMasses();
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            var sequenceMassBinNum = ProductScorerBasedOnDeconvolutedSpectra.GetBinNumber(sequenceMass);
            IList<int> ms2ScanNums;
            if (_sequenceMassBinToScanNumsMap.TryGetValue(sequenceMassBinNum, out ms2ScanNums)) return ms2ScanNums;
            return new int[0];
        }

        //public int NumDeisotopedPeaks { get; private set; }

        private readonly LcMsRun _run;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly double _ppmTolerance;
        private readonly Tolerance _tolerance;
        private readonly MzComparerWithPpmTolerance _comparer;

        private readonly Dictionary<int, IList<int>> _sequenceMassBinToScanNumsMap; // mass -> scan numbers

        public void PrecomputePossibleSequenceMasses(int numPeaks)
        {
            var ms2ScanNums = _run.GetScanNumbers(2);
            foreach (var ms2ScanNum in ms2ScanNums)
            {
                var possibleSequenceMasses = GetPossibleSequenceMasses(ms2ScanNum, numPeaks);
                var binNumHash = new HashSet<int>();
                foreach (var sequenceMass in possibleSequenceMasses)
                {
                    var dMass = sequenceMass * _ppmTolerance * 1e-6;
                    var minMass = sequenceMass - dMass;
                    var maxMass = sequenceMass + dMass;
                    var minBinNum = ProductScorerBasedOnDeconvolutedSpectra.GetBinNumber(minMass);
                    var maxBinNum = ProductScorerBasedOnDeconvolutedSpectra.GetBinNumber(maxMass);
                    for (var binNum = minBinNum; binNum <= maxBinNum; binNum++)
                    {
                        if (!binNumHash.Add(binNum)) continue;
                        IList<int> scanNums;
                        if (_sequenceMassBinToScanNumsMap.TryGetValue(binNum, out scanNums))
                        {
                            scanNums.Add(ms2ScanNum);
                        }
                        else
                        {
                            scanNums = new List<int> { ms2ScanNum };
                            _sequenceMassBinToScanNumsMap.Add(binNum, scanNums);
                        }
                    }
                }
            }
        }

        public IEnumerable<double> GetPossibleSequenceMasses(int ms2ScanNumber, int numPeaks)
        {
            var productSpec = _run.GetSpectrum(ms2ScanNumber) as ProductSpectrum;
            if (productSpec == null) return null;

            var isolationWindow = productSpec.IsolationWindow;
            var minMz = isolationWindow.MinMz;
            var maxMz = isolationWindow.MaxMz;

            var monoIsotopicPeaks = new List<DeisotopedPeak>();
            var precursorSpec = _run.GetSpectrum(_run.GetPrecursorScanNum(ms2ScanNumber));
            if (precursorSpec != null)
            {
                var peakList = precursorSpec.GetPeakListWithin(minMz, maxMz);
                var specWindow = precursorSpec.GetPeakListWithin(minMz - 1.0, maxMz + 1.0);
                monoIsotopicPeaks.AddRange(GetDeisotopedPeaks(specWindow, peakList, numPeaks));
            }

            var nextScanNum = _run.GetNextScanNum(ms2ScanNumber, 1);
            var nextSpec = _run.GetSpectrum(nextScanNum);
            if (nextSpec != null)
            {
                var peakList = nextSpec.GetPeakListWithin(minMz, maxMz);
                var specWindow = nextSpec.GetPeakListWithin(minMz - 1.0, maxMz + 1.0);
                monoIsotopicPeaks.AddRange(GetDeisotopedPeaks(specWindow, peakList, numPeaks));
            }

            return monoIsotopicPeaks.GroupBy(p => p.MonoIsotopicMassIndex)
                .Select(g => g.OrderByDescending(p => p.Score).First())
                .OrderByDescending(p => p.Score)
                .Take(numPeaks)
                .Select(p => p.MonoIsotopicMass);
        }

        private IEnumerable<DeisotopedPeak> GetDeisotopedPeaks(List<Peak> specWindow, IEnumerable<Peak> peakList, int numDeisotopedPeaksToGet)
        {
            var peakListSortedByIntensity = new List<Peak>(peakList);
            peakListSortedByIntensity.Sort(new IntensityComparer());
            var remainingPeakList = new LinkedList<Peak>(peakListSortedByIntensity);

            var deisotopedPeakSet = new SortedSet<DeisotopedPeak>();
            while (remainingPeakList.Any())
            {
                var peakWithHighestIntensity = remainingPeakList.First.Value;
                var peakMz = peakWithHighestIntensity.Mz;
                var score = new double[_maxCharge+1];
                for (var charge = _maxCharge; charge >= _minCharge; charge--)
                {
                    // check whether next isotope peak exists
                    var nextIsotopeMz = peakMz + Constants.C13MinusC12 / charge;
                    var nextIsotopePeak = PeakListUtils.FindPeak(specWindow, nextIsotopeMz, _tolerance);
                    if (nextIsotopePeak == null) continue;
                    var mostAbundantIsotopeMass = (peakMz - Constants.Proton) * charge;
                    var averagineIsoEnv = Averagine.GetIsotopomerEnvelope(mostAbundantIsotopeMass);
                    var approxMostAbundantIsotopeIndex = averagineIsoEnv.MostAbundantIsotopeIndex;
                    var monoIsotopicMass = mostAbundantIsotopeMass - approxMostAbundantIsotopeIndex * Constants.C13MinusC12;
                    var averagineIsotopeProfile = Averagine.GetTheoreticalIsotopeProfile(monoIsotopicMass, charge);
                    var corr = PeakListUtils.GetPearsonCorrelation(specWindow, averagineIsotopeProfile, _comparer);

                    score[charge] = corr;
                    var isValid = true;
                    for (var mult = 2; mult <= _maxCharge/charge; mult++)
                    {
                        var multiple = charge*mult;
                        if (score[multiple] > 0.8*corr)
                        {
                            isValid = false;
                            break;
                        }
                    }
                    if (!isValid) continue;
                    deisotopedPeakSet.Add(new DeisotopedPeak(monoIsotopicMass, charge, corr));
                    if (deisotopedPeakSet.Count > numDeisotopedPeaksToGet)
                    {
                        deisotopedPeakSet.Remove(deisotopedPeakSet.Min);
                    }
                }
                remainingPeakList.RemoveFirst();
            }
            return deisotopedPeakSet;
        }

        internal class DeisotopedPeak: IComparable<DeisotopedPeak>
        {
            internal DeisotopedPeak(double monoIsotopicMass, int charge, double score)
            {
                MonoIsotopicMass = monoIsotopicMass;
                MonoIsotopicMassIndex = ProductScorerBasedOnDeconvolutedSpectra.GetBinNumber(MonoIsotopicMass);
                Charge = charge;
                Score = score;
            }
            internal double MonoIsotopicMass { get; private set; }
            internal int MonoIsotopicMassIndex { get; private set; }
            internal int Charge { get; private set; }
            internal double Score { get; private set; }

            public int CompareTo(DeisotopedPeak other)
            {
                var comp = Score.CompareTo(other.Score);
                if (comp != 0) return comp;
                return MonoIsotopicMassIndex.CompareTo(other.MonoIsotopicMassIndex);
            }
        }
    }
}
