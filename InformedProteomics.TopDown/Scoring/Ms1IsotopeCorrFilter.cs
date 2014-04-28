using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.Scoring
{
    public class Ms1IsotopeCorrFilter: ISequenceFilter
    {
        public Ms1IsotopeCorrFilter(
            LcMsRun run, 
            int minCharge = 3, int maxCharge = 30,
            double ppmTolerance = 15,
            double corrThreshold = 0.5,
            int maxNumDeisotopedPeaksPerIsolationWindow = 40)
        {
            _run = run;
            _minCharge = minCharge;
            _maxCharge = maxCharge;
            _corrThreshold = corrThreshold;
            MaxNumDeisotopedPeaksPerIsolationWindow = maxNumDeisotopedPeaksPerIsolationWindow;
            _ppmTolerance = ppmTolerance;
            _tolerance = new Tolerance(_ppmTolerance);
            _comparer = new MzComparerWithPpmTolerance(ppmTolerance);
            _sequenceMassBinToScanNumsMap = new Dictionary<int, IList<int>>();
//            _topKFilter = new Ms1IsotopeTopKFilter(run, minCharge, maxCharge, ppmTolerance);
            PrecomputePossibleSequenceMasses();
        }

        public IEnumerable<int> GetMatchingMs2ScanNums(double sequenceMass)
        {
            var sequenceMassBinNum = ProductScorerBasedOnDeconvolutedSpectra.GetBinNumber(sequenceMass);
            IList<int> ms2ScanNums;
            if (_sequenceMassBinToScanNumsMap.TryGetValue(sequenceMassBinNum, out ms2ScanNums)) return ms2ScanNums;
            return new int[0];
        }

        public int MaxNumDeisotopedPeaksPerIsolationWindow { get; private set; }

        private readonly LcMsRun _run;
//        private readonly Ms1IsotopeTopKFilter _topKFilter;
        private readonly int _minCharge;
        private readonly int _maxCharge;
        private readonly double _corrThreshold;
        private readonly double _ppmTolerance;
        private readonly Tolerance _tolerance;
        private readonly MzComparerWithPpmTolerance _comparer;

        private readonly Dictionary<int, IList<int>> _sequenceMassBinToScanNumsMap; // mass -> scan numbers

        private void PrecomputePossibleSequenceMasses()
        {
            var ms2ScanNums = _run.GetScanNumbers(2);
            foreach (var ms2ScanNum in ms2ScanNums)
            {
                var possibleSequenceMasses = GetPossibleSequenceMasses(ms2ScanNum);
                var binNumHash = new HashSet<int>();
                foreach (var sequenceMass in possibleSequenceMasses)
                {
                    var dMass = sequenceMass*_ppmTolerance*1e-6;
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

        private static readonly MzComparerWithPpmTolerance SumComparer = new MzComparerWithPpmTolerance(5);
        public IEnumerable<double> GetPossibleSequenceMasses(int ms2ScanNumber)
        {
            var productSpec = _run.GetSpectrum(ms2ScanNumber) as ProductSpectrum;
            if (productSpec == null) return null;

            var isolationWindow = productSpec.IsolationWindow;
            var minMz = isolationWindow.MinMz;
            var maxMz = isolationWindow.MaxMz;

            var sequenceMassList = new List<double>();

            var summedSpec = _run.GetSummedMs1Spectrum(ms2ScanNumber, 4);
            if (summedSpec != null)
            {
                var peakList = summedSpec.GetPeakListWithin(minMz, maxMz);
                var specWindow = summedSpec.GetPeakListWithin(minMz - 1.0, maxMz + 1.0);
                ApplyDeconvolution(specWindow, peakList, ref sequenceMassList, MaxNumDeisotopedPeaksPerIsolationWindow);
            }

            //var precursorSpec = _run.GetSpectrum(_run.GetPrecursorScanNum(ms2ScanNumber));
            //if (precursorSpec != null)
            //{
            //    var peakList = precursorSpec.GetPeakListWithin(minMz, maxMz);
            //    var specWindow = precursorSpec.GetPeakListWithin(minMz - 1.0, maxMz + 1.0);
            //    ApplyDeconvolution(specWindow, peakList, ref sequenceMassList, maxNumPeaksToConsider);
            //}

            //var nextScanNum = _run.GetNextScanNum(ms2ScanNumber, 1);
            //var nextSpec = _run.GetSpectrum(nextScanNum);
            //if (nextSpec != null)
            //{
            //    var peakList = nextSpec.GetPeakListWithin(minMz, maxMz);
            //    var specWindow = nextSpec.GetPeakListWithin(minMz - 1.0, maxMz + 1.0);
            //    ApplyDeconvolution(specWindow, peakList, ref sequenceMassList, maxNumPeaksToConsider);
            //}

            return sequenceMassList;
        }

        private void ApplyDeconvolution(List<Peak> specWindow, List<Peak> peakList, ref List<double> deisotopedMassList, int numTrials)
        {
            var peakListSortedByIntensity = new List<Peak>(peakList);
            peakListSortedByIntensity.Sort(new IntensityComparer());
            var remainingPeakList = new LinkedList<Peak>(peakListSortedByIntensity);
            for (var i = 0; i < numTrials; i++)
            {
                ApplyDeconvolution(specWindow, ref remainingPeakList, ref deisotopedMassList);
            }
        }

        private void ApplyDeconvolution(List<Peak> specWindow, ref LinkedList<Peak> remainingPeakList, ref List<double> deisotopedMassList)
        {
            if (!remainingPeakList.Any()) return;

            var peakWithHighestIntensity = remainingPeakList.First.Value;
            var peakMz = peakWithHighestIntensity.Mz;
            var score = new double[_maxCharge + 1];
            for (var charge = _maxCharge; charge >= _minCharge; charge--)
            {
                // check whether next isotope peak exists
                var nextIsotopeMz = peakMz + Constants.C13MinusC12 / charge;
                var nextIsotopePeak = PeakListUtils.FindPeak(specWindow, nextIsotopeMz, _tolerance);
                if (nextIsotopePeak == null) continue;
                var mostAbundantIsotopeMass = (peakMz - Constants.Proton)*charge;
                var averagineIsoEnv = Averagine.GetIsotopomerEnvelope(mostAbundantIsotopeMass);
                var approxMostAbundantIsotopeIndex = averagineIsoEnv.MostAbundantIsotopeIndex;
                var monoIsotopicMass = mostAbundantIsotopeMass - approxMostAbundantIsotopeIndex*Constants.C13MinusC12;
                var averagineIsotopeProfile = Averagine.GetTheoreticalIsotopeProfile(monoIsotopicMass, charge);
                var corr = PeakListUtils.GetPearsonCorrelation(specWindow, averagineIsotopeProfile, _comparer);

                score[charge] = corr;
                var isValid = true;
                for (var mult = 2; mult <= _maxCharge / charge; mult++)
                {
                    var multiple = charge * mult;
                    if (score[multiple] > 0.8 * corr)
                    {
                        isValid = false;
                        break;
                    }
                }
                if (!isValid) continue;

                if(corr > _corrThreshold) deisotopedMassList.Add(monoIsotopicMass);
            }
            remainingPeakList.RemoveFirst();
        }

        //private void ApplyDeconvolutionOld(ref List<Peak> specWindow, ref List<double> deisotopedMassList, int numTrials)
        //{;
        //    for (var i = 0; i < numTrials; i++)
        //    {
        //        var isSuccess = ApplyDeconvolutionOld(ref specWindow, ref deisotopedMassList);
        //        if (!isSuccess) break;
        //    }
        //}

        //private bool ApplyDeconvolutionOld(ref List<Peak> specWindow, ref List<double> deisotopedMassList)
        //{
        //    // Find the peak with the maximum intensity
        //    var maxIntensityPeak = new Peak(0, 0);
        //    foreach (var p in specWindow)
        //    {
        //        if (p.Intensity > maxIntensityPeak.Intensity) maxIntensityPeak = p;
        //    }
        //    if (maxIntensityPeak.Intensity <= 0) return false;

        //    var bestCharge = -1;
        //    var selectedMonoIsotopeMass = 0.0;
        //    var bestScoreForAllCharges = _corrThreshold;
        //    List<double> deconvolutedMassesForAllCharges = null;

        //    var peakMz = maxIntensityPeak.Mz;
        //    for (var charge = _maxCharge; charge >= _minCharge; charge--)
        //    {
        //        if (bestCharge%charge == 0) continue;   // skip lower charges if its multiple has a good score
        //        var mass = maxIntensityPeak.Mz * charge;
        //        var mostAbundantIsotopeIndex = Averagine.GetIsotopomerEnvelope(mass).MostAbundantIsotopeIndex;
        //        var bestScore = 0.0;
        //        var bestMonoMass = 0.0;
        //        var deconvolutedMasses = new List<double>();

        //        for (var isotopeIndex = mostAbundantIsotopeIndex - _isotopeIndexOffsetThreshold; isotopeIndex <= mostAbundantIsotopeIndex + _isotopeIndexOffsetThreshold; isotopeIndex++)
        //        {
        //            var isotopeMz = peakMz - (isotopeIndex - mostAbundantIsotopeIndex)*Constants.C13MinusC12/charge;
        //            var deltaMass = isotopeMz*_ppmTolerance*1e-6;
        //            var peak = PeakListUtils.FindPeak(specWindow, isotopeMz-deltaMass, isotopeMz+deltaMass);
        //            if(peak == null) continue;
        //            var monoIsotopeMass = (peak.Mz - Constants.Proton) * charge - mostAbundantIsotopeIndex * Constants.C13MinusC12;
        //            if (monoIsotopeMass <= 0) continue;
        //            var theoIsoProfile = Averagine.GetTheoreticalIsotopeProfile(monoIsotopeMass, charge);
        //            var corr = PeakListUtils.GetPearsonCorrelation(specWindow, theoIsoProfile, _comparer);
        //            if (corr < _corrThreshold) continue;

        //            if (corr > bestScore)
        //            {
        //                bestScore = corr;
        //                bestMonoMass = monoIsotopeMass;
        //            }
        //            deconvolutedMasses.Add(monoIsotopeMass);
        //        }

        //        if (bestScore > bestScoreForAllCharges)
        //        {
        //            bestCharge = charge;
        //            selectedMonoIsotopeMass = bestMonoMass;
        //            bestScoreForAllCharges = bestScore;
        //            deconvolutedMassesForAllCharges = deconvolutedMasses;
        //        }
        //    }

        //    if (bestCharge < 0) return false;

        //    var selectedTheoIsoProfile = Averagine.GetTheoreticalIsotopeProfile(selectedMonoIsotopeMass, bestCharge);
        //    specWindow = PeakListUtils.GetExceptWith(specWindow, selectedTheoIsoProfile, _comparer);
        //    if(deconvolutedMassesForAllCharges != null) deisotopedMassList.AddRange(deconvolutedMassesForAllCharges);

        //    return true;
        //}
    }
}
