using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Config;
using InformedProteomics.Scoring.LikelihoodScoring.Data;
using InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class TrainingParameters
    {
        public TrainerConfiguration Config { get; set; }

        #region Constructors
        /// <summary>
        /// Constructor for creating empty TrainingParameters object.
        /// </summary>
        /// <param name="config">Configuration data read from configuration file.</param>
        public TrainingParameters(TrainerConfiguration config)
        {
            Config = config;
            _dataSet = new Queue<SpectrumMatch>();
            _massBins = new Dictionary<int, int>();
            _charges = new HashSet<int>();
            _massSorter = new Dictionary<int, Histogram<double>>();
            _rankTables = new Dictionary<int, List<RankTable>>();
            _drankTables = new Dictionary<int, List<RankTable>>();
            _ionProbabilities = new Dictionary<int, List<IonFrequencyTable>>();
            _massErrors = new Dictionary<int, List<MassErrorTable>>();
            _precursorOffsets = new Dictionary<int, List<PrecursorOffsets>>();
            _computed = false;
        }

        /// <summary>
        /// Constructor for creating TrainingParameters from existing data in
        /// a training parameter file.
        /// </summary>
        /// <param name="fileName">Name of the training parameter file.</param>
        /// <param name="config">Configuration data read from configuration file. Optional.</param>
        public TrainingParameters(string fileName, TrainerConfiguration config=null)
        {
            _charges = new HashSet<int>();
            _massSorter = new Dictionary<int, Histogram<double>>();
            _rankTables = new Dictionary<int, List<RankTable>>();
            _drankTables = new Dictionary<int, List<RankTable>>();
            _ionProbabilities = new Dictionary<int, List<IonFrequencyTable>>();
            _massErrors = new Dictionary<int, List<MassErrorTable>>();
            _precursorOffsets = new Dictionary<int, List<PrecursorOffsets>>();
            Config = config;
            ReadFromFile(fileName);
            _computed = true;
        }

        /// <summary>
        /// Constructor for creating TrainingParameters from resource file stream parameter file.
        /// </summary>
        /// <param name="resourceStream"></param>
        /// <param name="config"></param>
        public TrainingParameters(Stream resourceStream, TrainerConfiguration config = null)
        {
            _charges = new HashSet<int>();
            _massSorter = new Dictionary<int, Histogram<double>>();
            _rankTables = new Dictionary<int, List<RankTable>>();
            _drankTables = new Dictionary<int, List<RankTable>>();
            _ionProbabilities = new Dictionary<int, List<IonFrequencyTable>>();
            _massErrors = new Dictionary<int, List<MassErrorTable>>();
            _precursorOffsets = new Dictionary<int, List<PrecursorOffsets>>();
            Config = config;
            ReadFromResourceFile(resourceStream);
            _computed = true;
        }
        #endregion

        #region Getters
        public IonType[] GetIonTypes(int charge, double mass)
        {
            charge = GetCharge(charge);
            var index = _massSorter[charge].GetBinIndex(mass);
            return _rankTables[charge][index].IonTypes;
        }

        public RankTable GetRankTable(int charge, double mass, bool isDecoy)
        {
            charge = GetCharge(charge);
            var index = _massSorter[charge].GetBinIndex(mass);
            return (isDecoy ? _drankTables[charge][index] : _rankTables[charge][index]);
        }

        public IonFrequencyTable GetIonProbabilityTable(int charge, double mass)
        {
            charge = GetCharge(charge);
            var index = _massSorter[charge].GetBinIndex(mass);
            return _ionProbabilities[charge][index];
        }

        public MassErrorTable GetMassErrorTable(int charge, double mass)
        {
            charge = GetCharge(charge);
            var index = _massSorter[charge].GetBinIndex(mass);
            return _massErrors[charge][index];
        }
        #endregion

        #region Methods
        /// <summary>
        /// Adds lists of targets and decoys to the RankTables, IonProbabilityTables,
        /// and MassErrorTables. Partitions them by mass and precursor charge.
        /// </summary>
        /// <param name="targets">Target Peptide-Spectrum matches.</param>
        public void AddMatches(SpectrumMatchList targets)
        {
            _computed = false;
            foreach (var match in targets)
            {
                var charge = match.PrecursorCharge;
                if (!_charges.Contains(match.PrecursorCharge))
                {
                    _charges.Add(charge);
                    _rankTables.Add(charge, new List<RankTable>());
                    _drankTables.Add(charge, new List<RankTable>());
                    _ionProbabilities.Add(charge, new List<IonFrequencyTable>());
                    _massErrors.Add(charge, new List<MassErrorTable>());
                    _massSorter.Add(charge, new Histogram<double>());
                    _precursorOffsets.Add(charge, new List<PrecursorOffsets>());
                }
                _massSorter[charge].AddDatum(match.PrecursorComposition.Mass);
                _dataSet.Enqueue(match);
            }
        }

        public void WriteToFile(string fileName)
        {
            Compute();
            using (var file = new StreamWriter(fileName))
            {
                foreach (var charge in _charges)
                {
                    file.WriteLine("Charge\t"+charge);
                    file.Write("BinSize\t");
                    var massBins = _massSorter[charge].Bins;
                    for (var i = 0; i < _massBins[charge]; i++)
                    {
                        file.Write(massBins[i].Count);
                        if (i < _massBins[charge]-1) file.Write("\t");
                    }
                    file.WriteLine();
                    for (var i = 0; i < _massBins[charge]; i++)
                    {
                        file.Write("BinEdges\t{0}", _massSorter[charge].BinEdges[i]);
                        var max = Double.PositiveInfinity;
                        if (i < _massBins[charge]-1) max = _massSorter[charge].BinEdges[i+1];
                        file.Write("\t"+max);
                        file.WriteLine();
                        var ionTypes = _ionProbabilities[charge][i].SelectIonTypes(Config.SelectedIonThreshold);
                        file.WriteLine("RankProbabilities");
                        _rankTables[charge][i].Smooth(Config.SmoothingRanks, Config.SmoothingWindowSize);
                        _rankTables[charge][i].WriteToFile(file, ionTypes);
                        _drankTables[charge][i].Smooth(Config.SmoothingRanks, Config.SmoothingWindowSize);
                        _drankTables[charge][i].WriteToFile(file, ionTypes);
                        file.WriteLine("IonProbabilities");
                        _ionProbabilities[charge][i].WriteToFile(file);
                    }
                }
            }
        }
        #endregion

        #region Private Methods
        private void Read(StreamReader file)
        {
            var ionTypeFactory = new IonTypeFactory(2);
            var massBins = new List<double>();
            var binEdges = new List<double>();
            var charge = 0;
            var line = file.ReadLine();
            while (line != null)
            {
                var parts = line.Split('\t').ToList();
                var header = parts[0];
                parts.RemoveAt(0);
                if (header == "BinSize") massBins.AddRange(parts.Select(Convert.ToDouble));
                else if (header == "BinEdges")
                {
                    binEdges.Add(Convert.ToDouble(parts[0]));
                }
                else if (header == "Charge")
                {
                    if (binEdges.Count > 0)
                    {
                        _massSorter[charge].BinEdges = binEdges.ToArray();
                        binEdges.Clear();
                    }
                    charge = Convert.ToInt32(parts[0]);
                    if (!_charges.Contains(charge))
                    {
                        _charges.Add(charge);
                        _massSorter.Add(charge, new Histogram<double>());
                        _rankTables.Add(charge, new List<RankTable>());
                        _drankTables.Add(charge, new List<RankTable>());
                        _ionProbabilities.Add(charge, new List<IonFrequencyTable>());
                        _massErrors.Add(charge, new List<MassErrorTable>());
                    }
                }
                else if (header == "RankProbabilities")
                {
                    if (massBins.Count == 0 || charge == 0) throw new FormatException("Badly formatted training file.");
                    _rankTables[charge].Add(new RankTable(file, ionTypeFactory));
                    _drankTables[charge].Add(new RankTable(file, ionTypeFactory));
                }
                line = file.ReadLine();
            }
            _massSorter[charge].BinEdges = binEdges.ToArray();
        }

        private void ReadFromResourceFile(Stream resourceFile)
        {
            using (var fileStreamReader = new StreamReader(resourceFile))
            {
                Read(fileStreamReader);
            }
        }

        private void ReadFromFile(string fileName)
        {
            using (var fileStreamReader = new StreamReader(fileName))
            {
                Read(fileStreamReader);
            }
        }

        private void Compute()
        {
            if (_computed) return;
            Initialize();
            var total = _dataSet.Count;
            for (var i = 0; i < total; i++)
            {
                var target = _dataSet.Dequeue();
                var charge = target.PrecursorCharge;
                if (!_charges.Contains(charge)) continue;
                var decoy = new SpectrumMatch(target, true);
                var mass = target.PrecursorComposition.Mass;
                var index = _massSorter[charge].GetBinIndex(mass);
                ComputeMatch(target, charge, index);
                ComputeMatch(decoy, charge, index);
            }
            _computed = true;
        }

        private void ComputeMatch(SpectrumMatch match, int charge, int massIndex)
        {
            var isDecoy = match.Decoy;
            var massErrors = _massErrors[charge][massIndex];
            var rankTable = (isDecoy ? _drankTables[charge][massIndex] : _rankTables[charge][massIndex]);
            var ionFrequencies = _ionProbabilities[charge][massIndex];

            var acMatch = match;
            if (Config.AcquisitionMethod == AcquisitionMethod.Dia)
            {
                // filter out all peaks except ion peaks
                var ionPeakSpectrum = SpectrumFilter.FilterIonPeaks(match.Sequence, match.Spectrum,
                                                                    Config.IonTypes, Config.Tolerance);
                acMatch = new SpectrumMatch(acMatch.Sequence, ionPeakSpectrum, acMatch.ScanNum, acMatch.PrecursorCharge, acMatch.Decoy);
            }
            else
            {
                _precursorOffsets[charge][massIndex].AddMatch(match);
                // filter precursor peaks
                var precursorFilter = new PrecursorFilter(_precursorOffsets[charge][massIndex],
                    Config.MassErrorTolerance);
                acMatch = precursorFilter.Filter(acMatch);
            }
            rankTable.AddMatch(acMatch);
            if (isDecoy) return;

            massErrors.AddMatch(match);

            var filteredSpectrum = SpectrumFilter.GetFilteredSpectrum(match.Spectrum, Config.WindowWidth,
                                                                      Config.RetentionCount);
            var filteredMatch = new SpectrumMatch(match.Sequence, filteredSpectrum, match.PrecursorCharge);
            ionFrequencies.AddMatch(filteredMatch);
        }

        private void Initialize()
        {
            var validCharges = new HashSet<int>();
            // partition by mass
            foreach (var charge in _charges)
            {
                var chargeTableSize = _massSorter[charge].Total;
                if (chargeTableSize < Config.MassBinSize)
                {
                    continue;
//                    throw new Exception("Fewer than "+Config.MassBinSize+" items in charge "+charge+" data set.");
                }
                validCharges.Add(charge);
                _massBins.Add(charge, chargeTableSize/Config.MassBinSize);
                _massSorter[charge].Equalize(_massBins[charge], 0);
                var binCount = _massSorter[charge].BinEdges.Length;
                for (var i = 0; i < binCount; i++)
                {
                    _rankTables[charge].Add(new RankTable(Config.IonTypes, Config.Tolerance, Config.MaxRanks));
                    _drankTables[charge].Add(new RankTable(Config.IonTypes, Config.Tolerance, Config.MaxRanks));
                    _ionProbabilities[charge].Add(new IonFrequencyTable(Config.IonTypes, Config.Tolerance, Config.RelativeIntensityThreshold));
                    _massErrors[charge].Add(new MassErrorTable(Config.IonTypes, Config.Tolerance));
                    _precursorOffsets[charge].Add(new PrecursorOffsets(charge, Config.PrecursorOffsetWidth, Config.PrecursorOffsetThreshold));
                }
            }
            _charges = validCharges;
        }

        private int GetCharge(int charge)
        {
            if (_charges.Count == 0) throw new Exception("Empty training set.");
            var key = charge;
            var maxCharge = _charges.Max();
            while (!_charges.Contains(key)) key = (key + 1) % maxCharge;
            return key;
        }

        private readonly Queue<SpectrumMatch> _dataSet;
        private readonly Dictionary<int, int> _massBins;
        private HashSet<int> _charges;
        private readonly Dictionary<int, List<RankTable>> _rankTables;
        private readonly Dictionary<int, List<RankTable>> _drankTables;
        private readonly Dictionary<int, List<IonFrequencyTable>> _ionProbabilities;
        private readonly Dictionary<int, List<MassErrorTable>> _massErrors;
        private readonly Dictionary<int, List<PrecursorOffsets>> _precursorOffsets;
        private readonly Dictionary<int, Histogram<double>> _massSorter;
        private bool _computed;

        #endregion
    }
}
