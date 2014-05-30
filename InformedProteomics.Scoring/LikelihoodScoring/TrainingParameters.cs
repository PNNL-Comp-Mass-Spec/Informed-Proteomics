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
            _charges = new HashSet<int>();
            _massSorter = new Dictionary<int, Histogram<double>>();
            _rankTables = new Dictionary<int, List<RankTable>>();
            _drankTables = new Dictionary<int, List<RankTable>>();
            _ionProbabilities = new Dictionary<int, List<ProductIonFrequencyTable>>();
            _massErrors = new Dictionary<int, List<MassErrorTable>>();
            _precursorOffsets = new Dictionary<int, List<PrecursorOffsets>>();
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
            _ionProbabilities = new Dictionary<int, List<ProductIonFrequencyTable>>();
            _massErrors = new Dictionary<int, List<MassErrorTable>>();
            _precursorOffsets = new Dictionary<int, List<PrecursorOffsets>>();
            Config = config;
            ReadFromFile(fileName);
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

        public ProductIonFrequencyTable GetIonProbabilityTable(int charge, double mass)
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
        /// <param name="decoys">Decoy Peptide-Spectrum matches ordered the same as the target list.
        /// </param>
        public void AddMatches(SpectrumMatchList targets, SpectrumMatchList decoys)
        {
            if (targets.Count != decoys.Count) throw new ArgumentException("Unequal length target and decoy lists.");
            for (int i = 0; i < targets.Count; i++)
            {
                var target = targets[i];
                var decoy = decoys[i];
                var charge = target.PrecursorCharge;
                if (!_massSorter.ContainsKey(charge))
                {
                    Initialize(charge, targets);
                    _charges.Add(charge);
                }
                var mass = target.PrecursorComposition.Mass;
                var index = _massSorter[charge].GetBinIndex(mass);
                AddMatch(target, charge, index, false);
                AddMatch(decoy, charge, index, true);
            }
        }

        public void WriteToFile(string fileName)
        {
            using (var file = new StreamWriter(fileName))
            {
                foreach (var charge in _charges)
                {
                    file.WriteLine("Charge\t"+charge);
                    file.Write("BinSize\t");
                    var massBins = _massSorter[charge].Bins;
                    for (int i = 0; i < Config.MassBins; i++)
                    {
                        file.Write(massBins[i].Count);
                        if (i < Config.MassBins-1) file.Write("\t");
                    }
                    file.WriteLine();
                    for (int i = 0; i < Config.MassBins; i++)
                    {
                        file.Write("BinEdges\t{0}", _massSorter[charge].BinEdges[i]);
                        var max = Double.PositiveInfinity;
                        if (i < Config.MassBins-1) max = _massSorter[charge].BinEdges[i+1];
                        file.Write("\t"+max);
                        file.WriteLine();
                        var ionTypes = _ionProbabilities[charge][i].SelectIonTypess(Config.SelectedIonThreshold);
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

        private void ReadFromFile(string fileName)
        {
            var ionTypeFactory = new IonTypeFactory(2);
            using (var file = new StreamReader(fileName))
            {
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
                            _ionProbabilities.Add(charge, new List<ProductIonFrequencyTable>());
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
        }

        private void AddMatch(SpectrumMatch match, int charge, int massIndex, bool isDecoy)
        {
            var massErrors = _massErrors[charge][massIndex];
            var rankTable = (isDecoy ? _drankTables[charge][massIndex] : _rankTables[charge][massIndex]);
            var ionFrequencies = _ionProbabilities[charge][massIndex];

            var acMatch = match;
            if (Config.AcquisitionMethod == AcquisitionMethod.Dia)
            {
                // filter out all peaks except ion peaks
                var ionPeakSpectrum = SpectrumFilter.FilterIonPeaks(match.Sequence, match.Spectrum,
                                                                    Config.IonTypes, Config.Tolerance);
                acMatch = new SpectrumMatch(acMatch.Sequence, ionPeakSpectrum, acMatch.PrecursorCharge, acMatch.Decoy);
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

        private void Initialize(int charge, SpectrumMatchList matches)
        {
            var chargeMatches = matches.GetCharge(charge);
            _massSorter.Add(charge, new Histogram<double>());
            foreach (var match in chargeMatches)
            {
                _massSorter[charge].AddDatum(match.PrecursorComposition.Mass);
            }
            _massSorter[charge].Equalize(Config.MassBins, 0);

            _rankTables.Add(charge, new List<RankTable>());
            _drankTables.Add(charge, new List<RankTable>());
            _ionProbabilities.Add(charge, new List<ProductIonFrequencyTable>());
            _massErrors.Add(charge, new List<MassErrorTable>());
            _precursorOffsets.Add(charge, new List<PrecursorOffsets>());
            for (int i = 0; i < _massSorter[charge].BinEdges.Length; i++)
            {
                _rankTables[charge].Add(new RankTable(Config.IonTypes, Config.Tolerance, Config.MaxRanks));
                _drankTables[charge].Add(new RankTable(Config.IonTypes, Config.Tolerance, Config.MaxRanks));
                _ionProbabilities[charge].Add(new ProductIonFrequencyTable(Config.IonTypes, Config.Tolerance, Config.RelativeIntensityThreshold));
                _massErrors[charge].Add(new MassErrorTable(Config.IonTypes, Config.Tolerance));
                _precursorOffsets[charge].Add(new PrecursorOffsets(charge, Config.PrecursorOffsetWidth, Config.PrecursorOffsetThreshold));
            }
        }

        private int GetCharge(int charge)
        {
            if (_charges.Count == 0) throw new Exception("Empty training set.");
            int key = charge;
            int maxCharge = _charges.Max();
            while (!_charges.Contains(key)) key = (key + 1) % maxCharge;
            return key;
        }

        private readonly HashSet<int> _charges;
        private readonly Dictionary<int, List<RankTable>> _rankTables;
        private readonly Dictionary<int, List<RankTable>> _drankTables;
        private readonly Dictionary<int, List<ProductIonFrequencyTable>> _ionProbabilities;
        private readonly Dictionary<int, List<MassErrorTable>> _massErrors;
        private readonly Dictionary<int, List<PrecursorOffsets>> _precursorOffsets; 
        private readonly Dictionary<int, Histogram<double>> _massSorter;

        #endregion
    }
}
