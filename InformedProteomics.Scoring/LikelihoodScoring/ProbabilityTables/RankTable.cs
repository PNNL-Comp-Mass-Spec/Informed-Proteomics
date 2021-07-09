using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables
{
    public class RankTable
    {
        public int MaxRanks { get; private set; }
        public IonType[] IonTypes => _rankTable.Keys.ToArray();

        /// <summary>
        /// RankTable constructor for new RankTable without any previously calculated data.
        /// </summary>
        /// <param name="ionTypes">Ion types to calculate rank probabilities for.</param>
        /// <param name="tolerance">Tolerance for finding monoisotopic peak.</param>
        /// <param name="maxRanks">Maximum number of ranks in the RankTable. All ranks
        /// above maxRanks are put into the last rank.</param>
        public RankTable(IEnumerable<IonType> ionTypes, Tolerance tolerance, int maxRanks)
        {
            _tolerance = tolerance;
            MaxRanks = maxRanks;
            _rankTable = new Dictionary<IonType, double[]>();
            _rankTotals = new double[MaxRanks + 1];
            for (var i = 0; i < MaxRanks; i++)
            {
                _rankTotals[i] = 0.0;
            }

            foreach (var ionType in ionTypes)
            {
                _rankTable.Add(ionType, new double[MaxRanks + 1]);
                for (var i = 0; i < MaxRanks + 1; i++)
                {
                    _rankTable[ionType][i] = 0.0;
                }
            }
        }

        /// <summary>
        /// RankTable constructor for creating RankTable from training parameter file.
        /// </summary>
        /// <param name="file">Training data file past position of RankProbabilities label.</param>
        /// <param name="ionTypeFactory">IonTypeFactory object with all known possible ions in training
        /// parameter file.</param>
        public RankTable(TextReader file, IonTypeFactory ionTypeFactory)
        {
            MaxRanks = -1;
            _rankTable = new Dictionary<IonType, double[]>();
            ReadFromFile(file, ionTypeFactory);
        }

        /// <summary>
        /// Add single Peptide-Spectrum match to RankTable.
        /// </summary>
        /// <param name="match"></param>
        public void AddMatch(SpectrumMatch match)
        {
            var ranks = new RankedPeaks(match.Spectrum);
            for (var i = 0; i < ranks.Peaks.Length; i++)
            {
                var index = i;
                if (index >= MaxRanks)
                {
                    index = MaxRanks - 1;
                }

                _rankTotals[index]++;
                _rankTotals[MaxRanks]++;
            }
            foreach (var ionType in IonTypes)
            {
                var ions = match.GetCleavageIons(ionType);
                foreach (var ion in ions)
                {
                    var rank = ranks.RankIon(ion, _tolerance);
                    var rankIndex = GetRankIndex(rank);
                    _rankTable[ionType][rankIndex]++;
                }
            }
        }

        /// <summary>
        /// Add a collection of Peptide-Spectrum matches to RankTable.
        /// </summary>
        /// <param name="matches"></param>
        public void AddMatches(List<SpectrumMatch> matches)
        {
            foreach (var match in matches)
            {
                AddMatch(match);
            }
        }

        /// <summary>
        /// Smooth ranks by taking average of neighboring ranks of a given window size.
        /// </summary>
        /// <param name="smoothingRanks">Array of rank positions to smooth around.</param>
        /// <param name="smoothingWindowSize">Array of smoothing window sizes corresponding to
        /// the ranks in smoothingRanks.</param>
        public void Smooth(int[] smoothingRanks, int[] smoothingWindowSize)
        {
            if (smoothingRanks.Length != smoothingWindowSize.Length)
            {
                throw new ArgumentException("Unequal length smoothingRanks and smoothingWindowSize.");
            }

            foreach (var ionType in IonTypes)
            {
                var index = 0;
                var window = smoothingWindowSize[index];
                var total = 0.0;
                var count = 0;
                var endRank = 0;

                for (var startRank = 0; startRank < MaxRanks;)
                {
                    if (startRank >= smoothingRanks[index])
                    {
                        index++;
                        window = smoothingWindowSize[index];
                    }
                    endRank = Math.Min(endRank + window + 1, MaxRanks);
                    for (var i = startRank; i < endRank; i++)
                    {
                        total += _rankTable[ionType][i];
                        count++;
                    }
                    var average = total / count;
                    if (average > 0.0)
                    {
                        for (var i = startRank; i < endRank; i++)
                        {
                            _rankTable[ionType][i] = average;
                        }
                        total = 0.0;
                        count = 0;
                        startRank += window + 1;
                    }
                    else if (endRank == MaxRanks && startRank > 0)
                    {
                        startRank--;
                        total = 0.0;
                        count = 0;
                        continue;
                    }
                    if (endRank == MaxRanks)
                    {
                        break;
                    }
                }
            }
        }

        /// <summary>
        /// Get the number of peaks found for a given rank and ion type.
        /// </summary>
        /// <param name="rankNum"></param>
        /// <param name="ionType"></param>
        /// <returns>Number of matching peaks</returns>
        public Probability<IonType> GetRankProbability(int rankNum, IonType ionType)
        {
            var rankIndex = GetRankIndex(rankNum);
            var prob = _rankTable[ionType][rankIndex];
            var total = _rankTotals[rankIndex];
            return new Probability<IonType>(ionType, prob, total);
        }

        /// <summary>
        /// Creates a two dimensional array of all rank probabilities.
        /// </summary>
        /// <returns>Two dimensional array where the first dimension are ranks
        /// and the second dimension are ion types.</returns>
        public Probability<IonType>[,] GetProbabilities()
        {
            var ionTypes = IonTypes;
            var ranks = new Probability<IonType>[ionTypes.Length, MaxRanks];
            for (var i = 0; i < ionTypes.Length; i++)
            {
                for (var j = 0; j < MaxRanks + 1; j++)
                {
                    var ionType = ionTypes[i];
                    ranks[i, j] = new Probability<IonType>(ionTypes[i], _rankTable[ionType][j], _rankTotals[j]);
                }
            }
            return ranks;
        }

        /// <summary>
        /// Write RankTable to training parameter file.
        /// </summary>
        /// <param name="file">Open file to write RankTable to.</param>
        /// <param name="selectedIonTypes">The ion types to show in the file.
        /// All other ion types are averaged and displayed as "unexplained"</param>
        public void WriteToFile(StreamWriter file, IonType[] selectedIonTypes)
        {
            file.WriteLine("MaxRanks\t" + MaxRanks);
            var unselectedIonTypes = IonTypes.Except(selectedIonTypes).ToList();
            foreach (var ionType in selectedIonTypes)
            {
                file.Write(ionType.Name + "\t");
                for (var i = 0; i < MaxRanks + 1; i++)
                {
                    file.Write(Math.Round(_rankTable[ionType][i], 2) + "\t");
                }
                file.WriteLine();
            }
            var total = new double[MaxRanks + 1];
            var count = new double[MaxRanks + 1];
            var prob = new double[MaxRanks + 1];
            for (var i = 0; i < MaxRanks + 1; i++)
            {
                total[i] = 0.0;
                count[i] = 0.0;
                prob[i] = 0.0;
            }
            foreach (var ionType in unselectedIonTypes)
            {
                for (var i = 0; i < MaxRanks + 1; i++)
                {
                    total[i] += _rankTable[ionType][i];
                    count[i]++;
                }
            }
            for (var i = 0; i < MaxRanks + 1; i++)
            {
                prob[i] = Math.Round(total[i] / count[i], 2);
            }

            file.WriteLine("Unexplained\t" + string.Join("\t", prob));
            file.WriteLine("Total\t" + string.Join("\t", _rankTotals));
        }

        private void ReadFromFile(TextReader file, IonTypeFactory ionTypeFactory)
        {
            var line = file.ReadLine();
            while (line != null)
            {
                var parts = line.Split('\t').ToList();
                var header = parts[0];
                parts.RemoveAt(0);
                if (header == "MaxRanks")
                {
                    MaxRanks = Convert.ToInt32(parts[0]);
                }
                else if (header == "Total")
                {
                    if (MaxRanks < 0)
                    {
                        throw new FormatException("Badly formatted rank data.");
                    }

                    _rankTotals = new double[MaxRanks + 1];
                    for (var i = 0; i < MaxRanks + 1; i++)
                    {
                        _rankTotals[i] = Convert.ToInt32(parts[i]);
                    }
                    // reached end of rank table entry
                    break;
                }
                else if (header == "Unexplained")
                { }
                else
                {
                    if (MaxRanks < 0)
                    {
                        throw new FormatException("Badly formatted rank data.");
                    }

                    try
                    {
                        var ionType = ionTypeFactory.GetIonType(header);
                        _rankTable.Add(ionType, new double[MaxRanks + 1]);
                        for (var i = 0; i < MaxRanks + 1; i++)
                        {
                            _rankTable[ionType][i] = Convert.ToDouble(parts[i]);
                        }
                    }
                    catch (KeyNotFoundException)
                    {
                        throw new FormatException("Invalid ion type: " + header);
                    }
                }
                line = file.ReadLine();
            }
        }

        private int GetRankIndex(int rankNum)
        {
            if (rankNum < 1)
            {
                rankNum = MaxRanks + 1;
            }
            else if (rankNum > MaxRanks)
            {
                rankNum = MaxRanks;
            }

            return rankNum - 1;
        }

        private readonly Dictionary<IonType, double[]> _rankTable;
        private double[] _rankTotals;
        private readonly Tolerance _tolerance;
    }
}
