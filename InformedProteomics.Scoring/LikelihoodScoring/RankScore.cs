using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class RankScore
    {
        public RankScore(string fileName, int maxCharge=2)
        {
            _trainingSets = new Dictionary<int, Tuple<RankTable, RankTable>>();
            ReadTrainingSetFromFile(fileName, maxCharge);
        }

        public double GetScore(IonType ionType, int rankNum, int charge)
        {
            charge = GetCharge(charge);
            var rank = _trainingSets[charge].Item1.GetRank(rankNum);
            var decoyRank = _trainingSets[charge].Item2.GetRank(rankNum);
            var prob = rank.IonFrequencies[ionType];
            var decoyProb = decoyRank.IonFrequencies[ionType];

            //if (prob == 0) Console.WriteLine("*** Target: {0} {1}", ionType.Name, rankNum);
            //if (decoyProb == 0) Console.WriteLine("*** Decoy: {0} {1}", ionType.Name, rankNum);
            if (Math.Abs(prob) < 1e-10 || Math.Abs(decoyProb) < 1e-10) return 0;
            return (Math.Log(prob) - Math.Log(decoyProb));
        }

        public IonType[] GetIonTypes(int charge)
        {
            charge = GetCharge(charge);
            return _trainingSets[charge].Item1.GetBinEdges();
        }

        //public IonType[] IonTypes { get { return _trainingSets[GetCharge(0)].Item1.GetBinEdges();  }}

        private int GetCharge(int charge)
        {
            var keys = _trainingSets.Keys;
            if (keys.Count == 0) throw new Exception("Empty training set.");
            int key = charge;
            int maxCharge = keys.Max();
            while (!_trainingSets.ContainsKey(key)) key = (key + 1) % maxCharge;
            return key;
        }

        private void ReadTrainingSetFromFile(string fileName, int maxCharge)
        {
            var ionTypeFactory = new IonTypeFactory(maxCharge);
            int currCharge = 0;
            int ranks = 0;
            var targetProb = new List<RankProbability>();
            var decoyProb = new List<RankProbability>();
            var file = File.ReadAllLines(fileName).ToList();
            file.Add("Charge\t0");
            var ionTypes = new List<IonType>();
            foreach (var line in file)
            {
                var parts = line.Split('\t').ToList();
                if (parts.Count < 2) throw new FormatException("Invalid line in training file: "+line);
                var header = parts[0];
                parts.RemoveAt(0); 
                if (header == "Charge")
                {
                    if (currCharge != 0)
                    {
                        var unfoundProb = targetProb[ranks];
                        var decoyunfoundProb = decoyProb[ranks];
                        targetProb.RemoveAt(ranks);
                        decoyProb.RemoveAt(ranks);
                        var tarTable = new RankTable(targetProb, ionTypes.ToArray(), null, ranks) { NotFound = unfoundProb };
                        var decTable = new RankTable(decoyProb, ionTypes.ToArray(), null, ranks) { NotFound = decoyunfoundProb };
                        _trainingSets.Add(currCharge, new Tuple<RankTable, RankTable>(tarTable, decTable));
                        ionTypes = new List<IonType>();
                    }
                    currCharge = Convert.ToInt32(parts[0]);
                }
                else if (header == "Ranks")
                {
                    ranks = Convert.ToInt32(parts[0]);
                    targetProb = new List<RankProbability> { Capacity = ranks+1 };
                    decoyProb = new List<RankProbability> { Capacity = ranks+1 };
                    for (int i = 0; i < ranks+1; i++)
                    {
                        targetProb.Add(new RankProbability());
                        decoyProb.Add(new RankProbability());
                    }
                }
                else
                {
                    int index = header.IndexOf("-Decoy", StringComparison.Ordinal);
                    if (index < 0)
                    {
                        try
                        {
                            var ion = ionTypeFactory.GetIonType(header);
                            ionTypes.Add(ion);
                            for (int i = 0; i < ranks + 1; i++)
                                targetProb[i].AddIon(ion, Convert.ToDouble(parts[i]));
                        }
                        catch (KeyNotFoundException)
                        {
                            if (header == "Total")
                            {
                                for (int i = 0; i < ranks + 1; i++)
                                {
                                    targetProb[i].RankCount = Convert.ToInt32(parts[i]);
                                    decoyProb[i].RankCount = Convert.ToInt32(parts[i]);
                                }
                            }
                        }
                    }
                    else
                    {
                        var ionName = header.Substring(0, index);
                        try
                        {
                            var ion = ionTypeFactory.GetIonType(ionName);
                            for (int i = 0; i < ranks + 1; i++) decoyProb[i].AddIon(ion, Convert.ToDouble(parts[i]));
                        }
                        catch (KeyNotFoundException) {}
                    }
                }

            }
        }

        private readonly Dictionary<int, Tuple<RankTable, RankTable>> _trainingSets;

/*        private void ReadTrainingSetFromFile(string fileName, int maxCharge)
        {
            var ionTypeFactory = new IonTypeFactory(maxCharge);
            var parser = new TsvFileParser(fileName);
            var data = parser.GetAllData();
            var ranks = parser.GetData("Rank");
            var totalRanks = ranks.Count;
            var ionTypes = (from key in data.Keys
                            where (key != "Rank" && key != "Unexplained" &&
                                   key != "Total" && !key.Contains("De"))
                            select ionTypeFactory.GetIonType(key)).ToArray();

            var rankProbabilities = new List<RankProbability>();
            var decoyRankProbabilities = new List<RankProbability>();
            for (var i = 0; i < totalRanks; i++)
            {
                var rankProb = new RankProbability();
                var decoyRankProb = new RankProbability();
                foreach (var ionType in ionTypes)
                {
                    var target = Convert.ToDouble(data[ionType.Name][i]);
                    var decoy = Convert.ToDouble(data[ionType.Name + "-De"][i]);
                    rankProb.IonFrequencies[ionType] = target;
                    decoyRankProb.IonFrequencies[ionType] = decoy;
                }
                var total = Convert.ToInt32(data["Total"][i]);
                rankProbabilities.Add(rankProb);
                rankProbabilities[i].RankCount = total;
                decoyRankProbabilities.Add(decoyRankProb);
                decoyRankProbabilities[i].RankCount = total;
            }
            var unfoundProb = rankProbabilities[ranks.Count - 1];
            var decoyunfoundProb = decoyRankProbabilities[ranks.Count - 1];
            rankProbabilities.RemoveAt(ranks.Count - 1);
            decoyRankProbabilities.RemoveAt(ranks.Count - 1);

            var rankTable = new RankTable(rankProbabilities, ionTypes, null, totalRanks-1) {NotFound = unfoundProb};
            var decoyRankTable = new RankTable(decoyRankProbabilities, ionTypes, null, totalRanks-1)
                                { NotFound = decoyunfoundProb };
            _trainingSet = rankTable;
            _decoySet = decoyRankTable;
        }

        private RankTable _trainingSet;
        private RankTable _decoySet;*/
    }
}
