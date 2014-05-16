using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class RankScore
    {
        public RankScore(string fileName, int maxCharge=2)
        {
            ReadTrainingSetFromFile(fileName, maxCharge);
        }

        public double GetScore(IonType ionType, int rankNum)
        {
            var rank = _trainingSet.GetRank(rankNum);
            var decoyRank = _decoySet.GetRank(rankNum);
            var prob = rank.IonFrequencies[ionType];
            var decoyProb = decoyRank.IonFrequencies[ionType];

            //if (prob == 0) Console.WriteLine("*** Target: {0} {1}", ionType.Name, rankNum);
            //if (decoyProb == 0) Console.WriteLine("*** Decoy: {0} {1}", ionType.Name, rankNum);
            if (Math.Abs(prob) < 1e-10 || Math.Abs(decoyProb) < 1e-10) return 0;
            return (Math.Log(prob) - Math.Log(decoyProb));
        }

        public IonType[] IonTypes
        {
            get { return _trainingSet.GetBinEdges(); }
        }

        private void ReadTrainingSetFromFile(string fileName, int maxCharge)
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
                var rankProb = new RankProbability(ionTypes);
                var decoyRankProb = new RankProbability(ionTypes);
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
        private RankTable _decoySet;
    }
}
