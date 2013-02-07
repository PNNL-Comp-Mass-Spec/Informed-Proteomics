using System;
using InformedProteomics.Backend.Data.Results;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Scoring
{
    public class Scorer : IScorer
    {
        public float PrecursorIonScore { get; private set; }
        public float ProductIonScore { get; private set; }

        public Scorer(Sequence seq, DatabaseMultipleSubTargetResult matchedResult)
        {
            Seq = seq;
            MatchedResult = TrimXYData(matchedResult);
            //PrecursorIonScore = new PrecursorIonScorer(MatchedResult).Score;
            //ProductIonScore = new ProductIonScorer(MatchedResult).Score;
            //Console.WriteLine(PrecursorIonScore + "\t" + ProductIonScore);
            //PrecursorIonScore just use to choose best XIC. Then just use ProductIonScore... imp next.. TODO
            Score = ProductIonScore;
        }

        public Sequence Seq
        {
            get;
            private set;
        }

        public DatabaseMultipleSubTargetResult MatchedResult
        {
            get;
            private set;
        }

        public float Score
        {
            get;
            private set;
        }

        private DatabaseMultipleSubTargetResult TrimXYData(DatabaseMultipleSubTargetResult result)
        {
            var refXYData = result.PrecursorResultRep.XYData;
            foreach (var r in result.SubTargetResultList)
            {
                var aligned = DataUtil.AlignXYData(refXYData, r.XYData, 0);
                r.XYData.Xvalues = aligned.Xvalues;
                r.XYData.Yvalues = aligned.Yvalues;
            }
            foreach (var r in result.FragmentResultList)
            {
                var aligned = DataUtil.AlignXYData(refXYData, r.XYData, 1);
                r.XYData.Xvalues = aligned.Xvalues;
                r.XYData.Yvalues = aligned.Yvalues;
            }
            return result;
        }
    }
}
