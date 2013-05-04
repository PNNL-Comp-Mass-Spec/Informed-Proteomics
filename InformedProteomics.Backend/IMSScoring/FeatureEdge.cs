using System;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FeatureEdge
    {
        public FeatureNode LNode { get; private set; }
        public FeatureNode RNode { get; private set; }
        public double Weight { get; private set; } // used to calculate weight of a path
        public double Score { get; private set; } // used to calculate score
        private readonly int _ratio;
        private readonly double _lcCorrelation;
        private readonly double _imsCorrelation;

        public FeatureEdge(FeatureNode l, FeatureNode r)
        {
            LNode = l;
            RNode = r;
            var ri = r.Feature == null ? 0.0 : r.Feature.IntensityMax;
            var li = l.Feature == null ? 0.0 : l.Feature.IntensityMax;
            
            if (LNode is PrecursorFeatureNode) // TODO fix later when no summed intensity is used
                _ratio = GetRatioScore(ri, ri);
            else
                _ratio = GetRatioScore(li, ri);
            _lcCorrelation = StatisticsTools.GetLCCorrelation(l.Feature, r.Feature);
            _imsCorrelation = StatisticsTools.GetIMSCorrelation(l.Feature, r.Feature);
           
            Weight = GetWeight();
            Score = GetScore();
        }

        private double GetWeight()
        {
            var r = (FragmentFeatureNode)RNode;
            if (LNode is PrecursorFeatureNode)
                return SubScoreFactory.GetKLDivergence(r.FragmentIonClassBase, _ratio, _lcCorrelation, _imsCorrelation, r.GroupParameter); 
            var l = (FragmentFeatureNode) LNode;
            return SubScoreFactory.GetKLDivergence(l.FragmentIonClassBase, r.FragmentIonClassBase, _ratio, _lcCorrelation, _imsCorrelation, l.GroupParameter);
        }

        private double GetScore()
        { // contains the node score
            var score = RNode.Score;
            var r = (FragmentFeatureNode)RNode;
            if (LNode is PrecursorFeatureNode)
            {
                score += LNode.Score;
                score += SubScoreFactory.GetRatioScore(r.FragmentIonClassBase, _ratio, r.GroupParameter);
                score += SubScoreFactory.GetLCCorrelationScore(r.FragmentIonClassBase, _lcCorrelation, r.GroupParameter);
                score += SubScoreFactory.GetIMSCorrelationScore(r.FragmentIonClassBase, _imsCorrelation, r.GroupParameter);
                //Console.WriteLine("pre " + score);
            }
            else
            {
                var l = (FragmentFeatureNode) LNode;
                var rscore = SubScoreFactory.GetRatioScore(l.FragmentIonClassBase, r.FragmentIonClassBase, _ratio, r.GroupParameter);
                var lscore = SubScoreFactory.GetLCCorrelationScore(l.FragmentIonClassBase, r.FragmentIonClassBase, _lcCorrelation, r.GroupParameter);
                var iscore = SubScoreFactory.GetIMSCorrelationScore(l.FragmentIonClassBase, r.FragmentIonClassBase, _imsCorrelation, r.GroupParameter);
                score += rscore + lscore + iscore;
              //  Console.WriteLine(this + "\tfrac " + score + "\t" + rscore + "\t" + lscore + "\t" + iscore + "\t" + RNode.Score);
              //  Console.WriteLine("raw " + _lcCorrelation + "\t" + _imsCorrelation + "\t" + RNode.IsotopeCorrelation);
              //  Console.WriteLine(LNode.Feature);
              //  Console.WriteLine(RNode.Feature);
            }
            return score;
        }

        static public int GetRatioScore(double v1, double v2)
        {
            if (v1 <= 0)
            {
                if (v2 > 0) return -5;
                return -6;
            }
            if (v2 <= 0) return 5;
            var f = 1;
            if (v1 <= v2)
            {
                var tmp = v1;
                v1 = v2;
                v2 = tmp;
                f = -1;
            }
            var index = 0;
            for (; index < 5; index++)
            {
                v1 = v1*0.66;
                if (v1 < v2) break;
            }
            return index * f;
        }

        static public int[] GetAllRatioScores()
        {
            return new[] {-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 };
        }
    }
}
