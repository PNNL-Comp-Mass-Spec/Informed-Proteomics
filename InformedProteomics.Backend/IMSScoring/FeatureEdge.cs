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

            if (LNode is PrecursorFeatureNode) // TODO fix later when no summed intensity is used
                _ratio = GetRatioScore(r.Feature.IntensityMax, r.Feature.IntensityMax);
            else
                _ratio = GetRatioScore(l.Feature.IntensityMax, r.Feature.IntensityMax);
            _lcCorrelation = StatisticsTools.GetLCCorrelation(l.Feature, r.Feature);
            _imsCorrelation = StatisticsTools.GetIMSCorrelation(l.Feature, r.Feature);
           
            Weight = GetWeight();
            Score = GetScore();
        }

        private double GetWeight()
        {
            var r = (FragmentFeatureNode)RNode;
            if (LNode is PrecursorFeatureNode)
                return SubScoreFactory.GetKLDivergence(r.FragmentIonClassBase, _ratio, _lcCorrelation, _imsCorrelation, r.Parameter); 
            var l = (FragmentFeatureNode) LNode;
            return SubScoreFactory.GetKLDivergence(l.FragmentIonClassBase, r.FragmentIonClassBase, _ratio, _lcCorrelation, _imsCorrelation, l.Parameter);
        }

        private double GetScore()
        { // contains the node score
            var score = RNode.Score;
            var r = (FragmentFeatureNode)RNode;
            if (LNode is PrecursorFeatureNode)
            {
                score += LNode.Score;
                score += SubScoreFactory.GetRatioScore(r.FragmentIonClassBase, _ratio, r.Parameter);
                score += SubScoreFactory.GetLCCorrelationScore(r.FragmentIonClassBase, _lcCorrelation, r.Parameter);
                score += SubScoreFactory.GetIMSCorrelationScore(r.FragmentIonClassBase, _imsCorrelation, r.Parameter);
            }
            else
            {
                var l = (FragmentFeatureNode) LNode;
                score += SubScoreFactory.GetRatioScore(l.FragmentIonClassBase, r.FragmentIonClassBase, _ratio, r.Parameter);
                score += SubScoreFactory.GetLCCorrelationScore(l.FragmentIonClassBase, r.FragmentIonClassBase, _lcCorrelation, r.Parameter);
                score += SubScoreFactory.GetIMSCorrelationScore(l.FragmentIonClassBase, r.FragmentIonClassBase, _imsCorrelation, r.Parameter);
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
