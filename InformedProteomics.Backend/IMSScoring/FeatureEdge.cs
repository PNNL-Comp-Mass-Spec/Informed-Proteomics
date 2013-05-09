using System;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FeatureEdge
    {
        public FeatureNode LNode { get; private set; }
        public FeatureNode RNode { get; private set; }
        public double Weight { get; private set; } // used to calculate weight of a path
        public double RatioScore { get; private set; } // used to calculate score
        public double LcScore { get; private set; }
        public double ImsScore { get; private set; }
        public double NodeScore { get; private set; }

        private readonly int _ratio;
        private readonly double _lcCorrelation;
        private readonly double _imsCorrelation;

        private readonly SubScoreFactory _scoringParams;

        public FeatureEdge(FeatureNode l, FeatureNode r, SubScoreFactory scoringParams)
        {
            LNode = l;
            RNode = r;
            var ri = r.Feature == null ? 0.0 : r.Feature.IntensityMax;
            var li = l.Feature == null ? 0.0 : l.Feature.IntensityMax;
            _scoringParams = scoringParams;
            
            if (LNode is PrecursorFeatureNode) // TODO fix later when no summed intensity is used
                _ratio = GetRatioScore(ri, ri);
            else
                _ratio = GetRatioScore(li, ri);
            _lcCorrelation = StatisticsTools.GetLcCorrelation(l.Feature, r.Feature);
            _imsCorrelation = StatisticsTools.GetImsCorrelation(l.Feature, r.Feature);
           
            GetScore(); // TODO calculate when Score is needed..
            Weight = GetWeight();
        }

        private double GetWeight()
        {
            var r = (FragmentFeatureNode)RNode;
            if (LNode is PrecursorFeatureNode)
                return _scoringParams.GetKLDivergence(r.FragmentIonClassBase, _ratio, _lcCorrelation, _imsCorrelation, r.GroupParameter); 
            var l = (FragmentFeatureNode) LNode;
            return _scoringParams.GetKLDivergence(l.FragmentIonClassBase, r.FragmentIonClassBase, _ratio, _lcCorrelation, _imsCorrelation, l.GroupParameter);
        }

        private void GetScore()
        { // contains the node score
            NodeScore = RNode.Score;
            var r = (FragmentFeatureNode)RNode;
            if (LNode is PrecursorFeatureNode)
            {
                RatioScore = _scoringParams.GetRatioScore(r.FragmentIonClassBase, _ratio, r.GroupParameter);
                //Console.WriteLine("Prec : " + r.FragmentIonClassBase.Name +"\t" + _ratio + "\t" + RatioScore + "\t" + r.Feature);
                if (_lcCorrelation>=0)
                {
                    LcScore = _scoringParams.GetLcCorrelationScore(r.FragmentIonClassBase, _lcCorrelation, r.GroupParameter);
                    ImsScore = _scoringParams.GetImsCorrelationScore(r.FragmentIonClassBase, _imsCorrelation, r.GroupParameter);
                }
                else LcScore = ImsScore = 0;
                //Console.WriteLine("pre " + score);
            }
            else
            {
                var l = (FragmentFeatureNode) LNode;
                RatioScore = _scoringParams.GetRatioScore(l.FragmentIonClassBase, r.FragmentIonClassBase, _ratio, r.GroupParameter);
                if (_lcCorrelation >= 0)
                {
                    LcScore = _scoringParams.GetLcCorrelationScore(l.FragmentIonClassBase, r.FragmentIonClassBase, _lcCorrelation, r.GroupParameter);
                    ImsScore = _scoringParams.GetImsCorrelationScore(l.FragmentIonClassBase, r.FragmentIonClassBase, _imsCorrelation, r.GroupParameter);
                }
                else LcScore = ImsScore = 0;
                //  Console.WriteLine(this + "\tfrac " + score + "\t" + rscore + "\t" + lscore + "\t" + iscore + "\t" + RNode.Score);
              //  Console.WriteLine("raw " + _lcCorrelation + "\t" + _imsCorrelation + "\t" + RNode.IsotopeCorrelation);
              //  Console.WriteLine(LNode.Feature);
              //  Console.WriteLine(RNode.Feature);
            }
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
