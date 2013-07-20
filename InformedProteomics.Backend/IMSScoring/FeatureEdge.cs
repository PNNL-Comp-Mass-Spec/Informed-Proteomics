using System;

namespace InformedProteomics.Backend.IMSScoring
{
    public class FeatureEdge
    {
        public FeatureNode LNode { get; private set; }
        public FeatureNode RNode { get; private set; }
        public double Weight { get; private set; } // used to calculate weight of a path
        private double _ratioScore;
        private double _lcScore;
        private double _imsScore;
        private double _nodeScore;

        private bool _isScoreCalculated; // for speed-up
        private readonly int _ratio;
        private readonly double _lcCorrelation;
        private readonly double _imsCorrelation;

        private readonly SubScoreFactory _scoringParams;

        public FeatureEdge(FeatureNode l, FeatureNode r, PrecursorFeatureNode p, SubScoreFactory scoringParams)
        {
            LNode = l;
            RNode = r;
            var ri = r.Feature == null ? 0.0 : r.Feature.IntensityMax;
            var li = l.Feature == null ? 0.0 : l.Feature.IntensityMax;
            _scoringParams = scoringParams;
            
            if (LNode is PrecursorFeatureNode) // TODO fix later when no summed intensity is used
                _ratio = GetRatioIndex(ri, ri);
            else
                _ratio = GetRatioIndex(li, ri);
            _lcCorrelation = StatisticsTools.GetLcCorrelation(p.Feature, r.Feature);
            _imsCorrelation = StatisticsTools.GetImsCorrelation(p.Feature, r.Feature);
            
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

        private void CalculateScores()
        { // contains the node score
            if (_isScoreCalculated) return;
            _nodeScore = RNode.GetScore();
            var r = (FragmentFeatureNode)RNode;
            if (LNode is PrecursorFeatureNode)
            {
                _ratioScore = _scoringParams.GetRatioScore(r.FragmentIonClassBase, _ratio, r.GroupParameter);
                var rr = 0.0;
                if(LNode.Feature != null)
                    rr = LNode.Feature.IntensityMax;

                if (RNode.Feature != null) rr /= RNode.Feature.IntensityMax;
                else rr = 0.0;

                if (rr > 100 || rr < .01) _ratioScore -= 4; //TODO 
                //Console.WriteLine("Prec : " + r.FragmentIonClassBase.Name +"\t" + _ratio + "\t" + RatioScore + "\t" + r.Feature);

                //Console.WriteLine("pre " + _ratioScore);
            }
            else
            {
                var l = (FragmentFeatureNode) LNode;
                _ratioScore = _scoringParams.GetRatioScore(l.FragmentIonClassBase, r.FragmentIonClassBase, _ratio, r.GroupParameter);
               
            }
            if (_lcCorrelation >= 0)
            {
                _lcScore = _scoringParams.GetLcCorrelationScore(r.FragmentIonClassBase, _lcCorrelation, r.GroupParameter); //TODO
                _imsScore = _scoringParams.GetImsCorrelationScore(r.FragmentIonClassBase, _imsCorrelation, r.GroupParameter);
            }
            else _lcScore = _imsScore = 0;

            _isScoreCalculated = true;
        }

        public double GetNodeScore()
        {
            if (!_isScoreCalculated) CalculateScores();
            return _nodeScore;
        }

        public double GetRatioScore()
        {
            if (!_isScoreCalculated) CalculateScores(); 
            return _ratioScore;
        }

        public double GetLcScore()
        {
            if (!_isScoreCalculated) CalculateScores();
            return _lcScore;
        }
        
        public double GetImsScore()
        {
            if (!_isScoreCalculated) CalculateScores();
            return _imsScore;
        }

        static public int GetRatioIndex(double v1, double v2)
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

        static public int[] GetAllRatioIndices()
        {
            return new[] {-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 };
        }
    }
}
