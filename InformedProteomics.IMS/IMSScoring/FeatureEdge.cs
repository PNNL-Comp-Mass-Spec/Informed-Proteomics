namespace InformedProteomics.IMS.IMSScoring
{
    public class FeatureEdge
    {
        public FeatureNode LNode { get; private set; }
        public FeatureNode RNode { get; private set; }
        public double Weight { get; private set; } // used to calculate weight of a path
        private double _ratioScore;
        private double _nodeScore;

        private bool _isScoreCalculated; // for speed-up
        private readonly int _ratio;

        private readonly SubScoreFactory _scoringParams;

        public FeatureEdge(FeatureNode l, FeatureNode r, SubScoreFactory scoringParams)
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

            Weight = GetWeight();
        }

        private double GetWeight()
        {
            var r = (FragmentFeatureNode)RNode;
            if (LNode is PrecursorFeatureNode)
                return _scoringParams.GetKLDivergence(r.FragmentIonClassBase, _ratio, r.GroupParameter);
            var l = (FragmentFeatureNode) LNode;
            return _scoringParams.GetKLDivergence(l.FragmentIonClassBase, r.FragmentIonClassBase, _ratio, l.GroupParameter);
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

                //if (rr > 100 || rr < .01) _ratioScore -= 4; //TODO
                //Console.WriteLine("Prec : " + r.FragmentIonClassBase.Name +"\t" + _ratio + "\t" + RatioScore + "\t" + r.Feature);

                //Console.WriteLine("pre " + _ratioScore);
            }
            else
            {
                var l = (FragmentFeatureNode) LNode;
                _ratioScore = _scoringParams.GetRatioScore(l.FragmentIonClassBase, r.FragmentIonClassBase, _ratio, r.GroupParameter);
            }

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

        static public int GetRatioIndex(double v1, double v2)
        {
            //if (v2 <= 0) return 0;
            //if (v2 > v1) return -1;
            //return 1;
            if (v1 <= 0)
            {
                if (v2 > 0) return -4;
                return -5;
            }
            if (v2 <= 0) return 4;
            var f = 1;
            if (v1 <= v2)
            {
                var tmp = v1;
                v1 = v2;
                v2 = tmp;
                f = -1;
            }
            var index = 0;
            for (; index < 4; index++)
            {
                v1 = v1*0.66;
                if (v1 < v2) break;
            }
            return index * f;
            /*if (v1 <= 0)
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
            return index * f;*/
        }

        static public int[] GetAllRatioIndices()
        {
           // return new[] { -1, 0, 1 };
            return new[] { -5, -4, -3, -2, -1, 0, 1, 2, 3, 4};
           // return new[] {-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 };
        }
    }
}
