using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.IMS;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Backend.IMSScoring
{
    public class ImsScorer
    {
        //TODO training and scoirng should be consistant in terms of dealing with the missing peaks! Double check!!
        //TODO Precursor node score should be added later - not through the cutscore thing.
        private readonly ImsDataCached _imsData;
        private readonly Ion _precursorIon;
        private readonly SubScoreFactory _scoringParams;
        private FragmentFeatureGraph _graph;
        private Feature _previousPrecursorFeature; // only for speed-up
        public PrecursorFeatureNode PrecursorFeatureNode { get; private set; }
        public List<IonType> SupportingIonTypes;

        internal ImsScorer(ImsDataCached imsData, Ion precursorIon, SubScoreFactory scoringParams) // precursorComposition does not have protons; however, protons are incorporated when calculating mz
        {
            _imsData = imsData;
            _precursorIon = precursorIon;
            _scoringParams = scoringParams;
        }

        public ImsDataCached ImsData
        {
            get { return _imsData; }
        }
        
        private void UpdatePrecursorFeatureNode(Feature precursorFeature)
        {
            if (_previousPrecursorFeature != null && _previousPrecursorFeature == precursorFeature) return;
            var parameter = new GroupParameter(_precursorIon);
            PrecursorFeatureNode =
                new PrecursorFeatureNode(
                    IsotopomerFeatures.GetPrecursorIsotopomerFeatures(_imsData, _precursorIon, precursorFeature),
                    parameter, _scoringParams);
            _previousPrecursorFeature = precursorFeature;
        }

        public double GetPrecursorScore(Feature precursorFeature)
        {
            UpdatePrecursorFeatureNode(precursorFeature);
            return PrecursorFeatureNode.GetScore();
        }

        public double GetCutScore(char nTermAA, char cTermAA, Composition cutComposition, Feature precursorFeature)
        {
            UpdatePrecursorFeatureNode(precursorFeature);
            var parameter = new GroupParameter(cutComposition, nTermAA, cTermAA, _precursorIon);
            _graph = new FragmentFeatureGraph(_imsData, PrecursorFeatureNode, precursorFeature, _precursorIon,
                                                 cutComposition, parameter, _scoringParams);
            SupportingIonTypes = _graph.supportingIonTypes;
            return _graph.Score;
        }

        public double GetNodeScore() // for debug only
        {
            return _graph == null ? double.NaN : _graph.NodeScore;
        }

        public double GetRatioScore() // for debug only
        {
            return _graph == null ? double.NaN : _graph.RatioScore;
        }

    }
}
