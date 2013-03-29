using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.IMS;
using UIMFLibrary;

namespace InformedProteomics.Backend.IMSScoring
{
    public class ImsScorer
    {
        private ImsDataCached _imsData;
        private Composition _precursorComposition;

        public ImsScorer(ImsDataCached imsData, Composition precursorComposition)
        {
            _imsData = imsData;
            _precursorComposition = precursorComposition;
        }

        public double GetCutScore(char nTermAA, char cTermAA, Composition cutComposition, Feature precursorFeature)
        {
            throw new NotImplementedException();
        }
    }
}
