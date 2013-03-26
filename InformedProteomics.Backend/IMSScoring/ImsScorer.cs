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
        private Tolerance _fragmentTolerance = new Tolerance(20, DataReader.ToleranceType.PPM);

        public ImsScorer(ImsDataCached imsData)
        {
            _imsData = imsData;
        }
        //TODO

        //public double GetCutScore(Feature4D precursorFeature, Composition cutComposition, List<IonType> ionTypes)
        //{
        //    foreach (var ionType in ionTypes)
        //    {
        //        double fragmentMz = ionType.GetMz(cutComposition);
        //        var fragmentFeature = _imsData.GetFragmentXIC4D(precursorFeature, fragmentMz, _fragmentTolerance);
        //    }

        //    return 0;
        //}
    }
}
