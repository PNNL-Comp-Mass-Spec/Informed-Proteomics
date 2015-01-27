using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface IMs1FeaturePredictor
    {
        bool Predict(ChargeLcScanCluster cluster);
        double PredictProbability(ChargeLcScanCluster cluster);
    }
}
