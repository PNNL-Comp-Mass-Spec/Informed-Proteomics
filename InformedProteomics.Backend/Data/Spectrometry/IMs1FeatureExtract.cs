using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface IMs1FeatureExtract
    {
        string GetFeatureFile(string rawFilePath, double minMass = 3000, double maxMass = 50000);
    }
    
    public interface IMs1FeaturePredictor
    {
        bool Predict(ChargeLcScanCluster cluster);
        double PredictProbability(ChargeLcScanCluster cluster);
    }
}
