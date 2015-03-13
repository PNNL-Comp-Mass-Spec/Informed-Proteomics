using System.Collections.Generic;
//using LibSVMsharp;
//using LibSVMsharp.Extensions;
using System;
using InformedProteomics.Backend.Data.Spectrometry;

namespace ProMex
{
    /*
    public class Ms1FeatureSvmPredictor : IDisposable, IMs1FeaturePredictor
    {
        public Ms1FeatureSvmPredictor(string svmFilePath)
        {
            var svmModel = SVM.LoadModel(svmFilePath);
            if (svmModel == null) throw new Exception("Cannot load SVM Model");
            SetModel(svmModel);
        }

        public Ms1FeatureSvmPredictor(SVMModel model)
        {
            SetModel(model);
        }

        public SVMModel Model { get; private set; }
        public void SetModel(SVMModel model)
        {
            if (model != null)
            {
                Model = model.Clone();
                _ptrModel = SVMModel.Allocate(Model);
            }
        }

        public bool Predict(ChargeLcScanCluster cluster)
        {
            var nodes = new SVMNode[2 + ChargeLcScanScore.Count];
            nodes[0] = new SVMNode(1, cluster.EnvelopeCount);
            nodes[1] = new SVMNode(2, cluster.ScanLength);
            for (var i = 2; i < nodes.Length; i++)
                nodes[i] = new SVMNode(i + 1, cluster.GetScore(i - 2));

            if (nodes.Predict(_ptrModel) > 0) return true;

            return false;
        }

        public double PredictProbability(ChargeLcScanCluster cluster)
        {
            var nodes = new SVMNode[2 + ChargeLcScanScore.Count];
            nodes[0] = new SVMNode(1, cluster.EnvelopeCount);
            nodes[1] = new SVMNode(2, cluster.ScanLength);
            for (var i = 2; i < nodes.Length; i++)
                nodes[i] = new SVMNode(i + 1, cluster.GetScore(i-2));

            double[] estimations;
            var label = nodes.PredictProbability(_ptrModel, out estimations);

            return estimations[0];
        }


        private IntPtr _ptrModel = IntPtr.Zero;

        public double Predict(SVMNode[] x)
        {
            return x.Predict(_ptrModel);
        }
        public double PredictProbability(SVMNode[] x, out double[] estimations)
        {
            return x.PredictProbability(_ptrModel, out estimations);
        }
        public double PredictValues(SVMNode[] x, out double[] values)
        {
            return x.PredictValues(_ptrModel, out values);
        }
        public void Dispose()
        {
            SVMModel.Free(_ptrModel);
            _ptrModel = IntPtr.Zero;
            Model = null;
        }
    }

    public class Ms1FeatureLogisticRegPredictor : IMs1FeaturePredictor
    {
        public bool Predict(ChargeLcScanCluster cluster)
        {
            var prob = PredictProbability(cluster);

            return (prob > 0.5);
        }

        private static readonly double[] LogisticRegressionBetaHighMass = new double[]
        {
            -10.4321898371984,
            -0.00963637466315592,
            0.0157956024150966,
            0.00768770822894664,
            0.734055032398606,
            1.94789462501870,
            0.139890063816759,
            -0.0740586195021866,
            0.111790998170121,
            -0.0359792304409301,
            0.0274524869426300,
            -7.16014026630920,
            -32.5943322505819,
            -0.0428641979961823,
            2.76298440986434,
            1.70629209273113,
            13.8960083484043,
            -2.88365384325922
        };

        private static readonly double[] LogisticRegressionBetaLowMass = new double[]
        {
            -8.79560214946065,
            0.0286885357478392,
            0.0187266836377149,
            0.0348032419223919,
            2.55378705076427,
            -1.21797576452347,
            0.160182423693798,
            0.233556214398555,
            0.0168989700769169,
            0.0404199890240891,
            0.0311439422783994,
            -22.5120843419329,
            -8.15709551975472,
            0.0801906075471353,
            0.546285450031053,
            2.45833601052128,
            9.85177668060179,
            -0.732734888083536
        };

        public double PredictProbability(ChargeLcScanCluster cluster)
        {
            double[] logisticRegressionBetaVector = null;

            if (cluster.RepresentativeMass < 13000) logisticRegressionBetaVector = LogisticRegressionBetaLowMass;
            else logisticRegressionBetaVector = LogisticRegressionBetaHighMass;


            var eta = logisticRegressionBetaVector[0] +
                                 cluster.MinCharge * logisticRegressionBetaVector[1] +
                                 cluster.MaxCharge * logisticRegressionBetaVector[2] +
                                 cluster.EnvelopeCount * logisticRegressionBetaVector[3];

            for (var i = 4; i < logisticRegressionBetaVector.Length; i++)
                eta += cluster.GetScore(i - 4) * logisticRegressionBetaVector[i];

            var pi = Math.Exp(eta);
            return pi / (pi + 1);

            //return (cluster.RepresentativeMass < 13000) ? prob : prob * 10.0d;

        }

    }*/
}
