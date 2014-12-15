using System.Collections.Generic;
using LibSVMsharp;
using LibSVMsharp.Extensions;
using System;
using InformedProteomics.Backend.Data.Spectrometry;

namespace MS1FeatureExtractor
{
    public class Ms1FeaturePredictor : IDisposable
    {
        public Ms1FeaturePredictor() : this(null) { }
        public Ms1FeaturePredictor(SVMModel model)
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
            var x = cluster.GetInputValueForSvm();
            var nodes = new SVMNode[x.Length];
            for (var i = 0; i < x.Length; i++)
                nodes[i] = new SVMNode(i + 1, x[i]);

            if (nodes.Predict(_ptrModel) > 0) return true;

            return false;
        }

        public double PredictProbability(ChargeLcScanCluster cluster)
        {
            var x = cluster.GetInputValueForSvm();
            var nodes = new SVMNode[x.Length];
            for (var i = 0; i < x.Length; i++)
                nodes[i] = new SVMNode(i + 1, x[i]);

            return PredictProbability(nodes);
        }

        public double PredictProbability(SVMNode[] nodes)
        {
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
}
