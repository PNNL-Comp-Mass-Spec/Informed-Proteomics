using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.MassFeature;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;

namespace InformedProteomics.Graphics
{
    public class LcMsFeatureMap
    {
        public LcMsFeatureMap(IEnumerable<LcMsFeature> features)
        {
            _features = features;
            
            // Initialize x and y axes.
            _yAxis = new LinearAxis { Position = AxisPosition.Left, Title = "Monoisotopic Mass", StringFormat = "0.###" };
            _xAxis = new LinearAxis { Position = AxisPosition.Bottom, Title = "Elution Time", StringFormat = "0.###", };

            // Initialize feature map.
            _featureMap = new PlotModel { Title = "Feature Map" };
            _featureMap.Axes.Add(_xAxis);
            _featureMap.Axes.Add(_yAxis);
        }

        public LcMsFeatureMap(string ms1FtPath) : this(LcMsFeatureAlignment.LoadProMexResult(ms1FtPath))
        {

        }

        public void SaveImage(string imgPath)
        {
            
               
        }
           
        private readonly LinearAxis _xAxis;
        private readonly LinearAxis _yAxis;
        private PlotModel _featureMap;
        private readonly IEnumerable<LcMsFeature> _features;
    }
}
