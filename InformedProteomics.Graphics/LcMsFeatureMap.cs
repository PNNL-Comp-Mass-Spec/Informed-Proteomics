using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Security.Cryptography.X509Certificates;
using InformedProteomics.Backend.MassFeature;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.Wpf;
using LinearAxis = OxyPlot.Axes.LinearAxis;

namespace InformedProteomics.Graphics
{
    public class LcMsFeatureMap
    {
        public LcMsFeatureMap(IEnumerable<LcMsFeature> features, string title, double maxElutionTime, double minElutionTime = 0)
        {
            _features = features;

            // Initialize x and y axes.
            _yAxis = new LinearAxis { Position = AxisPosition.Left, Title = "Monoisotopic Mass", StringFormat = "0.###" };
            _xAxis = new LinearAxis { Position = AxisPosition.Bottom, Title = "Elution Time", StringFormat = "0.###", };
            _xAxis.Minimum = minElutionTime; //_features.Min(p => p.MinElutionTime) - 20;
            _xAxis.Maximum = maxElutionTime; //_features.Max(p => p.MaxElutionTime) + 10;

            // Initialize feature map.
            _featureMap = new PlotModel { Title = title };
            _featureMap.Axes.Add(_xAxis);
            _featureMap.Axes.Add(_yAxis);

        }

        public LcMsFeatureMap(string ms1FtPath, double maxElutionTime, double minElutionTime = 0)
            : this(LcMsFeatureAlignment.LoadProMexResult(ms1FtPath), Path.GetFileNameWithoutExtension(ms1FtPath), maxElutionTime, minElutionTime)
        {

        }

        public void SaveImage(string imgPath)
        {
            var maxAbundance = _features.Max(p => p.Abundance);
            var colorBreakValues = GetColorBreaksTable(_features, maxAbundance);
            foreach (var x in _features)
            {
                var line = new OxyPlot.Series.LineSeries();
                line.Points.Add(new DataPoint(x.MinElutionTime, x.Mass));
                line.Points.Add(new DataPoint(x.MaxElutionTime, x.Mass));
                var colorValue = x.Abundance / maxAbundance;
                colorValue = GetColorValue(colorValue, colorBreakValues);
                var r = GetRedValue(colorValue);
                var g = GetGreenValue(colorValue);
                var b = GetBlueValue(colorValue);
                line.Color = OxyColor.FromRgb(r, g, b);
                _featureMap.Series.Add(line);

            }
            using (var stream = File.Create(imgPath))
            {
                var pngExporter = new PngExporter();
                pngExporter.Export(_featureMap, stream);
            }

        }

        public byte GetBlueValue(double colorValue)
        {
            if (colorValue < .4) return Convert.ToByte(255);
            if (colorValue >= .4 && colorValue < .5) return Convert.ToByte(128);
            return Convert.ToByte(0);
        }

        public byte GetGreenValue(double colorValue)
        {
            if (colorValue < .2 || colorValue >= .9) return Convert.ToByte(0);
            if (colorValue >= .2 && colorValue < .3) return Convert.ToByte(128);
            if (colorValue >= .8 && colorValue < .9) return Convert.ToByte(128);
            return Convert.ToByte(255);
        }

        public byte GetRedValue(double colorValue)
        {
            if (colorValue < .6) return Convert.ToByte(0);
            if (colorValue >= .6 && colorValue < .7) return Convert.ToByte(128);
            return Convert.ToByte(255);
        }

        public double GetColorValue(double value, double[] breakTable)
        {
            for (int i = 0; i < breakTable.Length; i++)
            {
                if (value < breakTable[i])
                {
                    if (i == 0)
                    {
                        return 0;
                    }
                    if (value > breakTable[i - 1])
                    {
                        return i * .1;
                    }
                }
            }
            return 1;
        }

        public double[] GetColorBreaksTable(IEnumerable<LcMsFeature> features, double maxAbundance)
        {
            var colorBreaks = new double[9];
            var abundanceList = new List<double>();
            foreach (var x in features)
            {
                abundanceList.Add(x.Abundance);
            }
            abundanceList.Sort();
            var bins = abundanceList.Count;
            var percentage = .1;
            for (int i = 0; i < colorBreaks.Length; i++)
            {
                var bin = (int)Math.Floor(bins * percentage);
                colorBreaks[i] = abundanceList[bin] / maxAbundance;
                percentage += .1;
            }
            return colorBreaks;
        }

        private readonly LinearAxis _xAxis;
        private readonly LinearAxis _yAxis;
        private PlotModel _featureMap;
        private readonly IEnumerable<LcMsFeature> _features;
    }
}
