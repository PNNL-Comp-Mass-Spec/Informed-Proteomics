using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Security.Cryptography.X509Certificates;
using InformedProteomics.Backend.MassFeature;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;
//using OxyPlot.Wpf;

using LinearAxis = OxyPlot.Axes.LinearAxis;

namespace InformedProteomics.Graphics
{
    public class LcMsFeatureMap
    {
        public LcMsFeatureMap(IEnumerable<LcMsFeature> features, string title, double minTime, double maxTime, double minMass, double maxMass)
        {
            _features = features;
            // Initialize x and y axes.
            //var minMass = _features.Min(f => f.Mass);
            //var maxMass = _features.Max(f => f.Mass);
            _yAxis = new LinearAxis { Position = AxisPosition.Left, Title = "Monoisotopic Mass [Da]", StringFormat = "0.###", FontSize = 20};
            _xAxis = new LinearAxis { Position = AxisPosition.Bottom, Title = "Elution Time [Minute]", StringFormat = "0.###", FontSize = 20};
            _xAxis.Minimum = minTime; 
            _xAxis.Maximum = maxTime;
            _yAxis.Minimum = minMass;
            _yAxis.Maximum = maxMass;

            // Initialize feature map.
            _featureMap = new PlotModel { Title = title, TitleFontSize = 30 };
            _featureMap.Axes.Add(_xAxis);
            _featureMap.Axes.Add(_yAxis);

            var txtX = minTime + (maxTime - minTime) * 0.2;
            var txtY = maxMass - (maxMass - minMass) * 0.1;

            var annotation = new TextAnnotation
            {
                TextPosition = new DataPoint(txtX, txtY),
                Text = string.Format("Number of LCMS features = {0}", _features.Count()),
                FontSize = 25,
            };
            
            _featureMap.Annotations.Add(annotation);
        }

        public LcMsFeatureMap(string ms1FtPath, double minTime, double maxTime, double minMass, double maxMass)
            : this(LcMsFeatureAlignment.LoadProMexResult(ms1FtPath), Path.GetFileNameWithoutExtension(ms1FtPath), minTime, maxTime, minMass, maxMass)
        {

        }

        public void SaveImage(string imgPath, int nColors = 30)
        {
            var maxAbundance = _features.Max(p => p.Abundance);
            var colorBreakValues = GetColorBreaksTable(nColors, _features);
            foreach (var x in _features)
            {
                var line = new LineSeries();
                line.Points.Add(new DataPoint(x.MinElutionTime, x.Mass));
                line.Points.Add(new DataPoint(x.MaxElutionTime, x.Mass));
                var tableIndex = GetTableIndex(x.Abundance, colorBreakValues);
                line.Color = GetColor(tableIndex, colorBreakValues.Length);
                line.StrokeThickness = 0.3;
                _featureMap.Series.Add(line);
            }
            using (var stream = File.Create(imgPath))
            {
                //var pngExporter = new PngExporter();
                //pngExporter.Export(_featureMap, stream);
                OxyPlot.Wpf.PngExporter.Export(_featureMap, stream, 1200, 900, OxyColor.FromRgb(byte.MaxValue, byte.MaxValue, byte.MaxValue));
            }
        }

        public int GetTableIndex(double abundance, double[] breakTable)
        {

            if (abundance <= breakTable[0]) return 0;
            if (abundance > breakTable[breakTable.Count() - 1]) return -1;
            for (int i = 1; i < breakTable.Length; i++)
            {
                if (abundance > breakTable[i - 1] && abundance <= breakTable[i]) return i;
            }
            return breakTable.Length - 1;
        }

        public OxyColor GetColor(int value, int nColors)
        {

            if (value == -1) return OxyColor.FromHsv(0, 1, 1);
            var delta = 240.0 / (360.0 * (nColors));
            var hue = (240.0 / 360.0) - (delta * (value));
            return OxyColor.FromHsv(hue, 1, 1);
        }


        public double[] GetColorBreaksTable(int nColors, IEnumerable<LcMsFeature> features)
        {
            var abunanceArray = features.Select(f => f.Abundance).ToArray();
            Array.Sort(abunanceArray);
            var qIndex = (int)Math.Floor(abunanceArray.Count() / 4.0);
            var q1 = abunanceArray[qIndex];
            var q3 = abunanceArray[qIndex * 3];
            var iqr = q3 - q1;
            var outlierR = q3 + 1.5 * iqr;
            var nLeft = (int)Math.Floor(.75 * (nColors - 1));
            var nRight = (int)Math.Floor(.25 * (nColors - 1));
            var table = new double[nLeft + nRight];
            var deltaL = abunanceArray[qIndex * 2] / nLeft;
            var deltaR = (outlierR - abunanceArray[qIndex * 2]) / nRight;

            var val = 0d;
            for (int i = 0; i < nLeft; i++)
            {
                val += deltaL;
                table[i] = val;
            }
            val = abunanceArray[qIndex * 2];
            for (int i = nLeft; i < nRight + nLeft; i++)
            {
                val += deltaR;
                table[i] = val;
            }
            return table;
        }

        /*
        public void SaveImage(string imgPath)
        {
            var maxAbundance = _features.Max(p => p.Abundance);
            var colorBreakValues = GetColorBreaksTable(_features, maxAbundance);
            foreach (var x in _features)
            {
                var line = new LineSeries()
                {
                    StrokeThickness = 0.3,
                    MarkerSize = 0.3,
                };
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
                //var pngExporter = new PngExporter();
                OxyPlot.Wpf.PngExporter.Export(_featureMap, stream, 1200, 900, OxyColor.FromRgb(byte.MaxValue, byte.MaxValue, byte.MaxValue));
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
        */
        private readonly LinearAxis _xAxis;
        private readonly LinearAxis _yAxis;
        private PlotModel _featureMap;
        private readonly IEnumerable<LcMsFeature> _features;
    }
}
