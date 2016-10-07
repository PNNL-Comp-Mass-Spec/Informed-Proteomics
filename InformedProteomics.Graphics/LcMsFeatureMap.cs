using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using InformedProteomics.Backend.MassFeature;
using InformedProteomics.Backend.MassSpecData;
using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Axes;
using OxyPlot.Series;

using LinearAxis = OxyPlot.Axes.LinearAxis;

namespace InformedProteomics.Graphics
{
    public class LcMsFeatureMap
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="run">LcMsRun (for determining minimum and maximum elution time</param>
        /// <param name="features">List of features</param>
        /// <param name="title">Plot title</param>
        /// <param name="minMass">Minimum mass</param>
        /// <param name="maxMass">Maximum mass</param>
        public LcMsFeatureMap(LcMsRun run, IEnumerable<LcMsFeature> features, string title, double minMass, double maxMass) :
            this(features, title,
            minMass, maxMass,
            Math.Max(run.GetElutionTime(run.MinLcScan) - 5, 0),
            run.GetElutionTime(run.MaxLcScan) + 5)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="features">List of features</param>
        /// <param name="title">Plot title</param>
        /// <param name="minMass">Minimum mass</param>
        /// <param name="maxMass">Maximum mass</param>
        /// <param name="minTime">Minimum time (minutes)</param>
        /// <param name="maxTime">Maximum time (minutes)</param>
        public LcMsFeatureMap(IEnumerable<LcMsFeature> features, string title, double minMass, double maxMass, double minTime, double maxTime)
        {
            _features = features;

            // Initialize x and y axes.
            var yAxis = new LinearAxis
            {
                Position = AxisPosition.Left,
                Title = "Monoisotopic Mass [Da]",
                StringFormat = "0.###",
                FontSize = 20,
                Minimum = minMass,
                Maximum = maxMass
            };

            var xAxis = new LinearAxis
            {
                Position = AxisPosition.Bottom,
                Title = "Elution Time [Minute]",
                StringFormat = "0.###",
                FontSize = 20,
                Minimum = minTime,
                Maximum = maxTime
            };

            // Initialize feature map.
            _featureMap = new PlotModel { Title = title, TitleFontSize = 30 };
            _featureMap.Axes.Add(xAxis);
            _featureMap.Axes.Add(yAxis);

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

        public LcMsFeatureMap(LcMsRun run, string ms1FtPath, double minMass, double maxMass)
            : this(run,
            LcMsFeatureAlignment.LoadProMexResult(0, ms1FtPath, run, minMass, maxMass),
            Path.GetFileNameWithoutExtension(ms1FtPath),
            minMass, maxMass)
        {
        }

        public void SaveImage(string imgPath, int nColors = 30)
        {
            //var maxAbundance = _features.Max(p => p.Abundance);

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
                OxyPlot.Wpf.PngExporter.Export(_featureMap, stream, 1200, 900, OxyColor.FromRgb(byte.MaxValue, byte.MaxValue, byte.MaxValue));
            }
        }

        public int GetTableIndex(double abundance, double[] breakTable)
        {
            if (abundance <= breakTable[0]) return 0;
            if (abundance > breakTable[breakTable.Length - 1]) return -1;
            for (var i = 1; i < breakTable.Length; i++)
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
            var abundanceArray = features.Select(f => f.Abundance).ToArray();
            Array.Sort(abundanceArray);
            var qIndex = (int)Math.Floor(abundanceArray.Length / 4.0);
            var q1 = abundanceArray[qIndex];
            var q3 = abundanceArray[qIndex * 3];
            var iqr = q3 - q1;
            var outlierR = q3 + 1.5 * iqr;
            var nLeft = (int)Math.Floor(.75 * (nColors - 1));
            var nRight = (int)Math.Floor(.25 * (nColors - 1));
            var table = new double[nLeft + nRight];
            var deltaL = abundanceArray[qIndex * 2] / nLeft;
            var deltaR = (outlierR - abundanceArray[qIndex * 2]) / nRight;

            var val = 0d;
            for (var i = 0; i < nLeft; i++)
            {
                val += deltaL;
                table[i] = val;
            }

            val = abundanceArray[qIndex * 2];
            for (var i = nLeft; i < nRight + nLeft; i++)
            {
                val += deltaR;
                table[i] = val;
            }
            return table;
        }

        private readonly PlotModel _featureMap;
        private readonly IEnumerable<LcMsFeature> _features;
    }
}
