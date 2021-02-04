using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace InformedProteomics.TopDown.Scoring
{
    public class Ms1FtParser
    {
        public Ms1FtParser(string ms1FtFileName)
        {
            Ms1FtFileName = ms1FtFileName;
            if (!Read(ms1FtFileName))
            {
                Console.WriteLine("Error while reading {0}!", ms1FtFileName);
            }
        }

        public IEnumerable<ProMexFeature> GetAllFeatures()
        {
            return _featureArr.Where(feature => feature != null);
        }

        public ProMexFeature GetFeature(int featureId)
        {
            return _featureArr[featureId - 1];
        }

        public string Ms1FtFileName { get; }

        private static readonly string[] RequiredColumns = { "FeatureID" , "MinScan", "MaxScan", "MinCharge", "MaxCharge", "MonoMass", "Abundance"};

        private ProMexFeature[] _featureArr;
        private bool Read(string ms1FtFileName)
        {
            var featureList = new List<ProMexFeature>();

            int[] requiredColNum = null;
            var maxFeatureId = 0;
            foreach (var line in File.ReadLines(ms1FtFileName))
            {
                if (requiredColNum == null)
                {
                    requiredColNum = new int[RequiredColumns.Length];
                    for (var j = 0; j < requiredColNum.Length; j++)
                    {
                        requiredColNum[j] = -1;
                    }

                    var header = line.Split('\t');
                    for (var i = 0; i < header.Length; i++)
                    {
                        for (var j = 0; j < RequiredColumns.Length; j++)
                        {
                            if (header[i].Equals(RequiredColumns[j], StringComparison.InvariantCultureIgnoreCase))
                            {
                                requiredColNum[j] = i;
                                break;
                            }
                        }
                    }

                    if (requiredColNum.Any(colNum => colNum < 0))
                    {
                        return false;
                    }
                }

                var data = line.Split('\t');
                var featureId = Convert.ToInt32(data[requiredColNum[0]]);
                if (featureId > maxFeatureId)
                {
                    maxFeatureId = featureId;
                }

                var minScan = Convert.ToInt32(data[requiredColNum[1]]);
                var maxScan = Convert.ToInt32(data[requiredColNum[2]]);
                var minCharge = Convert.ToInt32(data[requiredColNum[3]]);
                var maxCharge = Convert.ToInt32(data[requiredColNum[4]]);
                var monoMass = Convert.ToDouble(data[requiredColNum[5]]);
                var abundance = Convert.ToDouble(data[requiredColNum[6]]);

                featureList.Add(new ProMexFeature(featureId, minScan, maxScan, minCharge, maxCharge, monoMass, abundance));
            }

            _featureArr = new ProMexFeature[maxFeatureId];
            foreach (var feature in featureList)
            {
                _featureArr[feature.FeatureId - 1] = feature;
            }

            return true;
        }
    }

    public class ProMexFeature
    {
        public ProMexFeature(int featureId, int minScan, int maxScan, int minCharge, int maxCharge, double monoMass, double abundance)
        {
            FeatureId = featureId;
            MinScan = minScan;
            MaxScan = maxScan;
            MinCharge = minCharge;
            MaxCharge = maxCharge;
            MonoMass = monoMass;
            Abundance = abundance;
        }

        public int FeatureId { get; }
        public int MinScan { get; }
        public int MaxScan { get; }
        public int MinCharge { get; }
        public int MaxCharge { get; }
        public double MonoMass { get; }
        public double Abundance { get; }
    }
}
