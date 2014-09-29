using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Security;
using System.Text;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Distributions;

namespace InformedProteomics.Backend.Data.Spectrometry
{

    public class ChargeLcScanCluster
    {
        public int MinCharge { get; private set; }
        public int MaxCharge { get; private set; }

        public int MinScanIndex { get; private set; }
        public int MaxScanIndex { get; private set; }

        public ChargeLcScanCluster(int minC, int maxC, int startScan, int endScan)
        {
            MinScanIndex = startScan;
            MaxScanIndex = endScan;
            MinCharge = minC;
            MaxCharge = maxC;
        }

    }

    class ChargeLcScanCell
    {
        public int Row;
        public int Col;

        public ChargeLcScanCell(int row, int col)
        {
            Row = row;
            Col = col;
        }

        public int ChargeIndex
        {
            get { return Row; }
        }

        public int ScanIndex
        {
            get { return Col;  }
        }
    }
    
    public class ChargeLcScanMatrix
    {
        public int NumberOfLcScans { get; private set; }
        public int MinCharge { get; private set; }
        public int MaxCharge { get; private set; }
        
        public double[,] Data { get; private set; }

        private List<double> _nonZeroIntensities;

        public double IntensityThreshold = 0.0;

        //private List<List<ChargeLcScanCell>> _clusters;
        //private List<KeyValuePair<double, List<ChargeLcScanCell>>> _clusters;

        public int NRows
        {
            get { return MaxCharge - MinCharge + 1; }
        }

        public int NColumns
        {
            get { return NumberOfLcScans;  }
        }

        public ChargeLcScanMatrix(int numberOfLcScans, int minCharge, int maxCharge)
        {
            NumberOfLcScans = numberOfLcScans;
            MinCharge = minCharge;
            MaxCharge = maxCharge;

            Data = new double[NRows, NColumns];
            _nonZeroIntensities = new List<double>();
        }

        public void SetIntensityData(int charge, IEnumerable<double> intensities)
        {
            var c = 0;
            
            foreach (var i in intensities)
            {
                Data[charge - MinCharge, c] = i;
                //ByteData[charge - MinCharge, c] = Convert.ToByte(i > IntensityThreshold);
                c++;
                if (i > 0) _nonZeroIntensities.Add(i);
            }

        }

        public IEnumerable<ChargeLcScanCluster> GetProbableChargeScanRegions()
        {
            double intensityThreshold = 0;

            if (_nonZeroIntensities.Count > 2)
                intensityThreshold = _nonZeroIntensities[_nonZeroIntensities.Count / 2];

            var clusters = Segmentation(intensityThreshold);

            for (var i = 0; i < clusters.Count; i++)
            {
                var largestSegment = clusters[i].Value;

                var startScanNum = NumberOfLcScans;
                var endScanNum = 0;

                var maxC = MinCharge;
                var minC = MaxCharge;

                foreach (var cell in largestSegment)
                {
                    if (cell.ScanIndex > endScanNum) endScanNum = cell.ScanIndex;
                    if (cell.ScanIndex < startScanNum) startScanNum = cell.ScanIndex;

                    if (cell.ChargeIndex + MinCharge > maxC) maxC = cell.ChargeIndex + MinCharge;
                    if (cell.ChargeIndex + MinCharge < minC) minC = cell.ChargeIndex + MinCharge;
                }

                if (i > 10 || (i > 1 && clusters[i].Key < 10)) yield break;

                yield  return new ChargeLcScanCluster(minC, maxC, startScanNum, endScanNum);
            }
            
        }

        private double ComputeScore(List<ChargeLcScanCell> cluster)
        {
            double score = cluster.Count;
            /*var startScan = NumberOfLcScans;
            var endScan = 0;

            foreach (var cell in cluster)
            {
                if (startScan > cell.Col) startScan = cell.Col;
                if (endScan < cell.Col) endScan = cell.Col;
            }
            */

            return score;
        }
        
        private List<KeyValuePair<double, List<ChargeLcScanCell>>> Segmentation(double intensityThreshold)
        {
            var byteData = new byte[NRows, NColumns];
            
            for(var i = 0; i < NRows; i++)
            {
                for (var j = 0; j < NColumns; j++)
                {
                    if (Data[i, j] >= intensityThreshold) byteData[i, j] = 1;

                }
            }

            //1) Label start as 1
	        //2) Get neighbor cells based on start point ->insert into a list L
	        //3) Label neighbor as 2
            var clusters = new List<KeyValuePair<double, List<ChargeLcScanCell>>>();

	        for (var i = 0; i < NRows; i++)
	        {
		        for (var j = 0; j < NColumns; j++)
		        {
		            if (byteData[i, j] != 1) continue;
		            
		            var collectedCells  = new List<ChargeLcScanCell>();
		            var neighbors       = new Queue<ChargeLcScanCell>();

		            // pick a seed cell (i, j)
		            neighbors.Enqueue(new ChargeLcScanCell(i, j));

		            do 
		            {
		                var cell = neighbors.Dequeue();
		                collectedCells.Add(cell);

		                byteData[cell.Row, cell.Col] = 2;  // mark it with a different bit value

		                //collect black neighbor cells around the focused pixel
		                for(var k = cell.Row - 1; k < cell.Row + 2; k++)
		                {
		                    for(var l = cell.Col - 1; l < cell.Col + 2; l++)
		                    {
		                        if (k < 0 || l < 0 || k >= NRows || l >= NColumns || byteData[k, l] != 1) continue;
		                        neighbors.Enqueue(new ChargeLcScanCell(k, l));
		                        byteData[k, l] = 2;
		                    }
		                }
		            } 
		            while(neighbors.Count > 0);

		            double clusterScore = ComputeScore(collectedCells);

                    clusters.Add(new KeyValuePair<double, List<ChargeLcScanCell>>(clusterScore, collectedCells));
		        }
	        }

            return clusters.OrderByDescending(x => x.Key).ToList();
        }

    }
}
