using System;
using System.Collections.Generic;
using System.Linq;
using System.Security;
using System.Text;
using MathNet.Numerics.Distributions;

namespace InformedProteomics.Test
{
    public class ChargeLcScanCell
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
    
    public class TestChargeLcScanMatrix
    {
        public int NumberOfLcScans { get; private set; }
        public int MinCharge { get; private set; }
        public int MaxCharge { get; private set; }
        
        public double[,] Data { get; private set; }
        public byte[,] ByteData { get; private set; }

        public double IntensityThreshold = 0.0;

        private List<List<ChargeLcScanCell>> _clusters;

        public int NRows
        {
            get { return MaxCharge - MinCharge + 1; }
        }

        public int NColumns
        {
            get { return NumberOfLcScans;  }
        }

        public TestChargeLcScanMatrix(int numberOfLcScans, int minCharge, int maxCharge)
        {
            NumberOfLcScans = numberOfLcScans;
            MinCharge = minCharge;
            MaxCharge = maxCharge;

            Data = new double[NRows, NColumns];
        }

        public void SetIntensityData(int charge, IEnumerable<double> intensities)
        {
            var c = 0;
            foreach (var i in intensities)
            {
                Data[charge - MinCharge, c] = i;
                ByteData[charge - MinCharge, c] = Convert.ToByte(i > IntensityThreshold);
                c++;
            }
        }

        public IEnumerable<List<ChargeLcScanCell>> GetSegment()
        {
            //return _segments.OrderByDescending(x => x.Count);
            return _clusters;
        }

        public void GetMostProbableScanNumbers(out int startScanNum, out int endScanNum)
        {
            var largestSegment = _clusters[0];
            startScanNum = NumberOfLcScans;
            endScanNum = 0;

            foreach (var cell in largestSegment)
            {
                if (cell.ScanIndex > endScanNum) endScanNum = cell.ScanIndex;
                if (cell.ScanIndex < startScanNum) startScanNum = cell.ScanIndex;
            }
        }
        
        public void Segmentation()
        {
	         //1) Label start as 1
	         //2) Get neighbor cells based on start point ->insert into a list L
	         //3) Label neighbor as 2

	        //std::queue<int_pair> neighbors; 
	        //int i = 0, j = 0;
            _clusters = new List<List<ChargeLcScanCell>>();

	        for (var i = 0; i < NRows; i++)
	        {
		        for (var j = 0; j < NColumns; j++)
		        {
		            if (ByteData[i, j] != 1) continue;
		            
		            var collectedCells  = new List<ChargeLcScanCell>();
		            var neighbors       = new Queue<ChargeLcScanCell>();

		            // pick a seed cell (i, j)
		            neighbors.Enqueue(new ChargeLcScanCell(i, j));

		            do 
		            {
		                var cell = neighbors.Dequeue();
		                collectedCells.Add(cell);

		                ByteData[cell.Row, cell.Col] = 2;  // mark it with a different bit value

		                //collect black neighbor cells around the focused pixel
		                for(var k = cell.Row - 1; k < cell.Row + 2; k++)
		                {
		                    for(var l = cell.Col - 1; l < cell.Col + 2; l++)
		                    {
		                        if (k < 0 || l < 0 || k >= NRows || l >= NColumns || ByteData[k, l] != 1) continue;
		                        neighbors.Enqueue(new ChargeLcScanCell(k, l));
		                        ByteData[k, l] = 2;
		                    }
		                }
		            } 
		            while(neighbors.Count > 0);

                    _clusters.Add(collectedCells);
		        }
	        }

            _clusters = _clusters.OrderByDescending(x => x.Count).ToList();

	        // restore
            for (var i = 0; i < NRows; i++)
                for (var j = 0; j < NColumns; j++)
                    if (ByteData[i, j] == 2) ByteData[i, j] = 1;
        }


    }
}
