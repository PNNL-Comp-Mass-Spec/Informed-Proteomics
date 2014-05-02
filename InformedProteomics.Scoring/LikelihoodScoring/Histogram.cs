using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class Histogram <T1>
    {
        private readonly IComparer<T1> _compare;
        private T1[] _binEdges;
        public List<List<T1>> Bins { get; private set; }

        public T1[] BinEdges
        {
            get { return _binEdges; }
            set
            {
                _binEdges = value;
                var datacollection = new List<T1>();
                foreach (var bin in Bins)
                {
                    datacollection.AddRange(bin);
                }
                Bins = new List<List<T1>>();
                Compute(datacollection);
            }
        }

        public int Total { get; private set; }

        public Histogram(List<T1> data, T1[] binEdges, IComparer<T1> compare=null)
        {
            _compare = compare;
            if (compare == null)
                _compare = Comparer<T1>.Default;
            Total = 0;
            Bins = new List<List<T1>>();
            for (int i = 0; i < binEdges.Length; i++)
            {
                Bins.Add(new List<T1>());
            }
            BinEdges = binEdges;
            Compute(data);
            Total += data.Count;
        }

        public Histogram(T1[] binEdges, IComparer<T1> compare = null)
        {
            _compare = compare;
            if (compare == null)
                _compare = Comparer<T1>.Default;
            Total = 0;
            Bins = new List<List<T1>>();
            for (int i = 0; i < binEdges.Length; i++)
            {
                Bins.Add(new List<T1>());
            }
            BinEdges = binEdges;
        }

        public Histogram(IComparer<T1> compare = null)
        {
            _compare = compare;
            if (compare == null)
                _compare = Comparer<T1>.Default;
            Total = 0;
            Bins = new List<List<T1>> {new List<T1>()};
        }

        public List<Probability<int>> Frequencies
        {
            get
            {
                return Bins.Select(bin => new Probability<int>(0, bin.Count, Total)).ToList();
            }
        }

        private void Equalize(int bins, List<T1> datacollection)
        {
            if (datacollection.Count == 0) return;
            datacollection.Sort(_compare);

            Bins = new List<List<T1>>();
            BinEdges = new T1[bins];
            for (int i = 0; i < bins; i++)
            {
                Bins.Add(new List<T1>());
            }
            int binSize = datacollection.Count / bins;
            for (int j = 0; j < bins; j++)
            {
                int min = j * binSize;
                Bins[j].AddRange(datacollection.GetRange(min, binSize));
                BinEdges[j] = datacollection[min];

                if (j == bins - 1)
                    Bins[j].AddRange(datacollection.GetRange(min + binSize, datacollection.Count - (min + binSize)));
            }
            Total = datacollection.Count;
        }

        public void Equalize(int bins)
        {
            var datacollection = new List<T1>();
            foreach (var bin in Bins)
            {
                datacollection.AddRange(bin);
            }
            Equalize(bins, datacollection);
        }

        public void Equalize(int bins, T1 minimum)
        {
            var datacollection = new List<T1>();
            foreach (var bin in Bins)
            {
                datacollection.AddRange(bin);
            }

            var excludeLow = new List<T1>();
            for (int i = 0; i < datacollection.Count; i++)
            {
                if (_compare.Compare(datacollection[i], minimum) >= 0)
                {
                    excludeLow.Add(datacollection[i]);
                }
            }
            Equalize(bins, excludeLow);
            if (bins > 0)
                BinEdges[0] = minimum;
        }

        public int GetBinIndex(T1 datum)
        {
            int i = 0;
            while (i < BinEdges.Length && (_compare.Compare(datum, BinEdges[i]) >= 0))
            {
                i++;
            }
            if (_compare.Compare(datum, BinEdges[BinEdges.Length - 1]) >= 0)
            {
                return BinEdges.Length - 1;
            }
            if (i != 0)
            {
                return i - 1;
            }
            return -1;
        }

        public void Compute(List<T1> data)
        {
            while (Bins.Count < BinEdges.Length)
            {
                Bins.Add(new List<T1>());
            }

            foreach (var datum in data)
            {
                int binIndex = GetBinIndex(datum);
                if (binIndex >= 0)
                {
                    Bins[binIndex].Add(datum);
                    Total++;
                }
            }
        }

        public void AddDatum(T1 datum)
        {
            if (BinEdges != null)
                Compute(new List<T1>{datum});
            else
            {
                Bins[0].Add(datum);
                Total++;
            }
        }

        public void AddData(List<T1> data)
        {
            if (BinEdges != null)
                Compute(data);
            else
            {
                Bins[0].AddRange(data);
                Total += data.Count;
            }
        }

        public int GetBinProbability(T1 value)
        {
            return GetBinIndex(value);
        }
    }

    public enum BinEdgeAlignment
    {
        Low,        // bin edge = low end of bin
        Center,     // bin edge = center of bin
        High        // bin edge = high end of bin
    };

    // bin edge alignment functions only for Histogram<double>, Histogram<float> and Histogram<int>
    static class NumericHistogramAlignment
    {
        public static double[] GetAlignedBinEdges(this Histogram<double> hist, BinEdgeAlignment align, double binWidth)
        {
            var binEdges = hist.BinEdges;
            var alignedEdges = new double[binEdges.Length];
            for (int i = 0; i < binEdges.Length; i++)
            {
                alignedEdges[i] = binEdges[i];
                if (align == BinEdgeAlignment.Center)
                {
                    var offset = binWidth/2;
                    alignedEdges[i] += offset;
                }
                else if (align == BinEdgeAlignment.High)
                {
                    alignedEdges[i] += binWidth;
                }
            }
            return alignedEdges;
        }
        public static float[] GetAlignedBinEdges(this Histogram<float> hist, BinEdgeAlignment align, float binWidth)
        {
            var binEdges = hist.BinEdges;
            var alignedEdges = new float[binEdges.Length];
            for (int i = 0; i < binEdges.Length; i++)
            {
                alignedEdges[i] = binEdges[i];
                if (align == BinEdgeAlignment.Center)
                {
                    var offset = binWidth / 2;
                    alignedEdges[i] += offset;
                }
                else if (align == BinEdgeAlignment.High)
                {
                    alignedEdges[i] += binWidth;
                }
            }
            return alignedEdges;
        }
        public static int[] GetAlignedBinEdges(this Histogram<int> hist, BinEdgeAlignment align, int binWidth)
        {
            var binEdges = hist.BinEdges;
            var alignedEdges = new int[binEdges.Length];
            for (int i = 0; i < binEdges.Length; i++)
            {
                alignedEdges[i] = binEdges[i];
                if (align == BinEdgeAlignment.Center)
                {
                    var offset = binWidth / 2;
                    alignedEdges[i] += offset;
                }
                else if (align == BinEdgeAlignment.High)
                {
                    alignedEdges[i] += binWidth;
                }
            }
            return alignedEdges;
        }
    }
}
