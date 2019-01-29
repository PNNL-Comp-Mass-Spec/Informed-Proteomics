using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.ProbabilityTables
{
    public class Histogram <T1>
    {
        public List<List<T1>> Bins { get; private set; }
        public int Total { get; private set; }

        /// <summary>
        /// Histogram constructor for initializing the histogram with data on creation.
        /// </summary>
        /// <param name="data">Data to insert into histogram when it is created.</param>
        /// <param name="binEdges">Array of minimum values for each histogram bin.</param>
        /// <param name="compare">Comparer to compare items in histogram.</param>
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

        /// <summary>
        /// Histogram constructor for initializing an empty histogram with
        /// predefined bin edges.
        /// </summary>
        /// <param name="binEdges">Array of minimum values for each histogram bin.</param>
        /// <param name="compare">Comparer to compare items in histogram.</param>
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

        /// <summary>
        /// Histogram constructor for initializing an empty histogram with only one
        /// bin.
        /// </summary>
        /// <param name="compare">Comparer to compare items in histogram.</param>
        public Histogram(IComparer<T1> compare = null)
        {
            _compare = compare;
            if (compare == null)
                _compare = Comparer<T1>.Default;
            Total = 0;
            Bins = new List<List<T1>> {new List<T1>()};
        }

        /// <summary>
        /// Edges of histogram bins. When it is set, data already in
        /// histogram is rearranged into new bins.
        /// </summary>
        public T1[] BinEdges
        {
            get { return _binEdges; }
            set
            {
                _binEdges = value;
                var dataCollection = new List<T1>();
                foreach (var bin in Bins)
                {
                    dataCollection.AddRange(bin);
                }
                Bins = new List<List<T1>>();
                Compute(dataCollection);
            }
        }

        /// <summary>
        /// Array of histogram probabilities. Each item corresponds to
        /// a bin in the histogram. Found = number of items in the bin.
        /// Total = number of items in the histogram.
        /// </summary>
        public Probability<T1>[] Frequencies
        {
            get
            {
                return BinEdges.Select((t, i) => new Probability<T1>(t, Bins[i].Count, Total)).ToArray();
            }
        }

        /// <summary>
        /// Get the index of the bin that a given item will be sorted into by the histogram.
        /// </summary>
        /// <param name="datum">Item to find the index of.</param>
        /// <returns>Index of the bin the item will be sorted into.</returns>
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

        /// <summary>
        /// Organize histogram into bins where each bin has the same number of items.
        /// Last bin contains the remainder.
        /// </summary>
        /// <param name="bins">Number of equal bins to generate.</param>
        public void Equalize(int bins)
        {
            var dataCollection = new List<T1>();
            foreach (var bin in Bins)
            {
                dataCollection.AddRange(bin);
            }
            Equalize(bins, dataCollection);
        }

        /// <summary>
        /// Organize histogram into bins where each bin has the same number of items.
        /// Last bin contains the remainder.
        /// </summary>
        /// <param name="bins">Number of equal bins to generate.</param>
        /// <param name="minimum">The lowest possible value of the lowest bin.</param>
        public void Equalize(int bins, T1 minimum)
        {
            var dataCollection = new List<T1>();
            foreach (var bin in Bins)
            {
                dataCollection.AddRange(bin);
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

        /// <summary>
        /// Add a single item into the histogram.
        /// </summary>
        /// <param name="datum">Item to add to the histogram.</param>
        public void AddDatum(T1 datum)
        {
            if (BinEdges != null)
                Compute(new List<T1> { datum });
            else
            {
                Bins[0].Add(datum);
                Total++;
            }
        }

        /// <summary>
        /// Add a list of data to the histogram.
        /// </summary>
        /// <param name="data">List of data to add to the histogram.</param>
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

        private void Equalize(int bins, List<T1> dataCollection)
        {
            if (dataCollection.Count == 0) return;
            dataCollection.Sort(_compare);

            Bins = new List<List<T1>>();
            BinEdges = new T1[bins];
            for (int j = 0; j < bins; j++)
            var binSize = dataCollection.Count / bins;
            {
                int min = j * binSize;
                Bins[j].AddRange(dataCollection.GetRange(min, binSize));
                BinEdges[j] = dataCollection[min];

                if (j == bins - 1)
                    Bins[j].AddRange(dataCollection.GetRange(min + binSize, dataCollection.Count - (min + binSize)));
            }
            Total = dataCollection.Count;
        }

        private void Compute(IEnumerable<T1> data)
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

        private readonly IComparer<T1> _compare;
        private T1[] _binEdges;
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
