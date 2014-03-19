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

        public List<Probability> Frequencies
        {
            get
            {
                return Bins.Select(bin => new Probability(bin.Count, Total)).ToList();
            }
        }

        public void Equalize(int bins)
        {
            var datacollection = new List<T1>();
            foreach (var bin in Bins)
            {
                datacollection.AddRange(bin);
            }

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
                int min = j*binSize;
                Bins[j].AddRange(datacollection.GetRange(min, binSize));
                BinEdges[j] = datacollection[min];

                if (j == bins-1)
                    Bins[j].AddRange(datacollection.GetRange(min+binSize, datacollection.Count-(min+binSize)));
            }
            Total = datacollection.Count;
        }

        public void Compute(List<T1> data)
        {
            int total = 0;

            while (Bins.Count < BinEdges.Length)
            {
                Bins.Add(new List<T1>());
            }

            foreach (var datum in data)
            {
                int i = 0;
                while (i < BinEdges.Length && (_compare.Compare(datum, BinEdges[i]) >= 0))
                {
                    i++;
                }
                if (_compare.Compare(datum, BinEdges[BinEdges.Length - 1]) >= 0)
                {
                    Bins[BinEdges.Length - 1].Add(datum);
                    total++;
                }
                else if (i != 0)
                {
                    Bins[i - 1].Add(datum);
                    total++;
                }
            }
            Total += total;
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
    }
}
