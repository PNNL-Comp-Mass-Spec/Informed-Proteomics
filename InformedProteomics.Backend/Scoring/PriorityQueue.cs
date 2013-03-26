using System;
using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Backend.Scoring
{
    public class PriorityQueue<TValue, TPriority> where TPriority : IComparable
    {
        private readonly SortedDictionary<TPriority, Queue<TValue>> _dict =
            new SortedDictionary<TPriority, Queue<TValue>>();

        public int Count { get; private set; }

        public bool Empty
        {
            get { return Count == 0; }
        }

        public void Enqueue(TValue val)
        {
            Enqueue(val, default(TPriority));
        }

        public void Enqueue(TValue val, TPriority pri)
        {
            ++Count;
            if (!_dict.ContainsKey(pri)) _dict[pri] = new Queue<TValue>();
            _dict[pri].Enqueue(val);
        }

        public TValue Dequeue()
        {
            --Count;
            var item = _dict.Last();
            if (item.Value.Count == 1) _dict.Remove(item.Key);
            return item.Value.Dequeue();
        }
    }

}
