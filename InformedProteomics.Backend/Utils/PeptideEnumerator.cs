using System;
using System.Collections;
using System.Collections.Generic;

namespace InformedProteomics.Backend.Utils
{
    public class PeptideEnumerator : IEnumerable<Sequence>
    {
        public IEnumerator<Sequence> GetEnumerator()
        {
            throw new NotImplementedException();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }
}
