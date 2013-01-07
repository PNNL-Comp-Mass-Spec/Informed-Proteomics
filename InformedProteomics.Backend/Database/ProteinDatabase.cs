using System;
using System.Collections;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Database
{
    public class ProteinEnumerator : IEnumerable<Sequence>
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
