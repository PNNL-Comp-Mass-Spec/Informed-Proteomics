using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Database
{
    public class AnnotationAndOffset
    {
        public AnnotationAndOffset(long offset, string sequence)
        {
            Offset = offset;
            Annotation = sequence;
        }

        public long Offset { get; private set; }
        public string Annotation { get; private set; }
    }
}
