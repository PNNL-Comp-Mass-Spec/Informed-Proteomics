using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Scoring;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics
{
    public class ICTopDown
    {
        public IEnumerable<Sequence> Sequences { get; set; }
        public int MinCharge { get; set; }
        public int MaxCharge { get; set; }
        //public ScoredSpectra Spectra { get; set; }
        //public ScoredSpectra Scorer { get; set; }

        public void ICTopdown()
        {
            foreach(Sequence sequence in Sequences)
            {
                SequenceVariantSet sv = new SequenceVariantSet();
            }

        }
        
    }
}
