using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.Backend.Quantification
{
    public class Ms2Feature
    {
        public int ScanNum { get; set; }
        public int Charge { get; set; }
        public double Mass { get; set; }

        public string DataSetId { get; set; }
        public string Sequence { get; set; }
    }

    public class MatchedMs1Feature : Ms1Feature
    {
        public MatchedMs1Feature(LcMsRun run, Ms2Feature ms2Feature) : base(run)
        {
            _ms2Feature = ms2Feature;
        }

        private readonly Ms2Feature _ms2Feature;
        //public List<Ms2Feature> Ms2Features { get; private set; }
    }
}
