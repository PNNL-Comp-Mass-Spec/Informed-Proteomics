using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Xic
    {
        public Xic(XicPeak[] xicPeaks)
        {
            XicPeaks = xicPeaks;
        }

        public XicPeak[] XicPeaks { get; private set; }
    }
}
