using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Utils
{
    public class NoiseFilter
    {
        public NoiseFilter(int signalToNoiseRatio)
        {
            SignalToNoiseRatio = signalToNoiseRatio;
        }

        public int SignalToNoiseRatio { get; private set; }

        //public Spectrum Filter(Spectrum spec)
        //{
            
        //}
    }
}
