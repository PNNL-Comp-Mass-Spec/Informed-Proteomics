using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.IMSTraining
{
    class MgfParser
    {
        public List<MSMSSpectrum> Spectra { get; private set; }
        public MgfParser(string fileName)
        {
            Spectra = new List<MSMSSpectrum>();
            //TODO
        }
    }
}
