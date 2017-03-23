using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using InformedProteomics.Backend.Data.Composition;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public interface IDeconvolutedSpectrum
    {
        DeconvolutedPeak FindPeak(Composition.Composition composition, Tolerance tolerance);
    }
}
