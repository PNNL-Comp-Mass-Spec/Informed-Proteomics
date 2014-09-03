using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    interface IChromatogramExtractor
    {
        // MS
        Xic GetFullPrecursorIonExtractedIonChromatogram(double mz, Tolerance tolerance);
        Xic GetFullPrecursorIonExtractedIonChromatogram(double minMz, double maxMz);

        // MS/MS
        Xic GetFullProductExtractedIonChromatogram(double mz, Tolerance tolerance, double precursorIonMz);
        Xic GetFullProductExtractedIonChromatogram(double minMz, double maxMz, double precursorIonMz);
    }
}
