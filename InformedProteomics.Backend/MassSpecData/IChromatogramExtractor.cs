using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Interface for objects that can supply chromatograms
    /// </summary>
    public interface IChromatogramExtractor
    {
        // MS

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z (using only MS1 spectra)
        /// XicPoint is created for every MS1 scan.
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <returns>XIC as an Xic object</returns>
        Xic GetFullPrecursorIonExtractedIonChromatogram(double mz, Tolerance tolerance);

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS1 spectra)
        /// XicPoint is created for every MS1 scan.
        /// </summary>
        /// <param name="minMz">min m/z</param>
        /// <param name="maxMz">max m/z</param>
        /// <returns>XIC as an Xic object</returns>
        Xic GetFullPrecursorIonExtractedIonChromatogram(double minMz, double maxMz);

        // MS/MS

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z range (using only MS2 spectra)
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <param name="precursorIonMz">precursor m/z of the precursor ion</param>
        /// <returns>XIC as an Xic object</returns>
        Xic GetFullProductExtractedIonChromatogram(double mz, Tolerance tolerance, double precursorIonMz);

        /// <summary>
        /// Returns a xic for the chosen range that covers the entire run.
        /// </summary>
        /// <param name="minMz"></param>
        /// <param name="maxMz"></param>
        /// <param name="precursorIonMz"></param>
        /// <returns></returns>
        Xic GetFullProductExtractedIonChromatogram(double minMz, double maxMz, double precursorIonMz);
    }
}
