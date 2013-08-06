using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public interface IRawDataParser
    {
        /// <summary>
        /// Gets the mass spectrum with the specified scanNum
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>mass spectrum</returns>
        Spectrum GetMassSpectrum(int scanNum);

        /// <summary>
        /// Gets all spectra
        /// </summary>
        /// <returns>all spectra</returns>
        IEnumerable<Spectrum> GetAllSpectra();

        /// <summary>
        /// Gets the scan numbers of the specified msLevel
        /// </summary>
        /// <param name="msLevel">MS level</param>
        /// <returns>scan numbers of the specified msLevel</returns>
        IEnumerable<int> GetScanNumbers(int msLevel);

        /// <summary>
        /// Gets the precursor information of the specified scan
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>precursor information</returns>
        PrecursorInfo GetPrecursorInfo(int scanNum);

        /// <summary>
        /// Gets the total number of spectra (all levels)
        /// </summary>
        /// <returns>total number of spectra</returns>
        int GetNumSpectra();

        /// <summary>
        /// Gets the total number of spectra of the specified ms level
        /// </summary>
        /// <returns>total number of spectra</returns>
        int GetNumSpectra(int msLevel);

    }
}
