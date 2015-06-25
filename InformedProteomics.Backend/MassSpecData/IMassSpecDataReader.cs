using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public interface IMassSpecDataReader
    {
        /// <summary>
        /// Gets all spectra
        /// </summary>
        /// <returns>all spectra</returns>
        IEnumerable<Spectrum> ReadAllSpectra();

        /// <summary>
        /// The number of spectra in the file.
        /// </summary>
        int NumSpectra { get; }

        /// <summary>
        /// Close the reader
        /// </summary>
        void Close();
    }
}
