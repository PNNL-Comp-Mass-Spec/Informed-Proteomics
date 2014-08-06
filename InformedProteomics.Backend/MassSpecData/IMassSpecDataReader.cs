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
        /// Close the reader
        /// </summary>
        void Close();
    }
}
