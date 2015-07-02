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
        /// Returns the spectrum specified by the scan number.
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns></returns>
        Spectrum ReadMassSpectrum(int scanNum);

        /// <summary>
        /// The number of spectra in the file.
        /// </summary>
        int NumSpectra { get; }

        /// <summary>
        /// Close the reader
        /// </summary>
        void Close();

        /// <summary>
        /// Try to make the reader random access capable
        /// </summary>
        /// <returns>true if is random access capable, false if not</returns>
        bool TryMakeRandomAccessCapable();
    }
}
