using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;
using PSI_Interface.CV;

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
        /// <param name="includePeaks"></param>
        /// <returns></returns>
        Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true);

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

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000767, native spectrum identifier format
        /// </summary>
        CV.CVID NativeIdFormat { get; }

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </summary>
        CV.CVID NativeFormat { get; }
    }
}
