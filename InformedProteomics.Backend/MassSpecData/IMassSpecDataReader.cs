using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;
using PSI_Interface.CV;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Interface for objects that can supply mass spectrometry data, usually from a file
    /// </summary>
    public interface IMassSpecDataReader : ISpectrumExtractor, IDisposable
    {
        /// <summary>
        /// Gets all spectra
        /// </summary>
        /// <param name="includePeaks"></param>
        /// <returns>all spectra</returns>
        IEnumerable<Spectrum> ReadAllSpectra(bool includePeaks = true);

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
        /// The Native Format of the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </summary>
        CV.CVID NativeFormat { get; }

        /// <summary>
        /// Path to the file; is <see cref="string.Empty"/> if the reader is in-memory
        /// </summary>
        string FilePath { get; }

        /// <summary>
        /// SHA-1 Checksum of the original input file (raw, mzML, .d folder, etc.) - lower case, hex characters only (no dashes)
        /// </summary>
        string SrcFileChecksum { get; }

        /// <summary>
        /// Version of the immediate prior input file (raw, mzML, .d folder, etc.)
        /// </summary>
        string FileFormatVersion { get; }
    }
}
