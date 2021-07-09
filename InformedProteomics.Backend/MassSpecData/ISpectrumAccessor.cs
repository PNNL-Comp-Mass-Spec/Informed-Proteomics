using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Interface for objects that can supply random access to spectra
    /// </summary>
    public interface ISpectrumAccessor : ISpectrumExtractor, IDisposable
    {
        /// <summary>
        /// Get the native ID of the specified scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns>Native ID string</returns>
        string GetNativeId(int scanNum);

        /// <summary>
        /// Index of first LC scan in the dataset
        /// </summary>
        int MinLcScan { get; }

        /// <summary>
        /// Index of last LC scan in the dataset
        /// </summary>
        int MaxLcScan { get; }

        /// <summary>
        /// The number of spectra in the file.
        /// </summary>
        int NumSpectra { get; }

        /// <summary>
        /// Get the elution time of the specified scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns>Elution time</returns>
        double GetElutionTime(int scanNum);

        /// <summary>
        /// Get the ion mobility drift time of the specified scan number (if data is ion mobility)
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns>Drift time</returns>
        double GetDriftTime(int scanNum);

        /// <summary>
        /// Gets the MS level of the specified scan
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <returns>MS level</returns>
        int GetMsLevel(int scanNum);

        /// <summary>
        /// Gets the greatest scan number smaller than ms2ScanNum
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <param name="msLevel">MS level</param>
        /// <returns>Previous scan number at the specified level</returns>
        int GetPrevScanNum(int scanNum, int msLevel);

        /// <summary>
        /// Gets the smallest scan number larger than ms2ScanNum
        /// </summary>
        /// <param name="scanNum">scan number</param>
        /// <param name="msLevel">MS level</param>
        /// <returns>Next scan number at the specified level</returns>
        int GetNextScanNum(int scanNum, int msLevel);

        /// <summary>
        /// Gets the scan numbers of the specified msLevel
        /// </summary>
        /// <param name="msLevel">MS level</param>
        /// <returns>Scan numbers of the specified msLevel</returns>
        IList<int> GetScanNumbers(int msLevel);

        /// <summary>
        /// Read and return the isolation window for the specified scan number
        /// </summary>
        /// <param name="scanNum"></param>
        /// <returns>Isolation window</returns>
        IsolationWindow GetIsolationWindow(int scanNum);
    }
}
