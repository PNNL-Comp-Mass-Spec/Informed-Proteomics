using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.MassSpecData;
using PSI_Interface.CV;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Implementation of a deconvoluter that can be used to create a deconvoluted implementation of <see cref="LcMsRun"/>
    /// </summary>
    public class LcmsRunDeconvoluter : IMassSpecDataReader
    {
        /// <summary>
        /// MassSpec data reader to read raw spectra from.
        /// </summary>
        private readonly IMassSpecDataReader dataReader;

        /// <summary>
        /// Spectrum deconvoluter.
        /// </summary>
        private readonly Deconvoluter deconvoluter;

        /// <summary>
        /// The MS levels (ex MS1, MS2, etc) to extract spectra for.
        /// </summary>
        private readonly HashSet<int> msLevelSet;

        /// <summary>
        /// The maximum number of threads the deconvoluter can use.
        /// </summary>
        private readonly int maxDegreeOfParallelism;

        /// <summary>
        /// Initializes a new instance of the <see cref="LcmsRunDeconvoluter" /> class.
        /// This constructor creates an instance with multiple MSLevels for default (MS1 and MS2).
        /// </summary>
        /// <param name="dataReader">MassSpec data reader to read raw spectra from.</param>
        /// <param name="deconvoluter">Spectrum deconvoluter.</param>
        /// <param name="msLevels">The MS levels (ex MS1, MS2, etc) to extract spectra for.</param>
        /// <param name="maxDegreeOfParallelism">The maximum number of threads the deconvoluter can use.</param>
        public LcmsRunDeconvoluter(
                    IMassSpecDataReader dataReader,
                    Deconvoluter deconvoluter,
                    IEnumerable<int> msLevels = null,
                    int maxDegreeOfParallelism = 1)
        {
            this.deconvoluter = deconvoluter;
            this.dataReader = dataReader;
            this.maxDegreeOfParallelism = maxDegreeOfParallelism;
            this.msLevelSet = msLevels == null ? new HashSet<int> { 1, 2 } : new HashSet<int>(msLevels);
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="LcmsRunDeconvoluter" /> class.
        /// This constructor creates an instance with a single MSLevel
        /// </summary>
        /// <param name="dataReader">MassSpec data reader to read raw spectra from.</param>
        /// <param name="deconvoluter">Spectrum deconvoluter.</param>
        /// <param name="msLevel">The MS level (ex MS1, MS2, etc) to extract spectra for.</param>
        /// <param name="maxDegreeOfParallelism">The maximum number of threads the deconvoluter can use.</param>
        public LcmsRunDeconvoluter(
            IMassSpecDataReader dataReader,
            Deconvoluter deconvoluter,
            int msLevel,
            int maxDegreeOfParallelism = 1)
            : this(dataReader, deconvoluter, new HashSet<int> { msLevel }, maxDegreeOfParallelism)
        {
        }

        /// <summary>
        /// Gets all spectra.
        /// Deconvolutes spectra in parallel as it reads them.
        /// </summary>
        /// <returns>Deconvoluted spectra.</returns>
        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            return this.dataReader.ReadAllSpectra()
                       .Where(spec => this.msLevelSet.Contains(spec.MsLevel))
                       .AsParallel()
                       .AsOrdered()
                       .WithDegreeOfParallelism(this.maxDegreeOfParallelism)
                       .Select(spec => this.deconvoluter.GetCombinedDeconvolutedSpectrum(spec));
        }

        /// <summary>
        /// Returns the spectrum specified by the scan number and deconvolutes it.
        /// </summary>
        /// <param name="scanNum">The scan to deconvolute.</param>
        /// <param name="includePeaks">Should it be deconvoluted?</param>
        /// <returns>Deconvoluted spectrum.</returns>
        public Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true)
        {
            var spectrum = this.dataReader.ReadMassSpectrum(scanNum, includePeaks);
            if (includePeaks)
            {
                return this.deconvoluter.GetCombinedDeconvolutedSpectrum(spectrum);
            }

            return spectrum;
        }

        /// <summary>
        /// Gets the number of spectra in the file.
        /// </summary>
        public int NumSpectra => this.dataReader.NumSpectra;

        /// <summary>
        /// Close the reader.
        /// </summary>
        public void Close()
        {
            this.dataReader.Close();
        }

        /// <summary>
        /// Cleans up the reader
        /// </summary>
        public void Dispose()
        {
            this.dataReader.Dispose();
        }

        /// <summary>
        /// Try to make the reader random access capable
        /// </summary>
        /// <returns>true if is random access capable, false if not</returns>
        public bool TryMakeRandomAccessCapable()
        {
            return this.dataReader.TryMakeRandomAccessCapable();
        }

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000767, native spectrum identifier format
        /// </summary>
        public CV.CVID NativeIdFormat => this.dataReader.NativeIdFormat;

        /// <summary>
        /// The Native Format of the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </summary>
        public CV.CVID NativeFormat => this.dataReader.NativeFormat;

        /// <summary>
        /// Path to the file; is <see cref="string.Empty"/> if the reader is in-memory
        /// </summary>
        public string FilePath => this.dataReader.FilePath;

        /// <summary>
        /// SHA-1 Checksum of the original input file (raw, mzML, .d folder, etc.)
        /// </summary>
        public string SrcFileChecksum => this.dataReader is IPbfLcMsRun run ? run.PbfFileChecksum : this.dataReader.SrcFileChecksum;

        /// <summary>
        /// Version of the immediate prior input file (raw, mzML, .d folder, etc.)
        /// </summary>
        public string FileFormatVersion => this.dataReader.FileFormatVersion;
    }
}
