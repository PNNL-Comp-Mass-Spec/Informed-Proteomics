using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.MassSpecData;
using PSI_Interface.CV;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class LcmsRunDeconvoluter : IMassSpecDataReader
    {
        private readonly IMassSpecDataReader dataReader;

        private readonly Deconvoluter deconvoluter;

        private readonly double corrScoreThreshold;

        private readonly List<int> msLevels;

        public LcmsRunDeconvoluter(
                    IMassSpecDataReader dataReader,
                    Deconvoluter deconvoluter,
                    double corrScoreThreshold = 0.7)
        {
            this.deconvoluter = deconvoluter;
            this.dataReader = dataReader;
            this.corrScoreThreshold = corrScoreThreshold;
        }

        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            var spectra = this.dataReader.ReadAllSpectra();
            return this.dataReader.ReadAllSpectra()
                    .Where(spec => spec.MsLevel == 2)
                    .AsParallel()
                    .AsOrdered()
                    .WithDegreeOfParallelism(6)
                    .Select(spec => this.deconvoluter.GetCombinedDeconvolutedSpectrum(spec));
        }

        public Spectrum ReadMassSpectrum(int scanNum, bool includePeaks = true)
        {
            var spectrum = this.dataReader.ReadMassSpectrum(scanNum, includePeaks);
            if (includePeaks)
            {
                return this.deconvoluter.GetCombinedDeconvolutedSpectrum(spectrum);
            }

            return spectrum;
        }

        public int NumSpectra => this.dataReader.NumSpectra;

        public void Close()
        {
        }

        public bool TryMakeRandomAccessCapable()
        {
            return true;
        }

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000767, native spectrum identifier format
        /// </summary>
        public CV.CVID NativeIdFormat => this.dataReader.NativeIdFormat;

        /// <summary>
        /// The NativeIdFormat stored/used by the source file - needed for tracking purposes.
        /// Child term of PSI-MS term MS:1000560, mass spectrometer file format
        /// </summary>
        public CV.CVID NativeFormat => this.dataReader.NativeFormat;
    }
}
