using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public class DataDependentAcquisitionRun
    {
        public DataDependentAcquisitionRun(IRawDataParser rawDataParser)
        {
            _rawDataParser = rawDataParser;

            _ms1List = new LinkedList<Spectrum>();
            // Parse all spectra
            foreach (var spec in _rawDataParser.GetAllSpectra())
            {
                if (spec.MsLevel == 1) _ms1List.AddLast(spec);  // MS1 spectrum
            }

        }

        /// <summary>
        /// Gets the extracted ion chromatogram of the specified m/z (using only MS1 spectra)
        /// </summary>
        /// <param name="mz">target m/z</param>
        /// <param name="tolerance">m/z tolerance</param>
        /// <returns>XIC as an array.</returns>
        public double[] GetExtractedIonChromatogram(double mz, Tolerance tolerance)
        {
            var xic = new double[_ms1List.Count];
            var index = 0;
            foreach (var spec in _ms1List)
            {
                var peak = spec.FindPeak(mz, tolerance);
                xic[index++] = peak != null ? peak.Intensity : 0;
            }
            return xic;
        }
        
        /// <summary>
        /// Gets MS/MS spectra whose isolation windows contain the most abundant peak of the precursorIon
        /// </summary>
        /// <param name="precursorIon"></param>
        /// <returns></returns>
        public IEnumerable<ProductSpectrum> GetMatchingMs2Spectra(Ion precursorIon)
        {
            throw new System.NotImplementedException();
        }

        private readonly Dictionary<int, Spectrum> _scanNumSpecMap;  // scan number -> spectrum
        private readonly LinkedList<Spectrum> _ms1List; // linked list of MS1 spectra
        private readonly IRawDataParser _rawDataParser;
        private readonly int[] _numSpectra;
    }
}
