using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;
using pwiz.ProteowizardWrapper;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// This class currently doesn't work because of an unresolved dll dependency.
    /// </summary>
    public sealed class ProteoWizardReader: IMassSpecDataReader, IDisposable
    {

        public ProteoWizardReader(string filePath)
        {
            _reader = new MsDataFileImpl(filePath);
        }

        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            throw new NotImplementedException();
        }

        public Spectrum ReadMassSpectrum(int scanIndex)
        {
            var pwizSpec = _reader.GetCentroidedSpectrum(scanIndex);

            var msLevel = _reader.GetMsLevel(scanIndex);
            if (msLevel > 1)
            {
                return new ProductSpectrum(pwizSpec.Mzs, pwizSpec.Intensities, scanIndex)
                {
                    ActivationMethod = ActivationMethod.Unknown,
                    IsolationWindow = null,
                    MsLevel = msLevel
                };
            }
            return new Spectrum(pwizSpec.Mzs, pwizSpec.Intensities, scanIndex);
        }

        public void Close()
        {
            throw new NotImplementedException();
        }

        private readonly MsDataFileImpl _reader;
        private readonly int _minLcScan;
        private readonly int _maxLcScan;
        private readonly Dictionary<int, int> _msLevel;

        public void Dispose()
        {
            if (_reader != null) _reader.Dispose();
        }
    }
}
