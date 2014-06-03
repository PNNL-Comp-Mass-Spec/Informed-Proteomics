using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.FileReaders
{
    class DiaTsvReader: IDataFileReader
    {
        public DiaTsvReader(string fileName, LazyLcMsRun lcms, bool decoy)
        {
            _fileName = fileName;
            _lcms = lcms;
            _decoy = decoy;
        }

        public IList<SpectrumMatch> Read()
        {
            var specMatches = new List<SpectrumMatch>();
            var tsvFile = new TsvFileParser(_fileName);
            var precursorCharges = tsvFile.GetData(PrecursorChargeHeader);
            var scans = tsvFile.GetData(ScanHeader);

            var peptides = tsvFile.GetData(TopDownPeptideHeader);
            if (peptides != null)
            {
                var peptideSet = new HashSet<string>();
                const double filterThreshold = QValueThreshold;
                var filterValues = tsvFile.GetData(QValueHeader);

                var aset = new AminoAcidSet();

                for (int i = 0; i < peptides.Count; i++)
                {
                    if (Convert.ToDouble(filterValues[i]) > filterThreshold || peptideSet.Contains(peptides[i])) continue;
                    peptideSet.Add(peptides[i]);
                    var scanNum = Convert.ToInt32(scans[i]);
                    int precursorCharge = Convert.ToInt32(precursorCharges[i]);
                    specMatches.Add(new SpectrumMatch(new Sequence(peptides[i], aset), _lcms, scanNum, precursorCharge, _decoy));
                }
            }
            return specMatches;
        }

        private readonly string _fileName;
        private readonly LazyLcMsRun _lcms;
        private readonly bool _decoy;

        private const string QValueHeader = "QValue";
        private const string TopDownPeptideHeader = "Sequence";
        private const double QValueThreshold = 0.01;
        private const string ScanHeader = "Scan";
        private const string PrecursorChargeHeader = "Charge";
    }
}
