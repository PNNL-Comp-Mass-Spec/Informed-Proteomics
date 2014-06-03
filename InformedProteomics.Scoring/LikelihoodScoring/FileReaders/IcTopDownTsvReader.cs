using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.FileReaders
{
    public class IcTopDownTsvReader: IDataFileReader
    {
        public IcTopDownTsvReader(string fileName, LazyLcMsRun lcms, bool decoy)
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
                var filterThreshold = QValueThreshold;
                var filterValues = tsvFile.GetData(QValueHeader);
                if (filterValues == null)
                {
                    filterValues = tsvFile.GetData(EvalueHeader);
                    filterThreshold = EValueThreshold;
                }

                var aset = new AminoAcidSet();

                for (int i = 0; i < peptides.Count; i++)
                {
                    if (Convert.ToDouble(filterValues[i]) > filterThreshold || peptideSet.Contains(peptides[i])) continue;
                    peptideSet.Add(peptides[i]);
                    var scanNum = Convert.ToInt32(scans[i]);
                    //                    var spectrum = lcms.GetSpectrum(scanNum);
                    //                    var spec = spectrum as ProductSpectrum;
                    //                    if (spec == null || spec.ActivationMethod != Act) continue;
                    int precursorCharge = Convert.ToInt32(precursorCharges[i]);
                    specMatches.Add(new SpectrumMatch(new Sequence(peptides[i], aset), _lcms, scanNum, precursorCharge, _decoy));
                }
            }
            return specMatches;
        }

        private readonly string _fileName;
        private readonly LazyLcMsRun _lcms;
        private readonly bool _decoy;

        private const string QValueHeader = "Qvalue";
        private const string EvalueHeader = "E-value";
        private const string TopDownPeptideHeader = "Annotation";
        private const double QValueThreshold = 0.01;
        private const double EValueThreshold = 0.1;
        private const string ScanHeader = "ScanNum";
        private const string PrecursorChargeHeader = "Charge";
    }
}
