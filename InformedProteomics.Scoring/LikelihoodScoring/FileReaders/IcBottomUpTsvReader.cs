using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Scoring.LikelihoodScoring.Data;

namespace InformedProteomics.Scoring.LikelihoodScoring.FileReaders
{
    public class IcBottomUpTsvReader : IDataFileReader
    {
        public IcBottomUpTsvReader(string fileName, LazyLcMsRun lcms, bool decoy)
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

            var peptides = tsvFile.GetData(BottomUpPeptideHeader);
            if (scans == null)
            {
                throw new FormatException();
            }

            var pepQValues = tsvFile.GetData(PepQValueHeader);
            var formulas = tsvFile.GetData(FormulaHeader);

            var peptideSet = new HashSet<string>();

            for (var i = 0; i < peptides.Count; i++)
            {
                if (Convert.ToDouble(pepQValues[i]) > PepQValueThreshold || peptideSet.Contains(peptides[i]))
                {
                    continue;
                }

                peptideSet.Add(peptides[i]);
                var scanNum = Convert.ToInt32(scans[i]);
                //                    var spectrum = lcms.GetSpectrum(scanNum);
                //                    var spec = spectrum as ProductSpectrum;
                //                    if (spec == null || spec.ActivationMethod != Act) continue;
                var precursorCharge = Convert.ToInt32(precursorCharges[i]);

                var formula = formulas?[i] != null ? formulas[i] : string.Empty;

                specMatches.Add(new SpectrumMatch(peptides[i], DataFileFormat.IcBottomUp, _lcms, scanNum, precursorCharge, _decoy, formula));
            }
            return specMatches;
        }

        private readonly string _fileName;
        private readonly LazyLcMsRun _lcms;
        private readonly bool _decoy;

        private const string BottomUpPeptideHeader = "Peptide";
        private const string PepQValueHeader = "PepQValue";
        private const string FormulaHeader = "Formula";
        private const double PepQValueThreshold = 0.01;
        private const string ScanHeader = "ScanNum";
        private const string PrecursorChargeHeader = "Charge";
    }
}
