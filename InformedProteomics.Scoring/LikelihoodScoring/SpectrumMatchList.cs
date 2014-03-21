using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class SpectrumMatchList
    {
        private const string TopDownScanHeader = "Scan(s)";
        private const string FdrHeader = "FDR";
        private const string EvalueHeader = "E-value";
        private const double FdrThreshold = 0.01;
        private const double EValueThreshold = 0.1;

        private const string BottomUpScanHeader = "ScanNum";
        private const string PepQValueHeader = "PepQValue";
        private const double PepQValueThreshold = 0.01;
        private const string FormulaHeader = "Formula";

        private const string PrecursorChargeHeader = "Charge";
        private const string PeptideHeader = "Peptide";
        private const int NumMutations = 3;

        public List<SpectrumMatch> Matches;

        private string Trim(string prot)
        {
            int start = prot.IndexOf('.') + 1;
            int length = prot.LastIndexOf('.') - start;
            return prot.Substring(start, length);
        }

        public SpectrumMatchList(LcMsRun lcms, TsvFileParser tsvFile, ActivationMethod act, bool useDecoy=false)
        {
            Matches = new List<SpectrumMatch>();

            var peptides = tsvFile.GetData(PeptideHeader);
            var precursorCharges = tsvFile.GetData(PrecursorChargeHeader);
            var scans = tsvFile.GetData(TopDownScanHeader);
            if (scans != null)
            {
                var filterThreshold = FdrThreshold;
                var filterValues = tsvFile.GetData(FdrHeader);
                if (filterValues == null)
                {
                    filterValues = tsvFile.GetData(EvalueHeader);
                    filterThreshold = EValueThreshold;
                }

                var aset = new AminoAcidSet();

                for (int i = 0; i < peptides.Count; i++)
                {
                    if (peptides[i].Contains('[') || peptides[i].Contains('U') || Convert.ToDouble(filterValues[i]) > filterThreshold) continue;
                    var spectrum = lcms.GetSpectrum(Convert.ToInt32(scans[i]));
                    var spec = spectrum as ProductSpectrum;
                    if (spec == null || spec.ActivationMethod != act) continue;
                    int precursorCharge = Convert.ToInt32(precursorCharges[i]);
                    peptides[i] = Trim(peptides[i]);
                    if (useDecoy)
                    {
                        var shuffled = SimpleStringProcessing.Shuffle(peptides[i]);
                        peptides[i] = SimpleStringProcessing.Mutate(shuffled, NumMutations);
                    }
                    Matches.Add(new SpectrumMatch(peptides[i], spectrum, precursorCharge, new Sequence(peptides[i], aset)));
                }
            }
            else
            {
                scans = tsvFile.GetData(BottomUpScanHeader);
                if (scans == null)
                    throw new FormatException();

                var pepQValues = tsvFile.GetData(PepQValueHeader);
                var formulas = tsvFile.GetData(FormulaHeader);

                for (int i = 0; i < peptides.Count; i++)
                {
                    if (Convert.ToDouble(pepQValues[i]) > PepQValueThreshold) continue;
                    var spectrum = lcms.GetSpectrum(Convert.ToInt32(scans[i]));
                    var spec = spectrum as ProductSpectrum;
                    if (spec == null || spec.ActivationMethod != act) continue;
                    int precursorCharge = Convert.ToInt32(precursorCharges[i]);
                    if (useDecoy)
                    {
                        var shuffled = SimpleStringProcessing.Shuffle(peptides[i]);
                        peptides[i] = SimpleStringProcessing.Mutate(shuffled, NumMutations);
                    }
                    Matches.Add(formulas[i] != null
                        ? new SpectrumMatch(peptides[i], spectrum, precursorCharge, formulas[i])
                        : new SpectrumMatch(peptides[i], spectrum, precursorCharge));
                }
            }
        }

        public List<SpectrumMatch> GetCharge(int charge)
        {
            return (from i in Matches
                    where i.PrecursorCharge == charge
                    select i).ToList();
        }
    }
}
