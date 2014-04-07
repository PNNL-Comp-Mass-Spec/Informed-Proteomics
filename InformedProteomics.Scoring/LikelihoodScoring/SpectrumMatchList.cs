using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
    public class SpectrumMatchList: List<SpectrumMatch>
    {
        private const string QValueHeader = "Qvalue";
        private const string EvalueHeader = "E-value";
        private const string TopDownPeptideHeader = "Annotation";
        private const double QValueThreshold = 0.01;
        private const double EValueThreshold = 0.1;

        private const string BottomUpPeptideHeader = "Peptide";
        private const string PepQValueHeader = "PepQValue";
        private const string FormulaHeader = "Formula";
        private const double PepQValueThreshold = 0.01;

        private const string ScanHeader = "ScanNum";
        private const string PrecursorChargeHeader = "Charge";
        private const int NumMutations = 3;

        public bool Decoy { get; private set; }
        public ActivationMethod Act { get; private set; }

        public SpectrumMatchList(LcMsRun lcms, TsvFileParser tsvFile, ActivationMethod act, bool useDecoy=false)
        {
            Decoy = useDecoy;
            Act = act;
            AddMatchesFromFile(lcms, tsvFile);
        }

        public SpectrumMatchList(ActivationMethod act, bool useDecoy=false)
        {
            Decoy = useDecoy;
            Act = act;
        }

        public void AddMatchesFromFile(LcMsRun lcms, TsvFileParser tsvFile)
        {
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
                    var spectrum = lcms.GetSpectrum(scanNum);
                    var spec = spectrum as ProductSpectrum;
                    if (spec == null || spec.ActivationMethod != Act) continue;
                    int precursorCharge = Convert.ToInt32(precursorCharges[i]);
                    if (Decoy)
                    {
                        var shuffled = SimpleStringProcessing.Shuffle(peptides[i]);
                        peptides[i] = SimpleStringProcessing.Mutate(shuffled, NumMutations);
                    }
                    Add(new SpectrumMatch(peptides[i], spectrum, scanNum, precursorCharge, new Sequence(peptides[i], aset)));
                }
            }
            else
            {
                peptides = tsvFile.GetData(BottomUpPeptideHeader);
                if (scans == null)
                    throw new FormatException();

                var pepQValues = tsvFile.GetData(PepQValueHeader);
                var formulas = tsvFile.GetData(FormulaHeader);

                var peptideSet = new HashSet<string>();

                for (int i = 0; i < peptides.Count; i++)
                {
                    if (Convert.ToDouble(pepQValues[i]) > PepQValueThreshold || peptideSet.Contains(peptides[i])) continue;
                    peptideSet.Add(peptides[i]);
                    var scanNum = Convert.ToInt32(scans[i]);
                    var spectrum = lcms.GetSpectrum(scanNum);
                    var spec = spectrum as ProductSpectrum;
                    if (spec == null || spec.ActivationMethod != Act) continue;
                    int precursorCharge = Convert.ToInt32(precursorCharges[i]);
                    if (Decoy)
                    {
                        var shuffled = SimpleStringProcessing.Shuffle(peptides[i]);
                        peptides[i] = SimpleStringProcessing.Mutate(shuffled, NumMutations);
                    }
                    Add((formulas != null && formulas[i] != null)
                        ? new SpectrumMatch(peptides[i], spectrum, scanNum, precursorCharge, formulas[i])
                        : new SpectrumMatch(peptides[i], spectrum, scanNum, precursorCharge));
                }
            }
        }

        public SpectrumMatch GetScan(int scanNum)
        {
            return (from i in this
                    where i.ScanNum == scanNum
                    select i).FirstOrDefault();
        }

        public List<SpectrumMatch> GetCharge(int charge)
        {
            return (from i in this
                    where i.PrecursorCharge == charge
                    select i).ToList();
        }
    }
}
