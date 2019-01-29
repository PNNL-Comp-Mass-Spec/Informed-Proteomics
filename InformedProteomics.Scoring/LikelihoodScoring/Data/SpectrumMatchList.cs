using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Scoring.LikelihoodScoring.FileReaders;

namespace InformedProteomics.Scoring.LikelihoodScoring.Data
{
    public class SpectrumMatchList: List<SpectrumMatch>
    {
        public int MaxCharge { get; }
        public bool Decoy { get; }
//        public ActivationMethod Act { get; private set; }
        public DataFileFormat SequenceFormat { get; }

        public SpectrumMatchList(LazyLcMsRun lcms, string tsvFile, DataFileFormat sequenceFormat,
                                 int maxCharge=0, bool useDecoy=false)
        {
            Decoy = useDecoy;
            MaxCharge = maxCharge;
            SequenceFormat = sequenceFormat;
            var fileReader = DataFileReaderFactory.Create(SequenceFormat, tsvFile, useDecoy, lcms);
            var specMatches = fileReader.Read();
            AddRange(from i in specMatches
                     where i.PrecursorCharge <= maxCharge
                     select i);
        }

        public SpectrumMatchList(string fileName, int maxCharge = 0, bool useDecoy = false)
        {
            Decoy = useDecoy;
            MaxCharge = maxCharge;
            SequenceFormat = DataFileFormat.Mgf;
            var fileReader = DataFileReaderFactory.Create(SequenceFormat, fileName, useDecoy);
            var specMatches = fileReader.Read();
            AddRange(from i in specMatches
                     where i.PrecursorCharge <= maxCharge
                     select i);
        }

        public SpectrumMatchList(int maxCharge = 0, bool useDecoy = false)
        {
            Decoy = useDecoy;
            MaxCharge = maxCharge;
            SequenceFormat = DataFileFormat.IcBottomUp;
        }

        public void AddMatch(SpectrumMatch match)
        {
            var newMatch = new SpectrumMatch(match, Decoy);
            Add(newMatch);
        }

        public void AddMatches(IEnumerable<SpectrumMatch> matches)
        {
            foreach (var match in matches)
            {
                AddMatch(match);
            }
        }

        public void FilterSpectra(double windowWidth=100, int retentionCount=6)
        {
            var filteredList = new SpectrumMatchList();
            foreach (var match in this)
            {
                var spectrum = SpectrumFilter.GetFilteredSpectrum(match.Spectrum, windowWidth, retentionCount);
                filteredList.Add(new SpectrumMatch(match.Peptide, SequenceFormat, spectrum, match.ScanNum, match.PrecursorCharge, match.Decoy));
            }
            Clear();
            AddRange(filteredList);
        }

        public SpectrumMatchList GetCharge(int charge)
        {
            var chargeMatchList = new SpectrumMatchList();
            chargeMatchList.AddRange(from i in this
                    where i.PrecursorCharge == charge
                    select i);
            return chargeMatchList;
        }
    }
}
