using System;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.PostProcessing;
using InformedProteomics.TopDown.Scoring;

namespace InformedProteomics.TopDown.TagBasedSearch
{
    public class FeatureBasedTagSearchEngine
    {
        public const int MaxNumProteinMatchesPerTag = 100;
        public const int MinTagLength = 5;

        public FeatureBasedTagSearchEngine(
            LcMsRun run,
            Ms1FtParser featureParser,
            SequenceTagParser tagParser,
            FastaDatabase fastaDb,
            Tolerance tolerance,
            AminoAcidSet aaSet,
            double maxSequenceMass = 50000.0,
            int minProductIonCharge = 1,
            int maxProductIonCharge = 20)
            : this(
                run,
                featureParser,
                null,
                tagParser,
                fastaDb,
                tolerance,
                aaSet,
                maxSequenceMass,
                minProductIonCharge,
                maxProductIonCharge)
        {
        }

        public FeatureBasedTagSearchEngine(
            LcMsRun run,
            Ms1FtParser featureParser,
            ProductScorerBasedOnDeconvolutedSpectra ms2Scorer,
            SequenceTagParser tagParser,
            FastaDatabase fastaDb,
            Tolerance tolerance,
            AminoAcidSet aaSet,
            double maxSequenceMass = 50000.0,
            int minProductIonCharge = 1,
            int maxProductIonCharge = 20)
        {
            _run = run;
            _ms2Scorer = ms2Scorer;
            _featureParser = featureParser;
            _ms1FtFilter = new Ms1FtFilter(run, tolerance, featureParser.Ms1FtFileName);
            _tagParser = tagParser;
            _fastaDb = fastaDb;
            _searchableDb = new SearchableDatabase(fastaDb);
            _tolerance = tolerance;
            _aaSet = aaSet;
            _maxSequenceMass = maxSequenceMass;
            _minProductIonCharge = minProductIonCharge;
            _maxProductIonCharge = maxProductIonCharge;
        }

        public void RunSearch()
        {
            Console.WriteLine("FeatureId\tSequence\tModifications\tMass\tCharge\tScore\tNTermScore\tCTermScore\tProteinName\tStartIndex\tEndIndex\tProteinLength");

            foreach (var feature in _featureParser.GetAllFeatures())
            {
                RunSearch(feature);
            }
        }

        public void RunSearch(ProMexFeature feature)
        {
            var ms2ScanNums = _ms1FtFilter.GetMatchingMs2ScanNums(feature.MonoMass)
                .Where(scanNum => scanNum > feature.MinScan && scanNum < feature.MaxScan)
                .ToArray();

            foreach (var ms2ScanNum in ms2ScanNums)
            {
                _tagParser.GetSequenceTags(ms2ScanNum);
            }
        }

        private readonly ILcMsRun _run;
        private readonly Ms1FtParser _featureParser;
        private readonly Ms1FtFilter _ms1FtFilter;

        private readonly ProductScorerBasedOnDeconvolutedSpectra _ms2Scorer;
        private readonly SequenceTagParser _tagParser;
        private readonly FastaDatabase _fastaDb;
        private readonly SearchableDatabase _searchableDb;

        private readonly Tolerance _tolerance;
        private readonly AminoAcidSet _aaSet;

        private readonly double _maxSequenceMass;
        private readonly int _minProductIonCharge;
        private readonly int _maxProductIonCharge;
    }
}
