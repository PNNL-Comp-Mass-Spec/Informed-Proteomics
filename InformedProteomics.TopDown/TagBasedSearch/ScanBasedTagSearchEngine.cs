using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Quantification;
using InformedProteomics.TopDown.PostProcessing;
using InformedProteomics.TopDown.Scoring;

namespace InformedProteomics.TopDown.TagBasedSearch
{
    public class ScanBasedTagSearchEngine
    {
        public const int MaxNumProteinMatchesPerTag = 100;
        public const int DefaultMinMatchedTagLength = 6;

        public ScanBasedTagSearchEngine(
            LcMsRun run, 
            SequenceTagParser tagParser, 
            FastaDatabase fastaDb, 
            Tolerance tolerance, 
            AminoAcidSet aaSet, 
            int minMatchedTagLength = DefaultMinMatchedTagLength,
            double maxSequenceMass = 50000.0, 
            int minProductIonCharge = 1, 
            int maxProductIonCharge = 20): this(
            run,
            null,
            tagParser,
            fastaDb,
            tolerance,
            aaSet,
            minMatchedTagLength,
            maxSequenceMass,
            minProductIonCharge,
            maxProductIonCharge)
        {
        }

        public ScanBasedTagSearchEngine(
            LcMsRun run,
            ProductScorerBasedOnDeconvolutedSpectra ms2Scorer,
            SequenceTagParser tagParser,
            FastaDatabase fastaDb,
            Tolerance tolerance,
            AminoAcidSet aaSet,
            int minMatchedTagLength = DefaultMinMatchedTagLength,
            double maxSequenceMass = 50000.0,
            int minProductIonCharge = 1,
            int maxProductIonCharge = 20)
        {
            _run = run;
            _ms2Scorer = ms2Scorer;
            _featureFinder = new TargetMs1FeatureMatrix(run);
            _tagParser = tagParser;
            _fastaDb = fastaDb;
            _searchableDb = new SearchableDatabase(fastaDb);
            _tolerance = tolerance;
            _aaSet = aaSet;
            _minMatchedTagLength = minMatchedTagLength;
            _maxSequenceMass = maxSequenceMass;
            _minProductIonCharge = minProductIonCharge;
            _maxProductIonCharge = maxProductIonCharge;
            MinScan = int.MinValue;
            MaxScan = int.MaxValue;
        }


        public int MinScan { get; set; }
        public int MaxScan { get; set; }

        public void RunSearch()
        {
            Console.WriteLine("Scan\tSequence\tModifications\tMass\tCharge\tScore\tNTermScore\tCTermScore\tProteinName\tStartIndex\tEndIndex\tProteinLength");
            foreach (var ms2ScanNum in _tagParser.GetScanNums())
            {
                if(ms2ScanNum >= MinScan && ms2ScanNum <= MaxScan) RunSearch(ms2ScanNum);
            }
        }

        public void RunSearch(int ms2ScanNum)
        {
            var spec = _run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            if (spec == null) return;

            var scorer = _ms2Scorer != null
                ? _ms2Scorer.GetMs2Scorer(ms2ScanNum)
                : new CorrMatchedPeakCounter(spec, _tolerance, _minProductIonCharge, _maxProductIonCharge);

            var tags = _tagParser.GetSequenceTags(ms2ScanNum);

            foreach (var tagSequenceMatch in GetMatches(tags, spec, scorer))
            {
                var tagMatch = tagSequenceMatch.TagMatch;
                Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}",
                    ms2ScanNum,
                    tagSequenceMatch.Sequence,
                    tagMatch.Modifications,
                    tagMatch.Mass,
                    tagMatch.Charge,
                    tagMatch.Score,
                    tagMatch.NTermScore,
                    tagMatch.CTermScore,
                    tagSequenceMatch.ProteinName,
                    tagMatch.StartIndex,
                    tagMatch.EndIndex,
                    _fastaDb.GetProteinLength(tagSequenceMatch.ProteinName)
                    );
            }
        }

        public IList<TagSequenceMatch> GetMatches(IEnumerable<SequenceTag> tags, ProductSpectrum spec, IScorer scorer)
        {
            // Match tags against the database
            var proteinsToTags = GetProteinToMatchedTagsMap(tags, _searchableDb, _aaSet, _tolerance, _tolerance);

            var tagSequenceMatchList = new List<TagSequenceMatch>();

            // Extend matches
            foreach (var entry in proteinsToTags)
            {
                var proteinName = entry.Key;
                var matchedTagSet = entry.Value;
                var proteinSequence = matchedTagSet.Sequence;

                var tagFinder = new TagMatchFinder(spec, scorer, _featureFinder, proteinSequence, _tolerance, _aaSet, _maxSequenceMass);

                foreach (var matchedTag in matchedTagSet.Tags)
                {
                    if (matchedTag.Length < _minMatchedTagLength) continue;
                    if (matchedTag.NTermFlankingMass == null && matchedTag.CTermFlankingMass == null) continue;

                    var matches = tagFinder.FindMatches(matchedTag).ToArray();
                    var prevScore = double.NegativeInfinity;
                    foreach (var match in matches.OrderByDescending(m => m.Score))
                    {
                        var sequence = proteinSequence.Substring(match.StartIndex, match.EndIndex - match.StartIndex);
                        var numMatches = matchedTag.Length * 2 + match.NTermScore + match.CTermScore;

                        // Poisson p-value score
                        //var n = (match.EndIndex - match.StartIndex - 1)*2;
                        //var lambda = numMatches / n;
                        //var pValue = 1 - Poisson.CDF(lambda, numMatches);
                        //var pScore = (pValue > 0) ? -Math.Log(pValue, 2) : 50.0;

                        if (numMatches < 10) break;
                        if (prevScore - numMatches > 2) break;
                        prevScore = numMatches;

                        tagSequenceMatchList.Add(new TagSequenceMatch(sequence, proteinName, match));
                    }
                }
            }
            return tagSequenceMatchList;
        }

        public class TagSequenceMatch
        {
            public TagSequenceMatch(string sequence, string proteinName, TagMatch tagMatch)
            {
                Sequence = sequence;
                ProteinName = proteinName;
                TagMatch = tagMatch;
            }

            public string Sequence { get; private set; }
            public string ProteinName { get; private set; }
            public TagMatch TagMatch { get; private set; }
        }

        public static Dictionary<string, MatchedTagSet> GetProteinToMatchedTagsMap(
            IEnumerable<SequenceTag> tags, 
            SearchableDatabase searchableDb, 
            AminoAcidSet aaSet, 
            Tolerance tolerance,
            Tolerance relaxedTolerance)
        {
            var fastaDb = searchableDb.FastaDatabase;
            var proteinsToTags = new Dictionary<string, MatchedTagSet>();
            foreach (var tag in tags)
            {
                var matchedIndices = searchableDb.FindAllMatchedSequenceIndices(tag.Sequence).ToArray();
                if (matchedIndices.Length > MaxNumProteinMatchesPerTag) continue;
                foreach (var index in matchedIndices)
                {
                    var proteinName = fastaDb.GetProteinName(index);
                    var startIndex = fastaDb.GetZeroBasedPositionInProtein(index);
                    var mass = aaSet.GetComposition(tag.Sequence).Mass;
                    var matchedTag = new MatchedTag(tag, startIndex) { Mass = mass };
                    MatchedTagSet existingMatchedTagSet;
                    if (proteinsToTags.TryGetValue(proteinName, out existingMatchedTagSet))
                    {
                        existingMatchedTagSet.Add(matchedTag);
                    }
                    else
                    {
                        var proteinSequence = fastaDb.GetProteinSequence(proteinName);
                        var matchedTagSet = new MatchedTagSet(proteinSequence, aaSet, tolerance, relaxedTolerance);
                        matchedTagSet.Add(matchedTag);
                        proteinsToTags.Add(proteinName, matchedTagSet);
                    }
                }
            }

            return proteinsToTags;
        }

        private readonly ILcMsRun _run;
        private readonly ProductScorerBasedOnDeconvolutedSpectra _ms2Scorer;
        private readonly SequenceTagParser _tagParser;
        private readonly FastaDatabase _fastaDb;
        private readonly SearchableDatabase _searchableDb;

        private readonly TargetMs1FeatureMatrix _featureFinder; 

        private readonly Tolerance _tolerance;
        private readonly AminoAcidSet _aaSet;

        private readonly double _maxSequenceMass;
        private readonly int _minProductIonCharge;
        private readonly int _maxProductIonCharge;

        private readonly int _minMatchedTagLength;
    }
}
