using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.SequenceTag;
using InformedProteomics.FeatureFinding.MassFeature;
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
            ISequenceTagFinder seqTagFinder,
            LcMsPeakMatrix featureFinder,
            FastaDatabase fastaDb,
            Tolerance tolerance,
            AminoAcidSet aaSet,
            CompositeScorerFactory ms2ScorerFactory = null,
            int minMatchedTagLength = DefaultMinMatchedTagLength,
            double maxSequenceMass = 50000.0,
            int minProductIonCharge = 1,
            int maxProductIonCharge = 20)
        {
            _run = run;
            _featureFinder = featureFinder;

            _searchableDb = new SearchableDatabase(fastaDb);

            _tolerance = tolerance;
            _aaSet = aaSet;
            _minMatchedTagLength = minMatchedTagLength;
            _maxSequenceMass = maxSequenceMass;
            _minProductIonCharge = minProductIonCharge;
            _maxProductIonCharge = maxProductIonCharge;
            MinScan = int.MinValue;
            MaxScan = int.MaxValue;
            _ms2ScorerFactory = ms2ScorerFactory;
            _seqTagFinder = seqTagFinder;
        }

        public int MinScan { get; private set; }
        public int MaxScan { get; private set; }
        private readonly CompositeScorerFactory _ms2ScorerFactory;
        private readonly ISequenceTagFinder _seqTagFinder;

        public long NumTags { get { return _seqTagFinder.NumberOfGeneratedTags();  } }

        public void SetDatabase(FastaDatabase fastaDb)
        {
            _searchableDb = new SearchableDatabase(fastaDb);
        }

        public FastaDatabase FastaDatabase { get { return _searchableDb.FastaDatabase; } }
        public IEnumerable<TagSequenceMatch> RunSearch(int ms2ScanNum)
        {
            var spec = _run.GetSpectrum(ms2ScanNum) as ProductSpectrum;
            var scorer = (_ms2ScorerFactory != null) ? _ms2ScorerFactory.GetMs2Scorer(ms2ScanNum) : new CorrMatchedPeakCounter(spec, _tolerance, _minProductIonCharge, _maxProductIonCharge);
            var tags = _seqTagFinder.GetAllSequenceTagString(ms2ScanNum);

            if (spec == null)
            {
                return Enumerable.Empty<TagSequenceMatch>();
            }
            return GetMatches(tags, spec, scorer);
        }

        public void RunSearch()
        {
            Console.WriteLine("Scan\tSequence\tModifications\tMass\tCharge\tScore\tNTermScore\tCTermScore\tProteinName\tStart\tEnd\tProteinLength");

            foreach(var ms2ScanNum in _run.GetScanNumbers(2))
            {
                if (ms2ScanNum >= MinScan && ms2ScanNum <= MaxScan)
                {
                    foreach (var tagSequenceMatch in RunSearch(ms2ScanNum))
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
                                                               FastaDatabase.GetProteinLength(tagSequenceMatch.ProteinName)
                                                               );
                    }
                }
            }
        }

        private IEnumerable<TagSequenceMatch> GetMatches(IEnumerable<SequenceTag> tags, ProductSpectrum spec, IScorer scorer)
        {
            // Match tags against the database
            var proteinsToTags = GetProteinToMatchedTagsMap(tags, _searchableDb, _aaSet, _tolerance, _tolerance);
            //var tagSequenceMatchList = new List<TagSequenceMatch>();

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
                    //var prevScore = double.NegativeInfinity;
                    //foreach (var match in matches.OrderByDescending(m => m.Score))
                    foreach(var match in matches)
                    {
                        var sequence = proteinSequence.Substring(match.StartIndex, match.EndIndex - match.StartIndex);
                        //re-scoring
                        var sequenceObj = Sequence.CreateSequence(sequence, match.ModificationText, _aaSet);
                        match.Score = sequenceObj.GetInternalCleavages().Sum(c => scorer.GetFragmentScore(c.PrefixComposition, c.SuffixComposition));

                        //var numMatches = matchedTag.Length * 2 + match.NTermScore + match.CTermScore;
                        //var score = match.NTermScore + match.CTermScore;
                        //score += (matchedTag.NumReliableNTermFlankingMasses > 0)
                          //  ? matchedTag.Length*CompositeScorer.ScoreParam.Prefix.ConsecutiveMatch
                            //: matchedTag.Length*CompositeScorer.ScoreParam.Suffix.ConsecutiveMatch;

                        // Poisson p-value score
                        //var n = (match.EndIndex - match.StartIndex - 1)*2;
                        //var lambda = numMatches / n;
                        //var pValue = 1 - Poisson.CDF(lambda, numMatches);
                        //var pScore = (pValue > 0) ? - Math.Log(pValue, 2) : 50.0;
                        //if (numMatches < 5) break;
                        //if (prevScore - numMatches > 2) break;
                        //prevScore = numMatches;

                        var pre = match.StartIndex == 0 ? '-' : proteinSequence[match.StartIndex - 1]; // startIndex is inclusive
                        var post = match.EndIndex >= proteinSequence.Length ? '-' : proteinSequence[match.EndIndex]; // endIndex is Exclusive

                        yield return new TagSequenceMatch(sequence, proteinName, match, pre, post);

                        //tagSequenceMatchList.Add(new TagSequenceMatch(sequence, proteinName, match, pre, post));
                    }
                }
            }
            //return tagSequenceMatchList;
        }

        public class TagSequenceMatch
        {
            public TagSequenceMatch(string sequence, string proteinName, TagMatch tagMatch, char pre, char post)
            {
                Sequence = sequence;
                ProteinName = proteinName;
                TagMatch = tagMatch;
                Pre = pre;
                Post = post;
            }

            public char Pre { get; private set; }
            public char Post { get; private set; }

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
                        if (proteinSequence == null)
                            proteinSequence = proteinName;

                        var matchedTagSet = new MatchedTagSet(proteinSequence, aaSet, tolerance, relaxedTolerance);
                        matchedTagSet.Add(matchedTag);
                        proteinsToTags.Add(proteinName, matchedTagSet);
                    }
                }
            }

            return proteinsToTags;
        }

        private readonly ILcMsRun _run;
        //private readonly ProductScorerBasedOnDeconvolutedSpectra _ms2Scorer;
        //private readonly SequenceTagParser _tagParser;
        private SearchableDatabase _searchableDb;

        private readonly LcMsPeakMatrix _featureFinder;
        //private readonly Ms1FtFilter _ms1FtFilter;

        private readonly Tolerance _tolerance;
        private readonly AminoAcidSet _aaSet;

        private readonly double _maxSequenceMass;
        private readonly int _minProductIonCharge;
        private readonly int _maxProductIonCharge;

        private readonly int _minMatchedTagLength;
    }
}
