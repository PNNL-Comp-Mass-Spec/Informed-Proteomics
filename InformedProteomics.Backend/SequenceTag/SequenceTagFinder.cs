using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.SequenceTag
{
    public class SequenceTagFinder : SequenceTagGraph<SequenceTagGraphEdge>
    {
        public SequenceTagFinder(ProductSpectrum spec, Tolerance tolerance, int minTagLength = 5, AminoAcid[] aminoAcidsArray = null) : base()
        {
            var baseIonTypes = spec.ActivationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;
            var ionTypeFactory = new IonTypeFactory(baseIonTypes, new List<NeutralLoss> { NeutralLoss.NoLoss }, MaxCharge);
            _ionTypes = ionTypeFactory.GetAllKnownIonTypes().ToArray();
            _aminoAcidsArray = aminoAcidsArray ?? AminoAcid.StandardAminoAcidArr;
            _tolerance = tolerance;
          
            if (_aminoAcidsArray.Length - 1 > Byte.MaxValue) throw new Exception("Too many amino acid types");

            _maxAminoAcidMass = 0d;
            _minAminoAcidMass = 10E4;
            foreach(var aa in _aminoAcidsArray)
            {
                if (aa.Composition.Mass > _maxAminoAcidMass) _maxAminoAcidMass = aa.Composition.Mass;
                if (aa.Composition.Mass < _minAminoAcidMass) _minAminoAcidMass = aa.Composition.Mass;

            }
            _minTagLength = minTagLength;

            _spectrum = spec;
            _deconvolutedPeaks = Deconvoluter.GetDeconvolutedPeaks(_spectrum, MinCharge, MaxCharge, IsotopeOffsetTolerance, FilteringWindowSize, _tolerance, 0.7);
            
            if (_deconvolutedPeaks.Count > 1000) _deconvolutedPeaks = Deconvoluter.FilterOut(_deconvolutedPeaks, 101, 20);

            SetNodeCount(_deconvolutedPeaks.Count);
            CollectSequenceTagGraphEdges();

            //Console.WriteLine("   {0} - {1}", GetNodeCount(), GetEdgeCount());

            _candidateSet = new Dictionary<string, SequenceTag>();

            _minIndexList = new List<int>();
            _maxIndexList = new List<int>();
            _minMassList = new List<double>();
            _maxMassList = new List<double>();
            _rankingList = new List<int[]>();
            MaxTagLen = 8;
        }
        
        public IEnumerable<SequenceTag> FindSequenceTags()
        {
            foreach(var node in Nodes())
            {
                var edges = OutEdges(node);
                if (edges.Count < 1) continue;
                FindPaths(node);

                if (_candidateSet.Count > 32767) break;
            }
            
            var candidateList = _candidateSet.Values.ToList();
            //candidateList.Sort();
            return candidateList;
        }

        public List<DeconvolutedPeak> DeconvolutedPeaks { get { return _deconvolutedPeaks; } }

        private readonly List<DeconvolutedPeak> _deconvolutedPeaks;
        private readonly double _maxAminoAcidMass;
        private readonly double _minAminoAcidMass;

        // constant varialbes
        private const int MaxCharge = 10;
        private const int MinCharge = 1;
        private const double FilteringWindowSize = 1.1;
        private const int IsotopeOffsetTolerance = 2;

        private readonly Tolerance _tolerance;
        private readonly AminoAcid[] _aminoAcidsArray;
        private readonly IonType[] _ionTypes;
        private readonly int _minTagLength;
        private readonly ProductSpectrum _spectrum;
        private readonly Dictionary<string, SequenceTag> _candidateSet;

        private void CollectSequenceTagGraphEdges()
        {
            const double massTh = 0.01;
            for (var i = 0; i < _deconvolutedPeaks.Count; i++)
            {
                //var massTh = Math.Min(0.05, _tolerance.GetToleranceAsTh(_deconvolutedPeaks[i].Mass));
                for (var j = i + 1; j < _deconvolutedPeaks.Count; j++)
                {
                    var peakGap = new SequenceTagGraphEdge(i, j, _deconvolutedPeaks[j].Mass - _deconvolutedPeaks[i].Mass);
                    if (peakGap.Mass > _maxAminoAcidMass + massTh) break;
                    if (peakGap.Mass < _minAminoAcidMass - massTh) continue;

                    foreach (var aa in _aminoAcidsArray)
                    {
                        var massError = Math.Abs(peakGap.Mass - aa.Composition.Mass);
                        if (massError < massTh) peakGap.AddMatchedAminoAcid(aa, massError);
                    }
                    if (peakGap.AminoAcidList.Count > 0) AddEdge(peakGap);
                }
            }
        }

        protected override void ProcessPath(IEnumerable<SequenceTagGraphEdge> edges)
        {
            var tag = new SequenceTag(edges);
            if (tag.Count < _minTagLength) return;

            tag.Score = GetRankSumScore(tag);
            if (_candidateSet.ContainsKey(tag.HashString))
            {
                if (_candidateSet[tag.HashString].Score < tag.Score) _candidateSet[tag.HashString] = tag;
            }
            else
            {
                _candidateSet.Add(tag.HashString, tag);
            }
        }

        private readonly List<double> _minMassList;
        private readonly List<double> _maxMassList;
        private readonly List<int> _minIndexList;
        private readonly List<int> _maxIndexList;
        private readonly List<int[]> _rankingList;

        private static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        static SequenceTagFinder()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
        }
        
        private double GetRankSumScore(SequenceTag seqTag)
        {
            var newRanking = true;
            var node1 = seqTag.First().Node1;
            var node2 = seqTag.Last().Node2;
            var massBuf = seqTag.TotalMass*0.1;
            var minMass = _deconvolutedPeaks[node1].Mass - massBuf;
            var maxMass = _deconvolutedPeaks[node2].Mass + massBuf;

            var minIndex = 0;
            var maxIndex = 0;
            int[] rankings = null;

            for (var i = 0; i < _minMassList.Count; i++)
            {
                if (Math.Abs(_minMassList[i] - minMass) < massBuf && Math.Abs(maxMass - _maxMassList[i]) < massBuf)
                {
                    newRanking = false;
                    minIndex = _minIndexList[i];
                    maxIndex = _maxIndexList[i];
                    rankings = _rankingList[i];
                    break;
                }
            }

            if (newRanking)
            {
                minIndex =_deconvolutedPeaks.BinarySearch(new DeconvolutedPeak(minMass, 0, 0));
                if (minIndex < 0) minIndex = ~minIndex;

                maxIndex = _deconvolutedPeaks.BinarySearch(new DeconvolutedPeak(maxMass, 0, 0));
                if (maxIndex < 0) maxIndex = ~maxIndex;
            
                var intensities = new List<double>();
                for (var i = minIndex; i < maxIndex; i++)
                {
                    intensities.Add(_deconvolutedPeaks[i].Intensity);
                }

                rankings = ArrayUtil.GetRankings(intensities);
                _minIndexList.Add(minIndex);
                _maxIndexList.Add(maxIndex);
                _minMassList.Add(minMass);
                _maxMassList.Add(maxMass);
                _rankingList.Add(rankings);                
            }

            var rankSum = 0d;
            var nMatchedPeaks = 0;
            foreach (var edge in seqTag)
            {
                if (edge.Node1 == node1 && edge.Node1 >= minIndex && edge.Node1 < maxIndex)
                {
                    nMatchedPeaks++;
                    rankSum += rankings[edge.Node1 - minIndex];
                }

                if (edge.Node2 >= minIndex && edge.Node2 < maxIndex)
                {
                    rankSum += rankings[edge.Node2 - minIndex];
                    nMatchedPeaks++;
                }
            }

            var pvalue = FitScoreCalculator.GetRankSumPvalue(maxIndex - minIndex, nMatchedPeaks, rankSum);
            if (pvalue > 0) return -Math.Log(pvalue, 2);
            return 50;
        }

        private double GetHyperGeometricScore(IEnumerable<SequenceTagGraphEdge> seqTag)
        {
            var sequenceTagGraphEdges = seqTag as IList<SequenceTagGraphEdge> ?? seqTag.ToList();
            var node1 = sequenceTagGraphEdges.First().Node1;
            var node2 = sequenceTagGraphEdges.Last().Node2;

            var massRange = _deconvolutedPeaks[node2].Mass - _deconvolutedPeaks[node1].Mass;
            var nObservedPeaks = node2 - node1 + 1;
            var nPossiblePeaks = (int)Math.Round(massRange / 0.01);

            var k1 = sequenceTagGraphEdges.Count + 1;
            var n1 = k1 * 2;

            var pvalue = FitScoreCalculator.GetHyperGeometricPvalue(nPossiblePeaks, nObservedPeaks, n1, k1);
            if (pvalue > 0) return -Math.Log(pvalue, 2);
            return 50;
        }

        public List<IdentifiedSequenceTag> ExtractExistingSequneceTags(Sequence sequence, int minTagLength = 5)
        {
            var cleavages = sequence.GetInternalCleavages();

            var bIonIndex = new List<int>(); // list of cleavage indices
            var yIonIndex = new List<int>();

            var bIonTypes = new List<IonType>(); // list of ion types
            var yIonTypes = new List<IonType>();

            int index = 0; // cleavage index 
            foreach (var c in cleavages)
            {
                foreach (var ionType in _ionTypes)
                {
                    Ion ion = null;
                    if (ionType.IsPrefixIon)
                        ion = ionType.GetIon(c.PrefixComposition);
                    else
                        ion = ionType.GetIon(c.SuffixComposition);

                    var matchedPeaks = _spectrum.ContainsIon(ion, _tolerance, 0.7);

                    if (matchedPeaks)
                    {
                        //string annotation = String.Format("{0}{1}{2}({3}+)", ionType.BaseIonType.Symbol, index, ionType.NeutralLoss.Name, ionType.Charge);
                        //Console.WriteLine(annotation);
                        if (ionType.IsPrefixIon)
                        {
                            if (!bIonIndex.Contains(index))
                            {
                                bIonIndex.Add(index);
                                bIonTypes.Add(ionType);
                            }
                        }
                        else
                        {
                            if (!yIonIndex.Contains(index))
                            {
                                yIonIndex.Add(index);
                                yIonTypes.Add(ionType);
                            }
                        }
                    }
                }
                index++;
            }

            var bTags = IdentifyTagCleavages(bIonIndex.ToArray(), bIonTypes, sequence);
            var seqTags = bTags.Where(t => t.GetLength() >= minTagLength).ToList();

            var yTags = IdentifyTagCleavages(yIonIndex.ToArray(), yIonTypes, sequence);
            seqTags.AddRange(yTags.Where(t => t.GetLength() >= minTagLength));

            return seqTags;
        }

        private IEnumerable<IdentifiedSequenceTag> IdentifyTagCleavages(int[] ionIndex, List<IonType> ionTypes, Sequence sequence)
        {
            var tags = new List<IdentifiedSequenceTag>();
            var tagBegin = 0;
            var tagEnd = 0;

            var tagBeginIdx = 0;

            for (var i = 0; i < ionIndex.Length; i++)
            {
                if (tagBegin == 0)
                {
                    tagBegin = ionIndex[i];
                    tagEnd = tagBegin;
                    tagBeginIdx = i;
                }
                else
                {
                    if (ionIndex[i] == tagEnd + 1)
                        tagEnd++;
                    else
                    {
                        if (tagEnd > tagBegin)
                            tags.Add(new IdentifiedSequenceTag(sequence, tagBegin, tagEnd, ionTypes.GetRange(tagBeginIdx, tagEnd - tagBegin).ToArray()));

                        tagBegin = ionIndex[i];
                        tagEnd = tagBegin;
                        tagBeginIdx = i;
                    }
                }
            }

            if (tagEnd > tagBegin)
                tags.Add(new IdentifiedSequenceTag(sequence, tagBegin, tagEnd, ionTypes.GetRange(tagBeginIdx, tagEnd - tagBegin).ToArray()));

            return tags;
        }

            
    }
}