using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Policy;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.MassFeature;

namespace InformedProteomics.Backend.SequenceTag
{
    public class SequenceTagFinder : SequenceTagGraph<SequenceTagGraphEdge>
    {
        public SequenceTagFinder(ProductSpectrum spec, Tolerance tolerance, int minTagLength = 5, int maxTagLength = 8, AminoAcid[] aminoAcidsArray = null)
            : base(maxTagLength)
        {
            var baseIonTypes = spec.ActivationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;
            var ionTypeFactory = new IonTypeFactory(baseIonTypes, new List<NeutralLoss> { NeutralLoss.NoLoss }, MaxCharge);
            _ionTypes = ionTypeFactory.GetAllKnownIonTypes().ToArray();
            _aminoAcidsArray = aminoAcidsArray ?? AminoAcid.StandardAminoAcidArr;
            _tolerance = tolerance;

            if (_aminoAcidsArray.Length - 1 > Byte.MaxValue) 
                throw new Exception("Too many amino acid types");

            _maxAminoAcidMass = 0d;
            _minAminoAcidMass = 10E4;
            foreach (var aa in _aminoAcidsArray)
            {
                if (aa.Composition.Mass > _maxAminoAcidMass) 
                    _maxAminoAcidMass = aa.Composition.Mass;
                
                if (aa.Composition.Mass < _minAminoAcidMass) 
                    _minAminoAcidMass = aa.Composition.Mass;

            }
            _minTagLength = minTagLength;

            _spectrum = spec;
            _deconvolutedPeaks = Deconvoluter.GetDeconvolutedPeaks(_spectrum, MinCharge, MaxCharge, IsotopeOffsetTolerance, FilteringWindowSize, _tolerance, 0.7);

            SetNodeCount(_deconvolutedPeaks.Count);
            CollectSequenceTagGraphEdges();

            //_candidateSet = new Dictionary<string, SequenceTag>();
            _seqTagSet = new HashSet<SequenceTagString>();
            NumberOfProcessedPaths = 0;
            MaxNumberOfProcessedPaths = 1024;
        }

        public IList<SequenceTagString> GetAllSequenceTagString()
        {
            /*
            var ret = new List<SequenceTagString>();
            foreach (var tag in FindSequenceTags())
            {
                var flankingMass = DeconvolutedPeaks[tag[0].Node1].Mass;
                foreach (var tagStr in tag.GetTagStrings())
                {
                    for (var k = 0; k < 2; k++)
                    {
                        if (k == 0)
                            ret.Add(new SequenceTagString(_spectrum.ScanNum, tagStr, true, flankingMass, _spectrum.ActivationMethod));
                        else
                            ret.Add(new SequenceTagString(_spectrum.ScanNum, SequenceTag.Reverse(tagStr), false, flankingMass, _spectrum.ActivationMethod));
                    }
                }
            }
            return ret;*/

            return FindSequenceTags().ToList();
        }

        public IEnumerable<SequenceTagString> FindSequenceTags()
        {
            var componentSet = ConnnectedComponents();
            //Console.WriteLine("\t# of nodes = {0}, # of edges = {1}, # of components = {2}", GetNodeCount(), GetEdgeCount(), componentSet.Count);
            var startNodeSet = new List<int>();
            foreach (var comp in componentSet)
            {
                if (comp.Count < _minTagLength) continue;

                startNodeSet.AddRange(comp);
            }

            foreach (var node in startNodeSet)
            {
                NumberOfProcessedPaths = 0;
                NumberOfAddedPaths = 0;
                StopFindPath = false;

                var edges = OutEdges(node);
                if (edges.Count < 1) continue;
                FindPaths(node);

                //if (_candidateSet.Count > 32767) break;

                if (_seqTagSet.Count > 32767) break;

            }

            //var candidateList = _candidateSet.Values.ToList();
            //candidateList.Sort();
            //return candidateList;

            return _seqTagSet;
        }

        public List<DeconvolutedPeak> DeconvolutedPeaks { get { return _deconvolutedPeaks; } }

        private readonly List<DeconvolutedPeak> _deconvolutedPeaks;
        private readonly double _maxAminoAcidMass;
        private readonly double _minAminoAcidMass;

        // constant varialbes
        private const int MaxCharge = 15;
        private const int MinCharge = 1;
        private const double FilteringWindowSize = 1.1;
        private const int IsotopeOffsetTolerance = 2;

        private readonly Tolerance _tolerance;
        private readonly AminoAcid[] _aminoAcidsArray;
        private readonly IonType[] _ionTypes;
        private readonly int _minTagLength;
        private readonly ProductSpectrum _spectrum;
        //private readonly Dictionary<string, SequenceTag> _candidateSet;

        private readonly HashSet<SequenceTagString> _seqTagSet;

        private void CollectSequenceTagGraphEdges()
        {
            //var edgeSet = new NodeSet<SequenceTagGraphEdge>();
            for (var i = 0; i < _deconvolutedPeaks.Count; i++)
            {
                var massTh = _tolerance.GetToleranceAsTh(_deconvolutedPeaks[i].Mass);

                for (var j = i + 1; j < _deconvolutedPeaks.Count; j++)
                {
                    if (_deconvolutedPeaks[i].PeakShare(_deconvolutedPeaks[j])) continue;

                    var massGap = _deconvolutedPeaks[j].Mass - _deconvolutedPeaks[i].Mass;
                    var maxMassGap = massGap + massTh;
                    var minMassGap = massGap - massTh;

                    var peakGap = new SequenceTagGraphEdge(i, j, massGap);
                    if (minMassGap > _maxAminoAcidMass) break;
                    if (maxMassGap < _minAminoAcidMass) continue;

                    foreach (var aa in _aminoAcidsArray)
                    {
                        var massError = Math.Abs(peakGap.Mass - aa.Composition.Mass);
                        if (minMassGap < aa.Composition.Mass && aa.Composition.Mass < maxMassGap) peakGap.AddMatchedAminoAcid(aa, massError);
                    }
                    if (peakGap.AminoAcidList.Count > 0) AddEdge(peakGap);
                    //if (peakGap.AminoAcidList.Count > 0) edgeSet.Add(peakGap);
                }
            }

            //var edgeComparer = new SequenceTagGraphEdgeComparer(_deconvolutedPeaks, _tolerance);
            //var clusteredEdges = edgeSet.ConnnectedComponents(edgeComparer);
        }
       
        protected override bool AlreadyUsedPeak(int node)
        {
            foreach (var e in EdgeList)
            {
                if (_deconvolutedPeaks[e.Node1].PeakShare(_deconvolutedPeaks[node])) return true;
                if (_deconvolutedPeaks[e.Node2].PeakShare(_deconvolutedPeaks[node])) return true;
            }
            return false;
        }

        protected int NumberOfAddedPaths;
        protected int MaxNumberOfProcessedPaths;
        protected int NumberOfProcessedPaths;

        protected override bool ProcessPath(IEnumerable<SequenceTagGraphEdge> edges)
        {
            var edgeList = edges.ToList();
            if (edgeList.Count < _minTagLength)
                return true;

            var tag = new SequenceTag(edgeList, _deconvolutedPeaks, _tolerance);

            //tag.Score = GetRankSumScore(tag);
            //if (tag.Score > 0.05) return true;
            //tag.Score = tag.GetTagStrings(_tolerance).Count;
            NumberOfProcessedPaths++;

            var added = false;
            foreach (var tagStr in GetSequenceTagStrings(tag))
            {
                if (_seqTagSet.Add(tagStr)) added = true; 
            }
            if (added) NumberOfAddedPaths++;

            /*
            if (_candidateSet.ContainsKey(tag.HashString))
            {
                //if (_candidateSet[tag.HashString].Score < tag.Score) _candidateSet[tag.HashString] = tag;
                _candidateSet[tag.HashString].Merge(tag.GetTagStrings());
            }
            else
            {
                _candidateSet.Add(tag.HashString, tag);
                NumberOfAddedPaths++;
            }
            */

            //if (NumberOfProcessedPaths > 100 && NumberOfProcessedPaths*0.3 > NumberOfAddedPaths) StopFindPath = true;
            if (NumberOfProcessedPaths > MaxNumberOfProcessedPaths) 
                StopFindPath = true;

            return true;
        }

        private IEnumerable<SequenceTagString> GetSequenceTagStrings(SequenceTag tag)
        {
            var flankingMass = _deconvolutedPeaks[tag[0].Node1].Mass;
            foreach (var tagStr in tag.GetTagStrings())
            {
                for (var k = 0; k < 2; k++)
                {
                    if (k == 0)
                        yield return new SequenceTagString(_spectrum.ScanNum, tagStr, true, flankingMass, _spectrum.ActivationMethod);
                    else
                        yield return new SequenceTagString(_spectrum.ScanNum, SequenceTag.Reverse(tagStr), false, flankingMass, _spectrum.ActivationMethod);
                }
    
            }
        }

        private static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        static SequenceTagFinder()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
        }

        public List<IdentifiedSequenceTag> ExtractExistingSequneceTags(Sequence sequence, int minTagLength = 5)
        {
            var cleavages = sequence.GetInternalCleavages();

            var bIonIndex = new List<int>(); // list of cleavage indices
            var yIonIndex = new List<int>();

            var bIonTypes = new List<IonType>(); // list of ion types
            var yIonTypes = new List<IonType>();

            var index = 0; // cleavage index 
            foreach (var c in cleavages)
            {
                foreach (var ionType in _ionTypes)
                {
                    Ion ion;
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
        /*
               internal class SequenceTagGraphEdgeComparer : INodeComparer<SequenceTagGraphEdge>
               {
                   internal SequenceTagGraphEdgeComparer(List<DeconvolutedPeak> deconvolutedPeaks, Tolerance tolerance)
                   {
                       _deconvolutedPeaks = deconvolutedPeaks;
                       _tolerance = tolerance;
                   }

                   public bool SameCluster(SequenceTagGraphEdge edge1, SequenceTagGraphEdge edge2)
                   {
                       var massTh = _tolerance.GetToleranceAsTh(_deconvolutedPeaks[edge1.Node1].Mass);
                       if (Math.Abs(edge1.Mass - edge2.Mass) > massTh) return false;
                       var startNodeMassDiff = Math.Abs(_deconvolutedPeaks[edge1.Node1].Mass - _deconvolutedPeaks[edge2.Node1].Mass);
                       if (startNodeMassDiff > massTh && Math.Abs(startNodeMassDiff - 1) > massTh) return false;
                       return true;
                   }

                   private List<DeconvolutedPeak> _deconvolutedPeaks;
                   private Tolerance _tolerance;
               }
               */

    }
}