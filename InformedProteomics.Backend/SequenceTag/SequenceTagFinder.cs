using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Biology;

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

            //Console.WriteLine("[{0}] # of mz peaks = {1}, # of ms peaks = {2}", _spectrum.ScanNum, _spectrum.Peaks.Length, _deconvolutedPeaks.Count);
            //if (_deconvolutedPeaks.Count > 1000) _deconvolutedPeaks = Deconvoluter.FilterOut(_deconvolutedPeaks, 101, 20);

            SetNodeCount(_deconvolutedPeaks.Count);
            CollectSequenceTagGraphEdges();

            _candidateSet = new Dictionary<string, SequenceTag>();

            //_minIndexList = new List<int>();
            //_maxIndexList = new List<int>();
            //_minMassList = new List<double>();
            //_maxMassList = new List<double>();
            //_rankingList = new List<int[]>();
            //MaxTagLen = 8;

            NumberOfProcessedPaths = 0;
            MaxNumberOfProcessedPaths = 1024;
        }
        
        public IEnumerable<SequenceTag> FindSequenceTags()
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

                if (_candidateSet.Count > 32767) break;

            }
            
            var candidateList = _candidateSet.Values.ToList();
            candidateList.Sort();
            return candidateList;
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
        private readonly Dictionary<string, SequenceTag> _candidateSet;

        private void CollectSequenceTagGraphEdges()
        {
            for (var i = 0; i < _deconvolutedPeaks.Count; i++)
            {
                var massTh = _tolerance.GetToleranceAsTh(_deconvolutedPeaks[i].Mass);

                for (var j = i + 1; j < _deconvolutedPeaks.Count; j++)
                {
                    if (_deconvolutedPeaks[i].PeakShare(_deconvolutedPeaks[j])) continue;
                    //if (Math.Abs(_deconvolutedPeaks[j].Charge - _deconvolutedPeaks[i].Charge) > 1) continue;
                    
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
                }
            }
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
            if (edgeList.Count < _minTagLength) return true;
            var tag = new SequenceTag(edges, _deconvolutedPeaks, _tolerance);
            
            //tag.Score = GetRankSumScore(tag);
            //if (tag.Score > 0.05) return true;
            //tag.Score = tag.GetTagStrings(_tolerance).Count;
            NumberOfProcessedPaths++;
            
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
            
            //if (NumberOfProcessedPaths > 100 && NumberOfProcessedPaths*0.3 > NumberOfAddedPaths) StopFindPath = true;
            if (NumberOfProcessedPaths > MaxNumberOfProcessedPaths) StopFindPath = true;

            return true;
        }
        
        private static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        static SequenceTagFinder()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
        }
        /*
        private readonly List<double> _minMassList;
        private readonly List<double> _maxMassList;
        private readonly List<int> _minIndexList;
        private readonly List<int> _maxIndexList;
        private readonly List<int[]> _rankingList;
        private double GetRankSumScore(SequenceTag seqTag)
        {
            var newRanking = true;
            var node1 = seqTag.First().Node1;
            var node2 = seqTag.Last().Node2;
            var massBuf = seqTag.TotalMass*0.1;
            //const double massBuf = 100d;

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
                minIndex =_deconvolutedPeaks.BinarySearch(new DeconvolutedPeak(minMass, 0));
                if (minIndex < 0) minIndex = ~minIndex;

                maxIndex = _deconvolutedPeaks.BinarySearch(new DeconvolutedPeak(maxMass, 0));
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

            var pvalue = FitScoreCalculator.GetRankSumPvalue(maxIndex - minIndex + 1, nMatchedPeaks, rankSum);

            return pvalue;
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
        }*/

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