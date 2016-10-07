using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Policy;
using System.Text;
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
            _deconvolutedPeaks = Deconvoluter.GetDeconvolutedPeaks(_spectrum.Peaks, MinCharge, MaxCharge, IsotopeOffsetTolerance, 1.1, _tolerance, 0.7);

            SetNodeCount(_deconvolutedPeaks.Count);
            CollectSequenceTagGraphEdges();

            _seqTagSet = new HashSet<SequenceTag>();
            NumberOfProcessedPaths = 0;
            MaxNumberOfProcessedPaths = 1024;
        }

        public IList<SequenceTag> GetAllSequenceTagString()
        {
            return FindSequenceTags().ToList();
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

                //if (_candidateSet.Count > 32767) break;
                if (_seqTagSet.Count > 32767) break;
            }

            //var candidateList = _candidateSet.Values.ToList();
            //candidateList.Sort();
            //return candidateList;

            return _seqTagSet;
        }

        //public List<DeconvolutedPeak> DeconvolutedPeaks { get { return _deconvolutedPeaks; } }

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
        private readonly HashSet<SequenceTag> _seqTagSet;

        private static readonly BaseIonType[] BaseIonTypesCID, BaseIonTypesETD;
        static SequenceTagFinder()
        {
            BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
        }

        private void CollectSequenceTagGraphEdges()
        {
            for (var i = 0; i < _deconvolutedPeaks.Count; i++)
            {
                var massTh = _tolerance.GetToleranceAsTh(_deconvolutedPeaks[i].Mass);

                for (var j = i + 1; j < _deconvolutedPeaks.Count; j++)
                {
                    //if (_deconvolutedPeaks[i].PeakShare(_deconvolutedPeaks[j])) continue;

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
       /*
        protected override bool AlreadyUsedPeak(int node)
        {
            foreach (var e in EdgeList)
            {
                if (_deconvolutedPeaks[e.Node1].PeakShare(_deconvolutedPeaks[node])) return true;
                if (_deconvolutedPeaks[e.Node2].PeakShare(_deconvolutedPeaks[node])) return true;
            }
            return false;
        }*/

        protected int NumberOfAddedPaths;
        protected int MaxNumberOfProcessedPaths;
        protected int NumberOfProcessedPaths;

        protected override bool ProcessPath(IEnumerable<SequenceTagGraphEdge> edges)
        {
            var edgeList = edges.ToList();
            if (edgeList.Count < _minTagLength)
                return true;

            NumberOfProcessedPaths++;

            var added = false;
            /*var tag = new SequenceTag(edgeList, _deconvolutedPeaks, _tolerance);
            foreach (var tagStr in GetSequenceTagStrings(tag))
            {
                if (_seqTagSet.Add(tagStr)) added = true;
            }*/
            foreach (var seqTag in EnumerateAllSequenceTags(edgeList))
            {
                if (_seqTagSet.Add(seqTag)) added = true;
            }

            if (added) NumberOfAddedPaths++;

            if (NumberOfProcessedPaths > MaxNumberOfProcessedPaths)
                StopFindPath = true;

            return true;
        }
        /*
        private IEnumerable<SequenceTag> GetSequenceTagStrings(SequenceTag tag)
        {
            var flankingMass = _deconvolutedPeaks[tag[0].Node1].Mass;
            foreach (var tagStr in tag.GetTagStrings())
            {
                for (var k = 0; k < 2; k++)
                {
                    if (k == 0)
                        yield return new SequenceTag(_spectrum.ScanNum, tagStr, true, flankingMass, _spectrum.ActivationMethod);
                    else
                        yield return new SequenceTag(_spectrum.ScanNum, SequenceTag.Reverse(tagStr), false, flankingMass, _spectrum.ActivationMethod);
                }
            }
        }*/

        private IEnumerable<SequenceTag> EnumerateAllSequenceTags(List<SequenceTagGraphEdge> edges)
        {
            var tagLength = edges.Count;
            var listOfAminoAcids = new List<List<AminoAcid>>(tagLength);

            for (var j = 0; j < tagLength; j++)
            {
                listOfAminoAcids.Add(edges[j].AminoAcidList);
            }

            var indexArray = new int[tagLength];
            var totalCombinations = listOfAminoAcids.Aggregate(1, (current, x) => current * x.Count);
            //TagStrings = new HashSet<string>();

            var massTh = _tolerance.GetToleranceAsTh(_deconvolutedPeaks[edges[0].Node1].Mass);
            var flankingMass = _deconvolutedPeaks[edges[0].Node1].Mass;

            for (var e = 0; e < totalCombinations; e++)
            {
                var sb = new StringBuilder();
                var mass = 0d;
                for (var i = 0; i < indexArray.Length; i++)
                {
                    sb.Append(listOfAminoAcids[i][indexArray[i]].Residue);
                    mass += listOfAminoAcids[i][indexArray[i]].Mass;
                }

                var massGap = _deconvolutedPeaks[edges[indexArray.Length - 1].Node2].Mass - _deconvolutedPeaks[edges[0].Node1].Mass;
                var massError = Math.Abs(massGap - mass);

                if (massError < massTh)
                {
                    //TagStrings.Add(sb.ToString());
                    var tagStr = sb.ToString();
                    yield return new SequenceTag(_spectrum.ScanNum, tagStr, true, flankingMass, _spectrum.ActivationMethod);
                    yield return new SequenceTag(_spectrum.ScanNum, Reverse(tagStr), false, flankingMass, _spectrum.ActivationMethod);
                }

                //increase indexArray
                for (var i = indexArray.Length - 1; i >= 0; i--)
                {
                    if (indexArray[i] == listOfAminoAcids[i].Count - 1) //reached the last item
                    {
                        indexArray[i] = 0;
                    }
                    else
                    {
                        indexArray[i]++;
                        break;
                    }
                }
            }
            //return TagStrings;
        }

        public static string Reverse(string s)
        {
            char[] charArray = s.ToCharArray();
            Array.Reverse(charArray);
            return new string(charArray);
        }

        /*
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
        */
    }
}