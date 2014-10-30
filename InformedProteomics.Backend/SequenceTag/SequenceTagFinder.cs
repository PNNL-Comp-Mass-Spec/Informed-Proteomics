using System;
using System.Collections.Generic;
using System.Drawing.Text;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Utils;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Biology;

using QuickGraph;
using QuickGraph.Algorithms;
using QuickGraph.Algorithms.ConnectedComponents;

namespace InformedProteomics.Backend.SequenceTag
{
    using SequenceTagGraphEdge  = TaggedEdge<int, double>;
    using SequenceTagGraph      = AdjacencyGraph<int, TaggedEdge<int, double>>; //This is actually the "Spectrum Graph"
    
    public class SequenceTagFinder
    {
        public List<DeconvolutedPeak> DeconvolutedPeaks { get; private set; }
        public int[] PeakRankings { get; private set; }
        
        public double[] MassArray { get; private set; }
        public double MaxAminoAcidMass { get; private set; }
        public double MinAminoAcidMass { get; private set; }

        public static int MaxCharge = 10;
        public static int MinCharge = 1;
        
        private const double FilteringWindowSize = 1.1;
        private const int IsotopeOffsetTolerance = 2;
        private const double ProductTolerancePpm = 10;

        private Tolerance _massTolerance = new Tolerance(5);

        public Spectrum SpectrumData { get; private set; }
        public SequenceTagGraph SequenceGraph { get; private set; }
        
        public AminoAcid[] AminoAcidsArray { get; private set; }
        public IonType[] IonTypes { get; private set; }
        
        private const int MaxCheckedCandidates = 1000000;
        private CandidatePathManger _candidatePathMgr;
        
        public SequenceTagFinder(Spectrum spectrum, AminoAcid[] aminoAcidsArray = null, IonType[] ionTypes = null)
        {
            AminoAcidsArray = aminoAcidsArray ?? AminoAcid.StandardAminoAcidArr;

            if (ionTypes == null)
            {
                var baseIonTypes = new List<BaseIonType> { BaseIonType.B, BaseIonType.Y };
                var ionTypeFactory = new IonTypeFactory(baseIonTypes, new List<NeutralLoss> { NeutralLoss.NoLoss }, MaxCharge);
                IonTypes = ionTypeFactory.GetAllKnownIonTypes().ToArray();
            }
            else
            {
                IonTypes = ionTypes;
            }
            
            SpectrumData = spectrum;
            MassArray = new double[AminoAcidsArray.Length];

            for (int i = 0; i < AminoAcidsArray.Length; i++)
            {
                MassArray[i] = AminoAcidsArray[i].Composition.Mass;
            }
            MaxAminoAcidMass = MassArray.Max() + 1;
            MinAminoAcidMass = MassArray.Min() - 1;
        }

        private const double RescalingConstantHighPrecision = Constants.RescalingConstantHighPrecision;
        public static int GetBinNumber(double mass)
        {
            return (int)Math.Round(mass * RescalingConstantHighPrecision);
        }


        public List<string> GetSequenceTags(int minimumTagLength = 5, int maximumTags = 100)
        {
            PreProcessingSpectrum();
            ConstructTagGraphs();

            var resultTags = new List<string>();
            var allPaths = new List<List<int>>();
            var tmp = new List<int>();
            var vertices = SequenceGraph.Vertices.ToArray();
            var _edges = SequenceGraph.Edges.ToArray();

            var candidates = GenerateAllPaths(maximumTags, minimumTagLength);

            int nSeqTags = 0;

            foreach (var candidate in candidates)
            {
                foreach (var s  in candidate.GetPossibleSequnceTag(AminoAcidsArray))
                {
                    resultTags.Add(s);
                    resultTags.Add(ReverseString(s));
                    nSeqTags += 2;
                }

                if (nSeqTags > maximumTags) break;
            }

            return resultTags;
        }

        private static string ReverseString(string s)
        {
            char[] charArray = s.ToCharArray();
            Array.Reverse(charArray);
            return new string(charArray);
        }

        public void PreProcessingSpectrum()
        {
            var deconvolutedPeaks = Deconvoluter.GetDeconvolutedPeaks(SpectrumData, MinCharge, MaxCharge, IsotopeOffsetTolerance, FilteringWindowSize, new Tolerance(ProductTolerancePpm), 0.7);
            
            var binHash = new HashSet<int>();
            DeconvolutedPeaks = new List<DeconvolutedPeak>();

            foreach (var deconvolutedPeak in deconvolutedPeaks)
            {
                var mass = deconvolutedPeak.Mass;
                var binNum = GetBinNumber(mass);
                if (!binHash.Add(binNum)) continue;

                DeconvolutedPeaks.Add(deconvolutedPeak);
            }

            // compute peak rankings for scoring
            var sorted = DeconvolutedPeaks.Select((x, i) => new KeyValuePair<DeconvolutedPeak, int>(x, i)).OrderByDescending(x => x.Key.Intensity).ToList();
            PeakRankings = new int[sorted.Count];
            for(var ranking = 1; ranking < sorted.Count; ranking++)
            {
                PeakRankings[sorted[ranking - 1].Value] = ranking;
            }
        }

        private IEnumerable<SequenceTagGraphEdge> CollectSequenceTagGraphEdges()
        {
            var edges = new List<SequenceTagGraphEdge>();

            for (var i = 0; i < DeconvolutedPeaks.Count; i++)
            {
                for (var j = i + 1; j < DeconvolutedPeaks.Count; j++)
                {
                    double massDiff = DeconvolutedPeaks[j].Mass - DeconvolutedPeaks[i].Mass;

                    if (massDiff > MaxAminoAcidMass) break;
                    if (massDiff < MinAminoAcidMass) continue;

                    if (MassArray.Any(t => Math.Abs(massDiff - t) < _massTolerance.GetToleranceAsTh(DeconvolutedPeaks[i].Mass)))
                    {
                        edges.Add(new SequenceTagGraphEdge(i, j, massDiff));
                    }
                }
            }
            return edges;
        }

        public void ConstructTagGraphs()
        {
            var edges = CollectSequenceTagGraphEdges();
            SequenceGraph = new SequenceTagGraph();

            foreach (var e in edges)
                SequenceGraph.AddVerticesAndEdge(e);
        }

        public IEnumerable<CandidatePath> GenerateAllPaths(int maxPaths = 10000, int minTagLength = 5)
        {
            _candidatePathMgr = new CandidatePathManger(maxPaths, SequenceGraph, new ScoringFunction(DeconvolutedPeaks), _massTolerance);
            var allPaths = new List<List<int>>();

            //foreach (var v in SequenceGraph.Vertices)
            foreach(var v in SequenceGraph.Roots())
            {
                if (SequenceGraph.IsOutEdgesEmpty(v)) continue;
                var tmp = new List<int>();

                FindAllPathsAt(v, ref allPaths, tmp, minTagLength);
            }
            FlushCandidates( ref allPaths );

            return _candidatePathMgr.GetRankedCandidates();
        }
        
        
        private void FlushCandidates(ref List<List<int>> allPaths)
        {
            foreach (var path in allPaths)
            {
                _candidatePathMgr.Add(path);
            }
            
            _candidatePathMgr.CutoffCandidates();
            allPaths.Clear();

            //if (CandidatePath.CountCreatedCandidates % 10000 == 0)
                //Console.WriteLine("Total Checked Candidates = {0}", CandidatePath.CountCreatedCandidates);
        }
        
        private void FindAllPathsAt(int v, ref List<List<int>> allPaths, List<int> tmp, int minTagLength = 5)
        {
            tmp.Add(v);
            
            if (SequenceGraph.IsOutEdgesEmpty(v))
            {
                if (tmp.Count >= minTagLength + 1)
                {
                    allPaths.Add(tmp);
                    if (allPaths.Count > 100000)
                    {
                        FlushCandidates(ref allPaths);
                        if (CandidatePath.CountCreatedCandidates > MaxCheckedCandidates) return;
                    }
                }

                return;
            }

            foreach (var e in SequenceGraph.OutEdges(v))
            {
                var tmp2 = new List<int>(tmp);
                int v2 = e.GetOtherVertex(v);
                FindAllPathsAt(v2, ref allPaths, tmp2, minTagLength);

                if (CandidatePath.CountCreatedCandidates > MaxCheckedCandidates) return;
            }
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
                foreach (var ionType in IonTypes)
                {
                    Ion ion = null;
                    if (ionType.IsPrefixIon)
                        ion = ionType.GetIon(c.PrefixComposition);
                    else
                        ion = ionType.GetIon(c.SuffixComposition);

                    var matchedPeaks = SpectrumData.ContainsIon(ion, new Tolerance(ProductTolerancePpm), 0.7);

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


        private static IEnumerable<IdentifiedSequenceTag> IdentifyTagCleavages(int[] ionIndex, List<IonType> ionTypes, Sequence sequence)
        {
            var tags = new List<IdentifiedSequenceTag>();
            int tagBegin = 0;
            int tagEnd = 0;

            int tagBeginIdx = 0;

            for (int i = 0; i < ionIndex.Length; i++)
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