using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;
using MathNet.Numerics;
using QuickGraph;
using QuickGraph.Algorithms;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.SequenceTag
{
    using SequenceTagGraph = AdjacencyGraph<int, TaggedEdge<int, double>>;

    public class CandidatePath : IComparable<CandidatePath>
    {
        private static double _eTrMax = -99999;
        private static double _eTrMin = 99999;
        private static double _eThMax = -99999;
        private static double _eThMin = 99999;

        public static int CountCreatedCandidates = 0;
        private static SequenceTagGraph _sequenceGraph;
        private static List<DeconvolutedPeak> _deconvolutedPeaks;
        private static Tolerance _tolerance;
        
        public List<int> Path;
        
        public double EtrScore { get; private set; }
        public double EthScore { get; private set; }
        public short[] IntegerMassArray { get; private set; }

        public string HashString { get; private set; }
        
        public static void Initialize(List<DeconvolutedPeak> deconvolutedPeaks, SequenceTagGraph seqGraph, Tolerance tol)
        {
            _eTrMax = -99999;
            _eTrMin = 99999;
            _eThMax = -99999;
            _eThMin = 99999;
            CountCreatedCandidates = 0;
            _sequenceGraph = seqGraph;
            _deconvolutedPeaks = deconvolutedPeaks;
            _tolerance = tol;

        }
        private CandidatePath(List<int> path, short[] intMassArraym, double etr = 0, double eth = 0)
        {
            Path = path;
            EtrScore = etr;
            EthScore = eth;
            IntegerMassArray = intMassArraym;
            var sb = new StringBuilder(1024);
            sb.Append(IntegerMassArray[0]);
            for (var i = 1; i < IntegerMassArray.Length; i++)
            {
                sb.Append(" ");
                sb.Append(IntegerMassArray[i]);
            }

            HashString = sb.ToString();            
        }
        
        public static CandidatePath Create(List<int> path, double etr, double eth)
        {
            var integerMassArray = new short[path.Count - 1];
            //_sequenceGraph = seqGraph;

            CountCreatedCandidates++;
            if (_eTrMax < etr) _eTrMax = etr;
            if (_eThMax < eth) _eThMax = eth;

            if (_eTrMin > etr) _eTrMin = etr;
            if (_eThMin > eth) _eThMin = eth;


            for (int j = 0; j < path.Count - 1; j++)
            {
                var source = path[j];
                var target = path[j + 1];
                var edges = _sequenceGraph.OutEdges(source).Where(outEdge => outEdge.Target.Equals(target)).ToArray();

                integerMassArray[j] = (byte) Math.Round(edges[0].Tag);
            }
            return new CandidatePath(path, integerMassArray, etr, eth);
        }

        public int CompareTo(CandidatePath other)
        {
            return -TotalScore.CompareTo(other.TotalScore);
        }

        public static double NormalizeScore(double s, double mx, double mn)
        {
            return (s - mn) / (mx - mn);
        }

        public static double RescaleScore(double ns, double mx, double mn)
        {
            return mn + ns*(mx - mn);
        }

        public double TotalScore
        {
            get
            {
                if (Math.Abs(_eTrMax - _eTrMin) < 10E-9) return 0;
                
                double score = 0.5*NormalizeScore(EtrScore, _eTrMax, _eTrMin) +
                               NormalizeScore(EthScore, _eThMax, _eThMin);
                return RescaleScore(score, _eThMax, _eThMin);
            }
        }

        public IEnumerable<string> GetPossibleSequnceTag(AminoAcid[] aminoAcidArray) 
        {
            var listOfAminoAcids = new List<char[]>(IntegerMassArray.Length);

            for (var j = 0; j < Path.Count - 1; j++)
            {
                var source = Path[j];
                var target = Path[j + 1];
                //var edges = _sequenceGraph.OutEdges(source).Where(outEdge => outEdge.Target.Equals(target)).ToArray();

                var massTol = _tolerance.GetToleranceAsTh(_deconvolutedPeaks[source].Mass);
                var massDiff = _deconvolutedPeaks[target].Mass - _deconvolutedPeaks[source].Mass;
                var residues = (from aa in aminoAcidArray where Math.Abs(aa.GetMass() - massDiff) < massTol select aa.Residue).ToList();
                listOfAminoAcids.Add(residues.ToArray());
            }

            var indexArray = new int[Path.Count - 1];

            int totalCombinations = listOfAminoAcids.Aggregate(1, (current, x) => current*x.Length);

            for (int e = 0; e < totalCombinations; e++)
            {
                var sb = new StringBuilder();

                for (int i = 0; i < indexArray.Length; i++)
                    sb.Append(listOfAminoAcids[i][indexArray[i]]);

                //increase indexArray
                for (int i = indexArray.Length - 1; i >= 0; i--)
                {
                    if (indexArray[i] == listOfAminoAcids[i].Length - 1) //reached the last item
                    {
                        indexArray[i] = 0;
                    }
                    else
                    {
                        indexArray[i]++;
                        break;
                    }
                }

                yield return sb.ToString();
            }
        }
    }

    public class CandidatePathManger 
    {
        private int _maxCandidates;
        private Dictionary<string, CandidatePath> _candidates;
        private ScoringFunction _scoringFunction;

        public CandidatePathManger(int maxCand, SequenceTagGraph sequenceGraph, ScoringFunction scoreFunc, Tolerance tol)
        {
            _maxCandidates = maxCand;
            CandidatePath.Initialize(scoreFunc.DeconvolutedPeaks ,sequenceGraph, tol);
            _candidates = new Dictionary<string, CandidatePath>();
            _scoringFunction = scoreFunc;
        }

        public void Add(List<int> path)
        {
            var eTh = _scoringFunction.ComputeScoreByHyperGeometricTest(path);
            var eTr = _scoringFunction.ComputeScoreByRankSumTest(path);

            var newCandidate = CandidatePath.Create(path, eTr, eTh);
            
            if (_candidates.ContainsKey(newCandidate.HashString))
            {
                if (_candidates[newCandidate.HashString].TotalScore < newCandidate.TotalScore)
                {
                    _candidates[newCandidate.HashString] = newCandidate;
                }
            }
            else
            {
                _candidates.Add(newCandidate.HashString, newCandidate);
            }
        }

        public void CutoffCandidates()
        {
            var candidateList = _candidates.Values.OrderByDescending(x => x.TotalScore).ToList();
            
            if (candidateList.Count > _maxCandidates)
                for (var i = _maxCandidates; i < candidateList.Count; i++) _candidates.Remove(candidateList[i].HashString);
        }

        public List<CandidatePath> GetRankedCandidates()
        {
            var candidateList = _candidates.Values.OrderByDescending(x => x.TotalScore).ToList();

            if (candidateList.Count > _maxCandidates)
            {
                for (var i = _maxCandidates; i < candidateList.Count; i++) _candidates.Remove(candidateList[i].HashString);
                candidateList.RemoveRange(_maxCandidates, candidateList.Count - _maxCandidates);
            }

            return candidateList;
        }
      
    }

}
