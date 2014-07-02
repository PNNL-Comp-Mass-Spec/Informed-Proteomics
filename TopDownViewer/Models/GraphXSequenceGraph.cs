using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Database;

namespace InformedProteomics.TopDownViewer.Models
{
    public class GraphXSequenceGraph: SequenceGraph
    {
        public DataGraph DataGraph { get; private set; }
        public static GraphXSequenceGraph Create(AminoAcidSet aaSet, string annotation)
        {
            const char delimiter = (char)FastaDatabase.Delimiter;
            if (annotation == null || !Regex.IsMatch(annotation, @"^[A-Z" + delimiter + @"]\.[A-Z]+\.[A-Z" + delimiter + @"]$")) return null;

            var nTerm = annotation[0] == FastaDatabase.Delimiter
                                  ? AminoAcid.ProteinNTerm
                                  : AminoAcid.PeptideNTerm;
            var cTerm = annotation[annotation.Length - 1] == FastaDatabase.Delimiter
                                  ? AminoAcid.ProteinCTerm
                                  : AminoAcid.PeptideCTerm;

            var sequence = annotation.Substring(2, annotation.Length - 4);
            return new GraphXSequenceGraph(aaSet, nTerm, sequence, cTerm);
        }

        public GraphXSequenceGraph(AminoAcidSet aminoAcidSet, AminoAcid nTerm, string sequence, AminoAcid cTerm):
                                          base(aminoAcidSet, nTerm, sequence, cTerm)
        {
            DataGraph = new DataGraph();
            BuildGraph();
        }

        private void BuildGraph()
        {
            var sequenceRev = _sequence.Reverse();
            var sequence = sequenceRev.Aggregate("", (current, aa) => current + aa);
            sequence = "\0" + sequence;
            var vertices = new DataVertex[_maxSeqIndex][];
            var mods = AminoAcidSet.GetModificationParams();
            // create vertices
            for (var si = _maxSeqIndex - 2; si > 1; si--)
            {
                var graphSi = si - 1;
                vertices[graphSi] = new DataVertex[_graph[si].Length];
                for (var mi = 0; mi < _graph[si].Length; mi++)
                {
                    var node = _graph[si][mi];
//                    SetSink(node.ModificationCombinationIndex);
                    var mod = mods.GetModificationCombination(node.ModificationCombinationIndex);
                    vertices[graphSi][mi] = new DataVertex
                    {
//                        PrefixComposition = GetComplementaryComposition(si, mi),
                        SuffixComposition = GetComposition(si, mi),
                        ModificationCombination = mod,
                        Text = ""
                    };
                    var vertex = vertices[graphSi][mi];
                    DataGraph.AddVertex(vertex);
                }
            }
            // connect vertices
            for (var si = _maxSeqIndex - 2; si > 2; si--)
            {
                var graphSi = si - 1;
                for (int mi = 0; mi < _graph[si].Length; mi++)
                {
                    var node = _graph[si][mi];
                    var currVertex = vertices[graphSi][mi];
                    foreach (var nextModIndex in node.GetPrevNodeIndices())
                    {
                        var nextVertex = vertices[graphSi - 1][nextModIndex];
                        var currVertexMods = currVertex.ModificationCombination.Modifications;
                        var nextVertexMods = nextVertex.ModificationCombination.Modifications;
                        var result = new List<Modification>(currVertexMods);
                        foreach (var mod in nextVertexMods)
                        {
                            if (result.Contains(mod)) result.Remove(mod);
                        }
                        var edgeModifications = new ModificationCombination(result);
                        var edge = new DataEdge(currVertex, nextVertex)
                        {
                            AminoAcid = AminoAcidSet.GetAminoAcid(sequence[graphSi]),
                            SequenceIndex = graphSi,
                            Modifications = edgeModifications
                        };
                        DataGraph.AddEdge(edge);
                    }
                }
            }
        }
    }
}
