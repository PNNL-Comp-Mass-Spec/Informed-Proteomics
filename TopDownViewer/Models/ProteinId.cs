using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.TopDownViewer.Models
{
    public class ProteinId: IIdData
    {
        public string Sequence { get; private set; }
        public GraphXSequenceGraph SequenceGraph { get; private set; }
        public List<ProteoformId> Proteoforms { get; private set; }
        public ProteinId(PrSm prsm)
        {
            Sequence = prsm.Protein;
            var aaSet = new AminoAcidSet(prsm.Config.Modifications, prsm.Config.MaxDynamicModificationsPerSequence);
            SequenceGraph = GraphXSequenceGraph.Create(aaSet, prsm.Annotation);
            Proteoforms = new List<ProteoformId>();
            Add(prsm);
        }

        public void Add(PrSm data)
        {
            var searchProteoform = new ProteoformId(data.SequenceText, data.AnnotatedCompositions);
            var pos = Proteoforms.BinarySearch(searchProteoform, new SequenceComparer());
            ProteoformId proteoform;
            if (pos < 0)
            {
                Proteoforms.Add(searchProteoform);
                proteoform = Proteoforms.Last();
                Proteoforms.Sort(new SequenceComparer());
            }
            else
            {
                proteoform = Proteoforms[pos];
            }
            proteoform.Add(data);
        }

        public PrSm GetHighestScoringPrSm()
        {
            PrSm highest = null;
            foreach (var proteoform in Proteoforms)
            {
                var pfHighest = proteoform.GetHighestScoringPrSm();
                if (highest == null || pfHighest.MatchedFragments >= highest.MatchedFragments)
                {
                    highest = pfHighest;
                }
            }
            return highest;
        }

        public List<PrSm> Scans
        {
            get
            {
                return (from proteoform in Proteoforms
                        from charge in proteoform.ChargeStates
                        from prsm in charge.PrSms select prsm).ToList();
            }
        }

        public List<ModifiedCompositionList> AnnotatedCompositions
        {
            get
            {
                var graph = new List<ModifiedCompositionList>();
                var mods = SequenceGraph.GetModificationCombinations();
                for (int i = 0; i < mods.Count(); i++)
                {
                    var fragComps = SequenceGraph.GetFragmentCompositions(i, -1).Reverse().ToArray();
                    graph.Add(new ModifiedCompositionList(Sequence, mods[i], fragComps));
                }
                return graph;
            }
        }
    }

    public class ModifiedCompositionList : List<Tuple<char, Composition>>
    {
        public ModificationCombination ModificationCombinations { get; private set; }
        public ModifiedCompositionList(string sequence,
                                       ModificationCombination modificationCombination,
                                       IList<Composition> compositions)
        {
            ModificationCombinations = modificationCombination;
            var revSequence = sequence.Reverse().ToArray();
            for (int i = 0; i < revSequence.Length; i++)
            {
                Add(new Tuple<char, Composition>(revSequence[i], compositions[i]));
            }
        }

        public ModificationCombination GetModificationCombinationDifference(ModificationCombination modificationCombination)
        {
            return new ModificationCombination(modificationCombination.Modifications.Except(ModificationCombinations.Modifications).ToList());
        }
    }

}
