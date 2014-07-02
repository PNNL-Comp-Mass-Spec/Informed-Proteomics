using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Composition;

namespace InformedProteomics.TopDownViewer.Models
{
    public class ProteoformId: IIdData
    {
        public string SequenceText { get; private set; }
        public List<Tuple<string,Composition>> AnnotatedCompositions { get; private set; }
        public List<ChargeStateId> ChargeStates { get; private set; }

        public ProteoformId(string sequenceText, List<Tuple<string, Composition>> annotatedCompositions)
        {
            SequenceText = sequenceText;
            AnnotatedCompositions = annotatedCompositions;
            ChargeStates = new List<ChargeStateId>();
        }

        public void Add(PrSm data)
        {
            var searchChargeState = new ChargeStateId(data.Charge);
            var pos = ChargeStates.BinarySearch(searchChargeState, new ChargeStateComparer());
            ChargeStateId chargeState;
            if (pos < 0)
            {
                ChargeStates.Add(searchChargeState);
                chargeState = ChargeStates.Last();
                ChargeStates.Sort(new ChargeStateComparer());
            }
            else
            {
                chargeState = ChargeStates[pos];
            }
            chargeState.Add(data);
        }

        public PrSm GetHighestScoringPrSm()
        {
            PrSm highest = null;
            foreach (var chargeState in ChargeStates)
            {
                var chargeStateHighest = chargeState.GetHighestScoringPrSm();
                if (highest == null || chargeStateHighest.MatchedFragments >= highest.MatchedFragments)
                {
                    highest = chargeStateHighest;
                }
            }
            return highest;
        }
    }

    internal class SequenceComparer : IComparer<ProteoformId>
    {
        public int Compare(ProteoformId x, ProteoformId y)
        {
            return (String.Compare(x.SequenceText, y.SequenceText, StringComparison.Ordinal));
        }
    }
}
