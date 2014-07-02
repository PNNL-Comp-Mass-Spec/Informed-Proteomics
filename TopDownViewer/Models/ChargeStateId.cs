using System;
using System.Collections.Generic;

namespace InformedProteomics.TopDownViewer.Models
{
    public class ChargeStateId: IIdData
    {
        public int Charge { get; private set; }
        public List<PrSm> PrSms { get; private set; }

        public ChargeStateId(int charge)
        {
            Charge = charge;
            PrSms = new List<PrSm>();
        }

        public void Add(PrSm data)
        {
            var pos = PrSms.BinarySearch(data);
            if (pos < 0)
            {
                PrSms.Add(data);
                PrSms.Sort();
            }
            else
            {
                throw new InvalidOperationException("Cannot insert duplicate identifications.");
            }
        }

        public PrSm GetHighestScoringPrSm()
        {
            PrSm highest = null;
            foreach (var prsm in PrSms)
            {
                if (highest == null || prsm.MatchedFragments >= highest.MatchedFragments)
                {
                    highest = prsm;
                }
            }
            return highest;
        }
    }

    class ChargeStateComparer : IComparer<ChargeStateId>
    {
        public int Compare(ChargeStateId x, ChargeStateId y)
        {
            return x.Charge.CompareTo(y.Charge);
        }
    }
}
