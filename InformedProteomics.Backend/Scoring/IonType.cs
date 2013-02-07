using System;

namespace InformedProteomics.Backend.Scoring
{
    public class IonType
    {
        public string Type { get; private set; }
        public string NeutralLoss { get; private set; }
        public bool IsPrefix { get; private set; }
        public int Charge { get; private set; }

        public IonType(string symbol, int charge)
        {
            Type = symbol.Substring(0, 1);
            IsPrefix = Type.ToCharArray()[0] < 'g';
            var i = Math.Max(symbol.LastIndexOf('+'), symbol.LastIndexOf('-'));
            NeutralLoss = "";
            if (i > 0)
            {
                NeutralLoss += symbol.Substring(i);
            }
            
            Charge = charge;
        }

        public override int GetHashCode()
        {
            return Type.GetHashCode()*NeutralLoss.GetHashCode()*Charge.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            if (this == obj) return true;
            var other = (IonType) obj;
            return Type.Equals(other.Type) && NeutralLoss.Equals(other.NeutralLoss) && Charge == other.Charge;

        }

        public override string ToString()
        {
            return Type + "\t" + NeutralLoss + "\t" + Charge;
        }
    }
}
