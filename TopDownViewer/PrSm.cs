using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.TopDownViewer
{
    public class PrSm
    {
        public string Sequence { get; private set; }
        public PrSm(string sequence)
        {
            Sequence = sequence;
            _chargeStates = new Dictionary<int, ChargeStateData>();
        }

        public List<int> Charges
        {
            get { return _chargeStates.Keys.ToList(); }
        }

        public ChargeStateData GetCharge(int charge)
        {
            return _chargeStates[charge];
        }

        public void AddCharge(ChargeStateData data)
        {
            if (!_chargeStates.ContainsKey(data.Charge))    _chargeStates.Add(data.Charge, data);
        }

        private readonly Dictionary<int, ChargeStateData> _chargeStates;
    }
}
