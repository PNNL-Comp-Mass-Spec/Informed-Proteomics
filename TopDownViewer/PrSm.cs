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
            _chargeStates = new Dictionary<int, Dictionary<int, IdData>>();
        }

        public List<int> Charges
        {
            get { return _chargeStates.Keys.ToList(); }
        }

        public List<int> GetScanNums(int charge)
        {
            return _chargeStates[charge].Keys.ToList();
        }

        public IdData GetData(int charge, int scanNum)
        {
            return _chargeStates[charge][scanNum];
        }

        public void AddId(IdData data)
        {

            var charge = data.Charge;
            var scanNum = data.Scan;
            if (!_chargeStates.ContainsKey(charge))
                _chargeStates.Add(charge, new Dictionary<int, IdData>());

            _chargeStates[charge].Add(scanNum, data);
        }

        public override string ToString()
        {
            return Sequence;
        }

        private readonly Dictionary<int, Dictionary<int, IdData>> _chargeStates;
    }
}
