using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class IonTypeFactory
    {
        private readonly IEnumerable<BaseIonType> _baseIons;
        private readonly IEnumerable<NeutralLoss> _neutralLosses;
        private readonly int _maxCharge;

        private Dictionary<string, IonType> _ionTypeMap;

        public IonTypeFactory()
            : this(
                BaseIonType.AllBaseIonTypes,    // a, b, c, x, y, z
                NeutralLoss.CommonNeutralLosses, // H2O and NH3 loss
                3   // up to charge 3
            )
        {
        }

        public IonTypeFactory(IEnumerable<BaseIonType> baseIons, IEnumerable<NeutralLoss> neutralLosses, int maxCharge)
        {
            _baseIons = baseIons;
            _neutralLosses = neutralLosses;
            _maxCharge = maxCharge;
            GenerateAllKnownIonTypes();
        }

        public IonType GetIonType(string name)
        {
            return _ionTypeMap[name];
        }

        public IEnumerable<IonType> GetAllKnownIonTypes()
        {
            return _ionTypeMap.Values.ToArray();
        }

        private void GenerateAllKnownIonTypes()
        {
            _ionTypeMap = new Dictionary<string, IonType>();

            for (int charge = 1; charge <= _maxCharge; charge++)
            {
                var chargeStr = charge == 1 ? "" : Convert.ToString(charge);
                foreach (var baseIonType in _baseIons)
                {
                    foreach (var neutralLossComposition in _neutralLosses)
                    {
                        string name = baseIonType.Symbol + chargeStr + neutralLossComposition.Name;
                        Composition offsetComposition = baseIonType.OffsetComposition -
                                                        neutralLossComposition.Composition;
                        _ionTypeMap[name] = new IonType(name, offsetComposition, charge, baseIonType.IsPrefix);
                    }
                }
            }
        }
    }
}
