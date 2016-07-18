using System;
using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    using InformedProteomics.Backend.Data.Sequence;

    public class IonTypeFactory
    {

        public IonTypeFactory()
            : this(
                BaseIonType.AllBaseIonTypes,    // a, b, c, x, y, z
                NeutralLoss.CommonNeutralLosses, // H2O and NH3 loss
                3   // up to charge 3
            )
        {
        }

        public IonTypeFactory(int maxCharge)
            : this(
                BaseIonType.AllBaseIonTypes,    // a, b, c, x, y, z
                NeutralLoss.CommonNeutralLosses, // H2O and NH3 loss
                maxCharge   // up to charge 3
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

        public static IonTypeFactory GetDeconvolutedIonTypeFactory(IEnumerable<BaseIonType> baseIons, IEnumerable<NeutralLoss> neutralLosses)
        {
            var ionTypeFactory = new IonTypeFactory(baseIons.Select(ion => ion.GetDeconvolutedIon()), neutralLosses, 1);
            return ionTypeFactory;
        }

        // e.g. y2-NH3
        public IonType GetIonType(string name)
        {
            return _ionTypeMap[name];
        }

        public IonType GetIonType(BaseIonType baseIonType, int charge, NeutralLoss neutralLoss)
        {
            var ionTypeName = string.Format("{0}{1}{2}", baseIonType.Symbol, charge, neutralLoss.Name);
            return _ionTypeMap[ionTypeName];
        }

        public IonType GetIonType(bool isPrefix, int charge, float offset)
        {
            //var intOffset = Convert.ToInt32(offset);
            //if (isPrefix)
            //{
            //    if (charge == 1)
            //    {
            //        if (intOffset == 1) return GetIonType("b");
            //        if (intOffset == -16) return GetIonType("b-NH3");
            //        if (intOffset == -17) return GetIonType("b-H2O");
            //        if (intOffset == -27) return GetIonType("a");
            //    }
            //}
            //else
            //{
            //    if (charge == 1)
            //    {
            //        if (intOffset == 19) return GetIonType("y");
            //        if (intOffset == 2) return GetIonType("y-NH3");
            //        if (intOffset == 1) return GetIonType("y-H2O");
            //        if (intOffset == 3) return GetIonType("z");
            //    }
            //}
            //Console.WriteLine("{0}_{1}_{2} doesn't exist!", isPrefix, charge, offset);
            throw new System.NotImplementedException();
        }

        public IEnumerable<IonType> GetAllKnownIonTypes()
        {
            return _ionTypeMap.Values.ToArray();
        }

        private readonly IEnumerable<BaseIonType> _baseIons;
        private readonly IEnumerable<NeutralLoss> _neutralLosses;
        private readonly int _maxCharge;

        private Dictionary<string, IonType> _ionTypeMap;

        private void GenerateAllKnownIonTypes()
        {
            _ionTypeMap = new Dictionary<string, IonType>();

            for (int charge = 1; charge <= _maxCharge; charge++)
            {
                var chargeStr = charge == 1 ? "" : Convert.ToString(charge);
                foreach (var baseIonType in _baseIons)
                {
                    foreach (var neutralLoss in _neutralLosses)
                    {
                        var name = baseIonType.Symbol + chargeStr + neutralLoss.Name;
                        var offsetComposition = baseIonType.OffsetComposition -
                                                        neutralLoss.Composition;
                        _ionTypeMap[name] = new IonType(name, offsetComposition, charge, baseIonType, neutralLoss);
                    }
                }
            }
        }

    }
}
