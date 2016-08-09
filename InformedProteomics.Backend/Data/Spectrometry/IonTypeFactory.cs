using System;
using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    using InformedProteomics.Backend.Data.Sequence;

    public class IonTypeFactory
    {

        private bool removeReduntantIonTypes;

        public IonTypeFactory(bool removeRedundantIonTypes = true)
            : this(
                BaseIonType.AllBaseIonTypes,    // a, b, c, x, y, z
                NeutralLoss.CommonNeutralLosses, // H2O and NH3 loss
                3   // up to charge 3
            )
        {
            this.removeReduntantIonTypes = removeRedundantIonTypes;
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

        /// <summary>
        /// Given a list of IonTypes of charge 1 and a charge value, a list of the same IonTypes of charge up to and
        /// including charge are outputted.
        /// </summary>
        /// <param name="ionTypes">List of IonTypes of charge 1</param>
        /// <param name="charge">IonTypes of charge up to and including charge will be outputted</param>
        /// <returns>List of charged IonTypes</returns>
        public static List<IonType> GetIonTypesFromDecharged(IEnumerable<IonType> ionTypes, int charge)
        {
            List<IonType> chargedIonTypes = new List<IonType>();
            foreach (IonType ionType in ionTypes)
            {
                for (int ch = 1; ch <= charge; ch++)
                {
                    IonType chargedIonType = 
                        new IonType(ionType.Name, ionType.OffsetComposition, ch, ionType.BaseIonType, ionType.NeutralLoss);
                    chargedIonTypes.Add(chargedIonType);
                }
            }
            return chargedIonTypes;
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
                        if (this.removeReduntantIonTypes &&
                            ((baseIonType == BaseIonType.Ar && neutralLoss == NeutralLoss.H2O) ||
                            (baseIonType == BaseIonType.Xr && neutralLoss == NeutralLoss.H2O) ||
                            (baseIonType == BaseIonType.YM1 && neutralLoss == NeutralLoss.NH3) ||
                            ((baseIonType == BaseIonType.D || baseIonType == BaseIonType.W ||
                            baseIonType == BaseIonType.V) && neutralLoss != NeutralLoss.NoLoss)))
                        {
                            continue;
                        }

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
