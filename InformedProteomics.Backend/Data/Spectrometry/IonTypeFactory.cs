using System;
using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    /// <summary>
    /// Factory for creating Ion types
    /// </summary>
    public class IonTypeFactory
    {
        private bool removeReduntantIonTypes;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="removeRedundantIonTypes"></param>
        public IonTypeFactory(bool removeRedundantIonTypes = true)
            : this(
                BaseIonType.AllBaseIonTypes,    // a, b, c, x, y, z
                NeutralLoss.CommonNeutralLosses, // H2O and NH3 loss
                3   // up to charge 3
            )
        {
            this.removeReduntantIonTypes = removeRedundantIonTypes;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="maxCharge"></param>
        public IonTypeFactory(int maxCharge)
            : this(
                BaseIonType.AllBaseIonTypes,    // a, b, c, x, y, z
                NeutralLoss.CommonNeutralLosses, // H2O and NH3 loss
                maxCharge   // up to charge 3
            )
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="baseIons"></param>
        /// <param name="neutralLosses"></param>
        /// <param name="maxCharge"></param>
        /// <param name="removeRedundantIonTypes"></param>
        public IonTypeFactory(
            IEnumerable<BaseIonType> baseIons,
            IEnumerable<NeutralLoss> neutralLosses,
            int maxCharge,
            bool removeRedundantIonTypes = true)
        {
            _baseIons = baseIons;
            _neutralLosses = neutralLosses;
            _maxCharge = maxCharge;
            this.removeReduntantIonTypes = removeRedundantIonTypes;
            GenerateAllKnownIonTypes();
        }

        /// <summary>
        /// Get a <see cref="IonTypeFactory"/> for the deconvoluted versions of <paramref name="baseIons"/> and <paramref name="neutralLosses"/>
        /// </summary>
        /// <param name="baseIons"></param>
        /// <param name="neutralLosses"></param>
        /// <returns></returns>
        public static IonTypeFactory GetDeconvolutedIonTypeFactory(IEnumerable<BaseIonType> baseIons, IEnumerable<NeutralLoss> neutralLosses)
        {
            var ionTypeFactory = new IonTypeFactory(baseIons.Select(ion => ion.GetDeconvolutedIon()), neutralLosses, 1);
            return ionTypeFactory;
        }

        /// <summary>
        /// Get an ion type by name, e.g. y2-NH3
        /// </summary>
        /// <param name="name"></param>
        /// <returns></returns>
        public IonType GetIonType(string name)
        {
            return _ionTypeMap[name];
        }

        /// <summary>
        /// Get ion type according to the parameters
        /// </summary>
        /// <param name="baseIonType"></param>
        /// <param name="charge"></param>
        /// <param name="neutralLoss"></param>
        /// <returns></returns>
        public IonType GetIonType(BaseIonType baseIonType, int charge, NeutralLoss neutralLoss)
        {
            var ionTypeName = string.Format("{0}{1}{2}", baseIonType.Symbol, charge, neutralLoss.Name);
            return _ionTypeMap[ionTypeName];
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

        /// <summary>
        /// Get all known ion types
        /// </summary>
        /// <returns></returns>
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
                            (baseIonType == BaseIonType.C && neutralLoss == NeutralLoss.NH3) ||
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
