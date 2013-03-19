using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend.Core;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Targets
{
    public class NominalMassTarget : TargetBase
    {
        public NominalMassTarget(int nominalMass, short chargeState, int msLevel) : base()
        {
            NominalMass = nominalMass;
            this.MZ = nominalMass/Constants.RescalingConstant;
            this.ChargeState = chargeState;
            this.MsLevel = msLevel;
        }

        public int NominalMass { get; private set; }

        public override string GetEmpiricalFormulaFromTargetCode()
        {
            throw new NotImplementedException();
        }

        public override string ToString()
        {
            return string.Format("ID: {0}, MZ: {1}, ChargeState: {2}", ID, MZ, ChargeState);
        }
    }
}
