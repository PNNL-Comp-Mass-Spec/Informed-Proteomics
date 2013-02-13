using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend.Core;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Targets
{
    class NominalMassTarget : TargetBase
    {
        public NominalMassTarget(int nominalMass, short chargeState) : base()
        {
            this.MZ = nominalMass/Constants.RescalingConstant;
            this.ChargeState = chargeState;
        }

        public override string GetEmpiricalFormulaFromTargetCode()
        {
            throw new NotImplementedException();
        }

        public override string ToString()
        {
            return string.Format("ID: {0}, MZ: {2}, ChargeState: {3}", ID, MZ, ChargeState);
        }
    }
}
