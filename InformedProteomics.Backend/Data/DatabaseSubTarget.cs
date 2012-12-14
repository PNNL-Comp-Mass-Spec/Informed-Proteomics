using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using DeconTools.Backend.Core;

namespace InformedProteomics.Backend.Data
{
	public class DatabaseSubTarget : TargetBase
	{
		public DatabaseSubTarget(String sequence, String empiricalFormula, double mass, double mz, short chargeState, float elutionTime) : base()
		{
			this.Code = sequence;
			this.EmpiricalFormula = empiricalFormula;
			this.MonoIsotopicMass = mass;
			this.MZ = mz;
			this.ChargeState = chargeState;
			ElutionTimeUnit = Globals.ElutionTimeUnit.NormalizedElutionTime;
			this.NormalizedElutionTime = elutionTime;
		}

		public override string GetEmpiricalFormulaFromTargetCode()
		{
			return this.EmpiricalFormula;
		}
	}
}
