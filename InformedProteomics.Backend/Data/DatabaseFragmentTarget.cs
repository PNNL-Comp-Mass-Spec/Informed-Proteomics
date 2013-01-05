using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data
{
	public class DatabaseFragmentTarget : DatabaseSubTarget
	{
		public Fragment Fragment { get; set; }

		public DatabaseFragmentTarget(Fragment fragment, String sequence, String empiricalFormula, float elutionTime) : base(sequence, empiricalFormula, fragment.Mass, fragment.Mz, (short)fragment.ChargeState, elutionTime)
		{
			this.Fragment = fragment;
			this.MsLevel = 2;
		}
	}
}
