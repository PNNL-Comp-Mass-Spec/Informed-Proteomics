using System;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.Data.Targets
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
