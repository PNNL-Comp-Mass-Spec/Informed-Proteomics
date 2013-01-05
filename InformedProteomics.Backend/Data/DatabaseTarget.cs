using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using DeconTools.Backend.Core;
using DeconTools.Backend.Utilities.IsotopeDistributionCalculation;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data
{
	public class DatabaseTarget
	{
		private Sequence _sequence;
		private double _minMz;
		private double _maxMz;
		private short _minChargeState;
		private short _maxChargeState;

		public DatabaseTarget(Sequence sequence, short minChargeState, short maxChargeState)
			: this(sequence, double.MinValue, double.MaxValue, minChargeState, maxChargeState)
		{
		}

		public DatabaseTarget(Sequence sequence, double minMz, double maxMz, short minChargeState, short maxChargeState)
		{
			this._sequence = sequence;
			this._minMz = minMz;
			this._maxMz = maxMz;
			this._minChargeState = minChargeState;
			this._maxChargeState = maxChargeState;
		}

		public List<TargetBase> CreatePrecursorTargets()
		{
			List<TargetBase> targetList = new List<TargetBase>();

			string sequenceString = this._sequence.SequenceString;
			string empiricalFormula = this._sequence.Composition.ToString();
			double monoIsotopicMass = this._sequence.GetMass();
			float normalizedElutionTime = PeptideUtil.CalculatePredictedElutionTime(sequenceString);

			// Create precursor targets for each charge state
			for (short chargeState = _minChargeState; chargeState <= _maxChargeState; chargeState++)
			{
				double mz = (monoIsotopicMass / chargeState) + Globals.PROTON_MASS;

				// Move on to the next charge state if the m/z value is too large
				if (mz > _maxMz) continue;

				// Stop creating targets if the m/z is too small (already created all possible targets)
				if (mz < _minMz) break;

				DatabaseSubTarget target = new DatabaseSubTarget(sequenceString, empiricalFormula, monoIsotopicMass, mz, chargeState, normalizedElutionTime, 1);
				targetList.Add(target);
			}

			// TODO: Create fragment targets. Create fragment targets based on precursor results. Possibly would need a seperate method for creating this after initial targets were searched.

			return targetList;
		}

		public List<TargetBase> CreateFragmentTargets(int maxChargeState)
		{
			List<TargetBase> targetList = new List<TargetBase>();

			string sequenceString = this._sequence.SequenceString;
			float normalizedElutionTime = PeptideUtil.CalculatePredictedElutionTime(sequenceString);

			IEnumerable<Fragment> fragmentList = FragmentationUtil.FindFragmentsForPeptide(sequenceString, maxChargeState);

			foreach (Fragment fragment in fragmentList)
			{
				string fragmentEmpFormula = IsotopicDistributionCalculator.Instance.GetAveragineFormulaAsString(fragment.Mass);

				DatabaseSubTarget target = new DatabaseSubTarget(sequenceString, fragmentEmpFormula, fragment.Mass, fragment.Mz, (short)fragment.ChargeState, normalizedElutionTime, 2);
				targetList.Add(target);
			}

			return targetList;
		} 
	}
}
