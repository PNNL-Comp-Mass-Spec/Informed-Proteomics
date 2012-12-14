using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using DeconTools.Backend.Core;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.Data
{
	public class DatabaseTarget
	{
		private Sequence _sequence;
		private List<TargetBase> _targetList;
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

			_targetList = CreateTargets();
		}

		private List<TargetBase> CreateTargets()
		{
			List<TargetBase> targetList = new List<TargetBase>();

			string sequenceString = this._sequence.SequenceStr;
			string empiricalFormula = this._sequence.GetComposition().ToString();
			double monoIsotopicMass = this._sequence.GetMass();
			float normalizedElutionTime = PeptideUtil.CalculatePredictedElutionTime(sequenceString);

			// Create precursor targets for each charge state
			for (short chargeState = _minChargeState; chargeState <= _maxChargeState; chargeState++)
			{
				double mz = (monoIsotopicMass / chargeState) + Globals.PROTON_MASS;

				// Move on to the next charge state if the m/z value is too small
				if (mz < _minMz) continue;

				// Stop creating targets if the m/z is too large (already created all possible targets)
				if (mz > _maxMz) break;

				DatabaseSubTarget target = new DatabaseSubTarget(sequenceString, empiricalFormula, monoIsotopicMass, mz, chargeState, normalizedElutionTime);
				targetList.Add(target);
			}

			// TODO: Create fragment targets. Create fragment targets based on precursor results. Possibly would need a seperate method for creating this after initial targets were searched.

			return targetList;
		}

	}
}
