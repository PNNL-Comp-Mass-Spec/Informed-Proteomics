using System;
using System.Collections.Generic;
using System.Linq;
using DeconTools.Backend;
using InformedProteomics.Backend.Data;
using MwtWinDll;

namespace InformedProteomics.Backend.Utils
{
	public class FragmentationUtil
	{
		private static MolecularWeightCalculator m_molecularWeightCalculator;
		private static MWPeptideClass.udtFragmentationSpectrumOptionsType m_fragmentationOptions;

		static FragmentationUtil()
		{
			m_molecularWeightCalculator = new MolecularWeightCalculator();
			m_molecularWeightCalculator.SetElementMode(MWElementAndMassRoutines.emElementModeConstants.emIsotopicMass);

			m_fragmentationOptions = new MWPeptideClass.udtFragmentationSpectrumOptionsType();
			m_fragmentationOptions.Initialize();

			MWPeptideClass.udtIonTypeOptionsType[] ionTypeOptions = new MWPeptideClass.udtIonTypeOptionsType[5];
			ionTypeOptions[(int)MWPeptideClass.itIonTypeConstants.itAIon] = new MWPeptideClass.udtIonTypeOptionsType { ShowIon = true, NeutralLossWater = false, NeutralLossAmmonia = false };
			ionTypeOptions[(int)MWPeptideClass.itIonTypeConstants.itBIon] = new MWPeptideClass.udtIonTypeOptionsType { ShowIon = true, NeutralLossWater = false, NeutralLossAmmonia = false };
			ionTypeOptions[(int)MWPeptideClass.itIonTypeConstants.itYIon] = new MWPeptideClass.udtIonTypeOptionsType { ShowIon = true, NeutralLossWater = false, NeutralLossAmmonia = false };
			//ionTypeOptions[(int)MWPeptideClass.itIonTypeConstants.itCIon] = new MWPeptideClass.udtIonTypeOptionsType { ShowIon = true };
			//ionTypeOptions[(int)MWPeptideClass.itIonTypeConstants.itZIon] = new MWPeptideClass.udtIonTypeOptionsType { ShowIon = true };

			m_fragmentationOptions.IonTypeOptions = ionTypeOptions;
			m_fragmentationOptions.DoubleChargeIonsShow = true;
			m_fragmentationOptions.TripleChargeIonsShow = true;
			m_fragmentationOptions.IntensityOptions.BYIonShoulder = 0;
			//m_fragmentationOptions.IntensityOptions.IonType[(int) MWPeptideClass.itIonTypeConstants.itAIon] = 100;
			m_fragmentationOptions.IntensityOptions.IonType[(int) MWPeptideClass.itIonTypeConstants.itBIon] = 100;
			m_fragmentationOptions.IntensityOptions.IonType[(int) MWPeptideClass.itIonTypeConstants.itYIon] = 100;
			//fragmentationOptions.IntensityOptions.IonType[(int)MWPeptideClass.itIonTypeConstants.itCIon] = 100;
			//fragmentationOptions.IntensityOptions.IonType[(int)MWPeptideClass.itIonTypeConstants.itZIon] = 100;
		}

		public static List<Fragment> FindFragmentsForPeptide(String peptideSequence, int maxChargeState)
		{
			List<Fragment> theoreticalFragments = new List<Fragment>();

			m_molecularWeightCalculator.Peptide.SetSequence(peptideSequence, MWPeptideClass.ntgNTerminusGroupConstants.ntgHydrogen, MWPeptideClass.ctgCTerminusGroupConstants.ctgHydroxyl, false);
			m_molecularWeightCalculator.Peptide.SetFragmentationSpectrumOptions(ref m_fragmentationOptions);

			MWPeptideClass.udtFragmentationSpectrumDataType[] fragmentationDataArray = new MWPeptideClass.udtFragmentationSpectrumDataType[0];
			int ionCount = m_molecularWeightCalculator.Peptide.GetFragmentationMasses(ref fragmentationDataArray);

			for (int i = 0; i < ionCount; i++)
			{
				MWPeptideClass.udtFragmentationSpectrumDataType fragmentationData = fragmentationDataArray[i];

				// Yes, the Molecular Weight Calculator returns the m/z value but calls it Mass
				if (fragmentationData.Charge <= maxChargeState)
				{
					Fragment theoreticalFragment = new Fragment
					                               	{
					                               		ChargeState = fragmentationData.Charge,
					                               		Mz = fragmentationData.Mass,
														Mass = (fragmentationData.Mass - Globals.PROTON_MASS) * fragmentationData.Charge,
					                               		IonType = fragmentationData.IonType.ToString().ElementAt(2).ToString().ToLower(),
					                               		ResidueNumber = fragmentationData.SourceResidueNumber,
					                               		IonSymbol = fragmentationData.Symbol
					                               	};
					theoreticalFragments.Add(theoreticalFragment);
				}
			}

			return theoreticalFragments;
		}
	}
}
