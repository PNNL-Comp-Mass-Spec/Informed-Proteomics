using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using DeconTools.Backend.Core;
using DeconTools.Backend.ProcessingTasks.ChromatogramProcessing;
using DeconTools.Workflows.Backend.Core;
using InformedProteomics.Backend.Data;
using InformedProteomics.Backend.Data.Targets;
using InformedProteomics.Backend.IMS;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test
{
	[TestFixture]
	public class TestDeconToolsWorkflow
	{
		[Test]
		public void TestUIMFTargetedMSMSWorkflowSingleTarget()
		{
			string datasetPath = @"\\protoapps\UserData\Slysz\Standard_Testing\Targeted_FeatureFinding\UIMF_Targeted_MSMS_Testing\RawData\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo.UIMF";
			string executorParameterFilename = @"\\protoapps\UserData\Slysz\Standard_Testing\Targeted_FeatureFinding\UIMF_Targeted_MSMS_Testing\Parameters\UIMFTargetedMSMSWorkflowExecutorParameters.xml";

			var executorParameters = new BasicTargetedWorkflowExecutorParameters();
			executorParameters.LoadParameters(executorParameterFilename);

			var executor = new BasicTargetedWorkflowExecutor(executorParameters, datasetPath);
			executor.Execute();
		}

		[Test]
		public void TestReadPeptideAndRunTargetedWorkflow()
		{
			short minChargeState = 1;
			short maxChargeState = 5;
			double minMz = 200;
			double maxMz = 2500;

			string datasetPath = @"\\protoapps\UserData\Slysz\Standard_Testing\Targeted_FeatureFinding\UIMF_Targeted_MSMS_Testing\RawData\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo.UIMF";

			var executorParameters = CreateWorkflowExecutorParameters();
			var workflowParameters = CreateWorkflowParameters();

			var executor = new UIMFTargetedWorkflowExecutor(executorParameters, workflowParameters, datasetPath);

			executor.Targets = new TargetCollection();

			AminoAcidSet aminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);

			// TODO: Read in a peptide list?
			List<string> peptideSequenceList = new List<string> { "IFFHLNAVALGDGGHYTCR" };
			foreach (string peptideSequence in peptideSequenceList)
			{
				IEnumerable<Composition> compositions = aminoAcidSet.GetCompositions(peptideSequence);

				foreach (Composition composition in compositions)
				{
					Composition peptideComposition = composition + Composition.H2O;

					// Create the database target
					Sequence sequence = new Sequence(peptideComposition, peptideSequence);

					DatabaseTarget databaseTarget = new DatabaseTarget(sequence, minMz, maxMz, minChargeState, maxChargeState);

					// Get the list of precursor targets from the DatabaseTarget object and pass them in
					executor.Targets.TargetList = databaseTarget.CreatePrecursorTargets();

					foreach (var target in executor.Targets.TargetList)
					{
						Console.WriteLine(target.ChargeState + "\t" + target.EmpiricalFormula + "\t" + target.Code + "\t" + target.MonoIsotopicMass + "\t" + target.MZ + "\t" + target.NormalizedElutionTime);
					}

					// Execute the workflow with the new targets
					executor.Execute();

					// TODO: Save results? Print results?
					List<TargetedResultBase> resultList = executor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
					foreach (var targetedResultBase in resultList)
					{
						TargetBase target = targetedResultBase.Target;

						Console.WriteLine("************************************************");
						Console.WriteLine("CS = " + target.ChargeState + "\tMass = " + target.MonoIsotopicMass + "\tm/z = " + target.MZ);

						if (targetedResultBase.ErrorDescription != null && targetedResultBase.ErrorDescription.Any())
						{
							Console.WriteLine(targetedResultBase.ErrorDescription);
						}
						if(targetedResultBase.FailedResult)
						{
							Console.WriteLine(targetedResultBase.FailureType);
						}
						else
						{
							foreach (var peakQualityData in targetedResultBase.ChromPeakQualityList)
							{
								// We only care about scoring peaks that have a possibility of representing the target
								if (peakQualityData.IsotopicProfileFound)
								{
									Console.WriteLine("NET = " + peakQualityData.Peak.NETValue + "\tScan = " + peakQualityData.Peak.XValue + "\tFit = " + peakQualityData.FitScore + "\tMass = " + peakQualityData.IsotopicProfile.MonoIsotopicMass);
								}
								else
								{
									Console.WriteLine("Isotopic Profile Not Found");
								}
							}
						}
					}
				}
			}
		}

		[Test]
		public void TestTargetFragments()
		{
			short minChargeState = 1;
			short maxChargeState = 5;
			double minMz = 200;
			double maxMz = 2500;

			string datasetPath = @"\\protoapps\UserData\Slysz\Standard_Testing\Targeted_FeatureFinding\UIMF_Targeted_MSMS_Testing\RawData\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo.UIMF";

			var executorParameters = CreateWorkflowExecutorParameters();
			var workflowParameters = CreateWorkflowParameters();

			var executor = new UIMFTargetedWorkflowExecutor(executorParameters, workflowParameters, datasetPath);

			executor.Targets = new TargetCollection();

			AminoAcidSet aminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);

			// TODO: Read in a peptide list?
			List<string> peptideSequenceList = new List<string> { "EYANQFMWEYSTNYGQAPLSLLVSYTK" };
			foreach (string peptideSequence in peptideSequenceList)
			{
				IEnumerable<Composition> compositions = aminoAcidSet.GetCompositions(peptideSequence);

				foreach (Composition composition in compositions)
				{
					Composition peptideComposition = composition + Composition.H2O;

					// Create the database target
					Sequence sequence = new Sequence(peptideComposition, peptideSequence);

					DatabaseTarget databaseTarget = new DatabaseTarget(sequence, minMz, maxMz, minChargeState, maxChargeState);

					// Get the list of precursor targets from the DatabaseTarget object and pass them in
					executor.Targets.TargetList = databaseTarget.CreateFragmentTargets(3);

					//foreach (var target in executor.Targets.TargetList)
					//{
					//    Console.WriteLine(target.ChargeState + "\t" + target.EmpiricalFormula + "\t" + target.Code + "\t" + target.MonoIsotopicMass + "\t" + target.MZ + "\t" + target.NormalizedElutionTime);
					//}

					// Execute the workflow with the new targets
					executor.Execute();

					// TODO: Save results? Print results?
					List<TargetedResultBase> resultList = executor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
					foreach (var targetedResultBase in resultList)
					{
						TargetBase target = targetedResultBase.Target;

						Console.WriteLine("************************************************");
						Console.WriteLine("CS = " + target.ChargeState + "\tMass = " + target.MonoIsotopicMass + "\tm/z = " + target.MZ);

						if (targetedResultBase.ErrorDescription != null && targetedResultBase.ErrorDescription.Any())
						{
							Console.WriteLine(targetedResultBase.ErrorDescription);
						}
						if (targetedResultBase.FailedResult)
						{
							Console.WriteLine(targetedResultBase.FailureType);
						}
						else
						{
							foreach (var peakQualityData in targetedResultBase.ChromPeakQualityList)
							{
								// We only care about scoring peaks that have a possibility of representing the target
								if (peakQualityData.IsotopicProfileFound)
								{
									Console.WriteLine("NET = " + peakQualityData.Peak.NETValue + "\tScan = " + peakQualityData.Peak.XValue + "\tFit = " + peakQualityData.FitScore + "\tMass = " + peakQualityData.IsotopicProfile.MonoIsotopicMass);
								}
								else
								{
									Console.WriteLine("Isotopic Profile Not Found");
								}
							}
						}
					}
				}
			}
		}

		private BasicTargetedWorkflowExecutorParameters CreateWorkflowExecutorParameters()
		{
			string currentDirectory = Directory.GetCurrentDirectory();

			var executorParameters = new BasicTargetedWorkflowExecutorParameters();

			executorParameters.CopyRawFileLocal = false;
			executorParameters.DeleteLocalDatasetAfterProcessing = false;
			executorParameters.LoggingFolder = currentDirectory + "/logs";
			executorParameters.ResultsFolder = currentDirectory + "/results";
			executorParameters.TargetedAlignmentIsPerformed = false;

			return executorParameters;
		}

		private UIMFTargetedMSMSWorkflowCollapseIMSParameters CreateWorkflowParameters()
		{
			var workflowParameters = new UIMFTargetedMSMSWorkflowCollapseIMSParameters();

			workflowParameters.AreaOfPeakToSumInDynamicSumming = 2;
			workflowParameters.ChromatogramCorrelationIsPerformed = false;
			workflowParameters.ChromGeneratorMode = Globals.ChromatogramGeneratorMode.MOST_ABUNDANT_PEAK;
			workflowParameters.ChromGenSourceDataPeakBR = 2;
			workflowParameters.ChromGenSourceDataSigNoise = 3;
			workflowParameters.ChromNETTolerance = 0.2;
			workflowParameters.ChromPeakDetectorPeakBR = 1;
			workflowParameters.ChromPeakDetectorSigNoise = 1;
			workflowParameters.ChromPeakSelectorMode = Globals.PeakSelectorMode.SmartUIMF;
			workflowParameters.ChromSmootherNumPointsInSmooth = 9;
			workflowParameters.ChromToleranceInPPM = 25;
			workflowParameters.MaxScansSummedInDynamicSumming = 100;
			workflowParameters.MSPeakDetectorPeakBR = 1.3;
			workflowParameters.MSPeakDetectorSigNoise = 3;
			workflowParameters.MSToleranceInPPM = 25;
			workflowParameters.MultipleHighQualityMatchesAreAllowed = true;
			workflowParameters.NumMSScansToSum = 1;
			workflowParameters.NumChromPeaksAllowedDuringSelection = int.MaxValue;
			workflowParameters.ProcessMsMs = true;
			workflowParameters.ResultType = Globals.ResultType.BASIC_TARGETED_RESULT;
			workflowParameters.SummingMode = SummingModeEnum.SUMMINGMODE_STATIC;

			return workflowParameters;
		}
	}
}
