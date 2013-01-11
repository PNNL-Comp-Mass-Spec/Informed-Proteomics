using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using DeconTools.Backend.Algorithms;
using DeconTools.Backend.Core;
using DeconTools.Backend.ProcessingTasks.ChromatogramProcessing;
using DeconTools.Workflows.Backend.Core;
using InformedProteomics.Backend.Data;
using InformedProteomics.Backend.Data.Results;
using InformedProteomics.Backend.Data.Sequence;
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
		public void TestTargetPeptideAndFragments()
		{
			short minChargeState = 1;
			short maxChargeState = 5;
			double minMz = 200;
			double maxMz = 2500;

			string datasetPath = @"\\protoapps\UserData\Slysz\Standard_Testing\Targeted_FeatureFinding\UIMF_Targeted_MSMS_Testing\RawData\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo.UIMF";

			var executorParameters = CreateWorkflowExecutorParameters();
			var workflowParameters = CreateWorkflowParameters();

			// Create workflow for precursors
			var precursorExecutor = new UIMFTargetedWorkflowExecutor(executorParameters, workflowParameters, datasetPath);
			precursorExecutor.Targets = new TargetCollection();

			// Create workflow for fragments
			var fragmentWorkflowParameters = CreateWorkflowParameters();
			fragmentWorkflowParameters.ChromNETTolerance = 0.005;
			var fragmentExecutor = new UIMFTargetedWorkflowExecutor(executorParameters, fragmentWorkflowParameters, datasetPath);
			fragmentExecutor.Targets = new TargetCollection();

			double slope;
			double intercept;
			double rSquared;

			AminoAcidSet aminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);

			// TODO: Read in a peptide list?
			List<string> peptideSequenceList = new List<string> { "EYANQFMWEYSTNYGQAPLSLLVSYTK" };
			foreach (string peptideSequence in peptideSequenceList)
			{
				IEnumerable<Composition> compositions = aminoAcidSet.GetCompositions(peptideSequence);

				foreach (Composition composition in compositions)
				{
					List<DatabaseMultipleSubTargetResult> resultList = new List<DatabaseMultipleSubTargetResult>();
					Composition peptideComposition = composition + Composition.H2O;

					// Create the database target
					Sequence sequence = new Sequence(peptideComposition, peptideSequence);

					DatabaseTarget databaseTarget = new DatabaseTarget(sequence, minMz, maxMz, minChargeState, maxChargeState);

					// Get the list of precursor targets from the DatabaseTarget object and pass them in
					precursorExecutor.Targets.TargetList = databaseTarget.CreatePrecursorTargets();

					foreach (var target in precursorExecutor.Targets.TargetList)
					{
						Console.WriteLine(target.ChargeState + "\t" + target.EmpiricalFormula + "\t" + target.Code + "\t" +
										  target.MonoIsotopicMass + "\t" + target.MZ + "\t" + target.NormalizedElutionTime);
					}

					// Execute the workflow with the new targets
					precursorExecutor.Execute();

					List<TargetedResultBase> precursorResultList = precursorExecutor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
					UIMFTargetedMSMSWorkflowCollapseIMS precursorWorkflow = (UIMFTargetedMSMSWorkflowCollapseIMS)precursorExecutor.TargetedWorkflow;
					Dictionary<ChromPeak, XYData> precursorChromPeakToXYDataMap = precursorWorkflow.ChromPeakToXYDataMap;

					foreach (var precursorTargetedResultBase in precursorResultList)
					{
						if (precursorTargetedResultBase.ErrorDescription != null && precursorTargetedResultBase.ErrorDescription.Any())
						{
							continue;
						}
						if (precursorTargetedResultBase.FailedResult)
						{
							continue;
						}

						DatabaseSubTarget precursorTarget = precursorTargetedResultBase.Target as DatabaseSubTarget;

						foreach (var precursorPeakQualityData in precursorTargetedResultBase.ChromPeakQualityList)
						{
							if (precursorPeakQualityData.FitScore > 0.3) continue;

							ChromPeak precursorChromPeak = precursorPeakQualityData.Peak;
							XYData precursorProfileData = precursorChromPeakToXYDataMap[precursorChromPeak];
							precursorProfileData.NormalizeYData();
							XICProfile precursorXicProfile = new XICProfile(precursorPeakQualityData, precursorPeakQualityData.Peak);

							// TODO: Remove this hard-coded filter
							//if (precursorPeakQualityData.IsotopicProfile.ChargeState != 3 || precursorChromPeak.NETValue < 0.58 || precursorChromPeak.NETValue > 0.59) continue;

							DatabaseSubTargetResult precursorResult = new DatabaseSubTargetResult(precursorTarget, databaseTarget, precursorProfileData, precursorXicProfile);

							// Find out if this result relates to a previous result
							bool foundMatch = false;
							foreach (DatabaseMultipleSubTargetResult result in resultList)
							{
								if(result.DoesNewResultBelong(precursorResult))
								{
									result.AddNewResult(precursorResult);
									foundMatch = true;
									break;
								}
							}

							if (!foundMatch)
							{
								DatabaseMultipleSubTargetResult newResult = new DatabaseMultipleSubTargetResult(precursorResult);
								resultList.Add(newResult);
							}
						}
					}

					foreach (DatabaseMultipleSubTargetResult result in resultList)
					{
						XYData precursorProfileData = result.PrecursorResultRep.XYData;

						int chargeStateToSearch = result.ChargeStateList.Max();
						double elutionTime = result.ElutionTime;

						Console.WriteLine("Targeting NET = " + (float)elutionTime + "\tCS = " + chargeStateToSearch + "\tFit = " + result.PrecursorResultRep.XICProfile.DeconToolsFitScore);

						// Create fragment targets
						fragmentExecutor.Targets.TargetList = databaseTarget.CreateFragmentTargets((float)elutionTime, chargeStateToSearch);
						fragmentExecutor.Execute();

						List<TargetedResultBase> fragmentResultList = fragmentExecutor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
						UIMFTargetedMSMSWorkflowCollapseIMS fragmentWorkflow = (UIMFTargetedMSMSWorkflowCollapseIMS)fragmentExecutor.TargetedWorkflow;
						Dictionary<ChromPeak, XYData> fragmentChromPeakToXYDataMap = fragmentWorkflow.ChromPeakToXYDataMap;

						foreach (var fragmentTargetedResultBase in fragmentResultList)
						{
							if (fragmentTargetedResultBase.ErrorDescription != null && fragmentTargetedResultBase.ErrorDescription.Any())
							{
								continue;
							}
							if (fragmentTargetedResultBase.FailedResult)
							{
								continue;
							}

							DatabaseFragmentTarget fragmentTarget = fragmentTargetedResultBase.Target as DatabaseFragmentTarget;

							foreach (var fragmentPeakQualityData in fragmentTargetedResultBase.ChromPeakQualityList)
							{
								if (fragmentPeakQualityData.FitScore > 0.3) continue;

								ChromPeak fragmentChromPeak = fragmentPeakQualityData.Peak;
								XYData fragmentProfileData = fragmentChromPeakToXYDataMap[fragmentChromPeak];
								fragmentProfileData.NormalizeYData();
								XICProfile fragmentXicProfile = new XICProfile(fragmentPeakQualityData, fragmentPeakQualityData.Peak);

								DataUtil.CorrelateXYData(precursorProfileData, fragmentProfileData, 1, out slope, out intercept, out rSquared);
								//Console.WriteLine(fragmentTarget.Fragment.ToString());
								//Console.WriteLine("NET = " + fragmentChromPeak.NETValue + "\tScan = " + fragmentChromPeak.XValue + "\tPeakWidth = " + fragmentChromPeak.Width + "\tFit = " + fragmentPeakQualityData.FitScore + "\tMass = " + fragmentPeakQualityData.IsotopicProfile.MonoIsotopicMass);
								//Console.WriteLine("Slope = " + slope + "\tIntercept = " + intercept + "\tRSquared = " + rSquared);
								//fragmentProfileData.Display();

								if (rSquared >= 0.5)
								{
									Console.WriteLine(fragmentTarget.Fragment.ToString());
									DatabaseFragmentTargetResult fragmentResult = new DatabaseFragmentTargetResult(fragmentTarget, fragmentProfileData, fragmentXicProfile, fragmentPeakQualityData);
									result.FragmentResultList.Add(fragmentResult);
								}
							}
						}

						foreach (DatabaseFragmentTargetResult fragmentResult in result.FragmentResultList)
						{
							List<MSPeak> isotopicPeakList = fragmentResult.PeakQualityData.IsotopicProfile.Peaklist;
							foreach (MSPeak msPeak in isotopicPeakList)
							{
								Console.WriteLine(msPeak.XValue + "\t\t" + msPeak.Height);
							}
						}
					}
				}
			}
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
					executor.Targets.TargetList = databaseTarget.CreatePrecursorTargets();

					foreach (var target in executor.Targets.TargetList)
					{
						Console.WriteLine(target.ChargeState + "\t" + target.EmpiricalFormula + "\t" + target.Code + "\t" + target.MonoIsotopicMass + "\t" + target.MZ + "\t" + target.NormalizedElutionTime);
					}

					// Execute the workflow with the new targets
					executor.Execute();

					// TODO: Save results? Print results?
					List<TargetedResultBase> resultList = executor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
					UIMFTargetedMSMSWorkflowCollapseIMS workflow = (UIMFTargetedMSMSWorkflowCollapseIMS)executor.TargetedWorkflow;
					PrintResults(resultList, workflow.ChromPeakToXYDataMap);
				}
			}
		}

		[Test]
		public void TestTargetFragments()
		{
			short minChargeState = 1;
			short maxChargeState = 3;
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
					executor.Targets.TargetList = databaseTarget.CreateFragmentTargets(0.584f, maxChargeState);

					//foreach (var target in executor.Targets.TargetList)
					//{
					//    Console.WriteLine(target.ChargeState + "\t" + target.EmpiricalFormula + "\t" + target.Code + "\t" + target.MonoIsotopicMass + "\t" + target.MZ + "\t" + target.NormalizedElutionTime);
					//}

					// Execute the workflow with the new targets
					executor.Execute();

					// TODO: Save results? Print results?
					List<TargetedResultBase> resultList = executor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
					UIMFTargetedMSMSWorkflowCollapseIMS workflow = (UIMFTargetedMSMSWorkflowCollapseIMS)executor.TargetedWorkflow;
					PrintResults(resultList, workflow.ChromPeakToXYDataMap);
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

		private void PrintResults(IEnumerable<TargetedResultBase> resultList, Dictionary<ChromPeak, XYData> chromPeakToXYDataMap)
		{
			foreach (var targetedResultBase in resultList)
			{
				TargetBase target = targetedResultBase.Target;

				DatabaseFragmentTarget fragmentTarget = target as DatabaseFragmentTarget;
				if(fragmentTarget != null)
				{
					Console.WriteLine("************************************************");
					Console.WriteLine(fragmentTarget.Fragment.ToString());
				}
				else
				{
					Console.WriteLine("************************************************");
					Console.WriteLine("CS = " + target.ChargeState + "\tMass = " + target.MonoIsotopicMass + "\tm/z = " + target.MZ);
				}
				

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
							ChromPeak peak = peakQualityData.Peak;

							XYData profileData = chromPeakToXYDataMap[peak];
							XICProfile xicProfile = new XICProfile(peakQualityData, peakQualityData.Peak);
							Console.WriteLine("NET = " + peak.NETValue + "\tScan = " + peak.XValue + "\tPeakWidth = " + peak.Width + "\tFit = " + peakQualityData.FitScore + "\tMass = " + peakQualityData.IsotopicProfile.MonoIsotopicMass);
							profileData.Display();
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
