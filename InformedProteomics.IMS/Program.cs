using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using DeconTools.Backend;
using DeconTools.Backend.Core;
using DeconTools.Backend.ProcessingTasks.ChromatogramProcessing;
using DeconTools.Workflows.Backend.Core;
using InformedProteomics.Backend.Data.Results;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Targets;
using InformedProteomics.Backend.IMS;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.IMS
{
    class Program
    {
        public static void Main(string[] args)
        {
			//TextWriter textWriter = new StreamWriter("profiles.csv");

			//string line = "ion";
			//for(int i = 509; i <= 527; i+=2)
			//{
			//    line += "," + i;
			//}
			//textWriter.WriteLine(line);

			short minChargeState = 1;
			short maxChargeState = 5;
			double minMz = 200;
			double maxMz = 2500;

			//string datasetPath = @"\\protoapps\UserData\Slysz\Standard_Testing\Targeted_FeatureFinding\UIMF_Targeted_MSMS_Testing\RawData\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo.UIMF";
        	string datasetPath = args[0];

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

			TextWriter dtaWriter = new StreamWriter("test_dta.txt");

			AminoAcidSet aminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);

			// TODO: Read in a peptide list?
			//List<string> peptideSequenceList = new List<string> { "EYANQFMWEYSTNYGQAPLSLLVSYTK" };
			IEnumerable<string> peptideSequenceList = ReadPeptideList(args[1]);

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

					//foreach (var target in precursorExecutor.Targets.TargetList)
					//{
					//    Console.WriteLine(target.ChargeState + "\t" + target.EmpiricalFormula + "\t" + target.Code + "\t" +
					//                      target.MonoIsotopicMass + "\t" + target.MZ + "\t" + target.NormalizedElutionTime);
					//}

					// Execute the workflow with the new targets
					precursorExecutor.Execute();

					List<TargetedResultBase> precursorResultList = precursorExecutor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
					UIMFTargetedMSMSWorkflowCollapseIMS precursorWorkflow = (UIMFTargetedMSMSWorkflowCollapseIMS)precursorExecutor.TargetedWorkflow;
					Dictionary<ChromPeak, XYData> precursorChromPeakToXYDataMap = precursorWorkflow.ChromPeakToXYDataMap;
					int maxChargeStateOfPrecursor = 0;

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
							int chargeState = precursorPeakQualityData.IsotopicProfile.ChargeState;

							// I don't believe any peaks that are not wide enough to be valid peptides
							if (precursorChromPeak.Width < 5 || precursorChromPeak.Width > 20) continue;

							DatabaseSubTargetResult precursorResult = new DatabaseSubTargetResult(precursorTarget, databaseTarget, precursorProfileData, precursorXicProfile, precursorPeakQualityData.FitScore);

							if (chargeState > maxChargeStateOfPrecursor) maxChargeStateOfPrecursor = chargeState;

							// Find out if this result relates to a previous result
							bool foundMatch = false;
							foreach (DatabaseMultipleSubTargetResult result in resultList)
							{
								if (result.DoesNewResultBelong(precursorResult))
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

						Console.WriteLine("Targeting Sequence = " + result.PrecursorResultRep.DatabaseSubTarget.Code + " NET = " + (float)elutionTime + "\tCS = " + chargeStateToSearch + "\tFit = " + result.PrecursorResultRep.XICProfile.DeconToolsFitScore);
						//precursorProfileData.Display();

						// Create fragment targets
						fragmentExecutor.Targets.TargetList = databaseTarget.CreateFragmentTargets((float)elutionTime, chargeStateToSearch);
						fragmentExecutor.Execute();

						List<TargetedResultBase> fragmentResultList = fragmentExecutor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
						UIMFTargetedMSMSWorkflowCollapseIMS fragmentWorkflow = (UIMFTargetedMSMSWorkflowCollapseIMS)fragmentExecutor.TargetedWorkflow;
						Dictionary<ChromPeak, XYData> fragmentChromPeakToXYDataMap = fragmentWorkflow.ChromPeakToXYDataMap;

						foreach (TargetedResultBase fragmentTargetedResultBase in fragmentResultList)
						{
							DatabaseFragmentTarget fragmentTarget = fragmentTargetedResultBase.Target as DatabaseFragmentTarget;
							//Console.WriteLine(fragmentTarget.Fragment.ToString());

							if (fragmentTargetedResultBase.ErrorDescription != null && fragmentTargetedResultBase.ErrorDescription.Any())
							{
								//Console.WriteLine(fragmentTargetedResultBase.ErrorDescription);
								continue;
							}
							if (fragmentTargetedResultBase.FailedResult)
							{
								//Console.WriteLine("Failed Result.");
								continue;
							}

							foreach (ChromPeakQualityData fragmentPeakQualityData in fragmentTargetedResultBase.ChromPeakQualityList)
							{
								//Console.WriteLine("Fit Score = " + fragmentPeakQualityData.FitScore);
								if (fragmentPeakQualityData.FitScore > 0.3) continue;

								ChromPeak fragmentChromPeak = fragmentPeakQualityData.Peak;
								XYData fragmentProfileData = fragmentChromPeakToXYDataMap[fragmentChromPeak];
								fragmentProfileData.NormalizeYData();
								XICProfile fragmentXicProfile = new XICProfile(fragmentPeakQualityData, fragmentPeakQualityData.Peak);

								//Console.WriteLine(fragmentTarget.Fragment.ToString());
								//Console.WriteLine("NET = " + fragmentChromPeak.NETValue + "\tScan = " + fragmentChromPeak.XValue + "\tPeakWidth = " + fragmentChromPeak.Width + "\tFit = " + fragmentPeakQualityData.FitScore + "\tMass = " + fragmentPeakQualityData.IsotopicProfile.MonoIsotopicMass);
								DataUtil.CorrelateXYData(precursorProfileData, fragmentProfileData, 1, out slope, out intercept, out rSquared);
								//Console.WriteLine("Slope = " + slope + "\tIntercept = " + intercept + "\tRSquared = " + rSquared);
								//fragmentProfileData.Display();

								//String label = fragmentTarget.Fragment.IonSymbol;
								//for (int i = 0; i < fragmentTarget.Fragment.ChargeState; i++)
								//{
								//    label += "+";
								//}

								//WriteFragmentProfileToFile(fragmentProfileData, label, textWriter, 509, 527);

								if (rSquared >= 0.4)
								{
									//Console.WriteLine(fragmentTarget.Fragment.ToString());
									DatabaseFragmentTargetResult fragmentResult = new DatabaseFragmentTargetResult(fragmentTarget, fragmentProfileData, fragmentXicProfile, fragmentPeakQualityData);
									result.FragmentResultList.Add(fragmentResult);
								}
							}
						}

						DtaUtil.AppendToDtaFile(result, dtaWriter);
					}
				}
			}

			dtaWriter.Close();
			//textWriter.Close();
        }

		private static BasicTargetedWorkflowExecutorParameters CreateWorkflowExecutorParameters()
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

		private static UIMFTargetedMSMSWorkflowCollapseIMSParameters CreateWorkflowParameters()
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

		private static IEnumerable<string> ReadPeptideList(string fileLocation)
		{
			List<string> peptideList = new List<string>();

			TextReader textReader = new StreamReader(fileLocation);

			string line = "";
			while ((line = textReader.ReadLine()) != null)
			{
				peptideList.Add(line);
			}

			textReader.Close();

			return peptideList;
		}
    }
}
