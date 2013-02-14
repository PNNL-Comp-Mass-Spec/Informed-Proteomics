using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using DeconTools.Backend;
using DeconTools.Backend.Core;
using DeconTools.Backend.ProcessingTasks.ChromatogramProcessing;
using DeconTools.Workflows.Backend.Core;
using InformedProteomics.Backend.Data.Results;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Targets;
using InformedProteomics.Backend.IMS;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestUtils
    {
        [Test]
        public void TestCompositions()
        {
            const string sequence = "PEPTIDE";
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var compositions = aaSet.GetCompositions(sequence);

            Composition[] compositionArr = compositions.ToArray();
            Assert.AreEqual(compositionArr.Count(), 1);
            Console.WriteLine(compositionArr[0]);
            Console.WriteLine(compositionArr[0].GetMass());
            Console.WriteLine(compositionArr[0].GetNominalMass());
            Assert.AreEqual(compositionArr[0].ToString(), "C34H51N7O14S0");
        }

        [Test]
        public void TestIQMillion()
        {
            //const short minChargeState = 1;
            //const short maxChargeState = 4;
            //const int minPrecursorNominalMass = 400;
            //const int maxPrecursorNominalMass = 2500;

            //const string datasetPath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";

            //var executorParameters = CreateWorkflowExecutorParameters();
            //var precursorWorkflowParameters = CreateWorkflowParameters();

            //// Create workflow for precursors
            //var precursorExecutor = new UIMFTargetedWorkflowExecutor(executorParameters, precursorWorkflowParameters, datasetPath)
            //                            {
            //                                Targets = new TargetCollection()
            //                            };

            //for (int precursorNominalMass = minPrecursorNominalMass;
            //     precursorNominalMass <= maxPrecursorNominalMass;
            //     precursorNominalMass++)
            //{
            //    // creates nominal mass targets
            //    precursorExecutor.Targets.TargetList = 
                
            //}

            //// Create workflow for fragments
            //var fragmentWorkflowParameters = CreateWorkflowParameters();
            //fragmentWorkflowParameters.ChromNETTolerance = 0.005;
            //var fragmentExecutor = new UIMFTargetedWorkflowExecutor(executorParameters, fragmentWorkflowParameters, datasetPath)
            //                           {
            //                               Targets = new TargetCollection()
            //                           };

            //TextWriter dtaWriter = new StreamWriter("test_NM_dta.txt");

            //var aminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);

            ////IEnumerable<string> peptideSequenceList = ReadPeptideList(@"../../../TestFiles/BSA_ST.txt");
            ////var peptideSequenceList = new List<string> { "CCAADDKEACFAVEGPK" };

            //{
            //    var resultList = new List<DatabaseMultipleSubTargetResult>();
            //    Composition peptideComposition = composition + Composition.H2O;

            //    // Create the database target
            //    precursorExecutor.Targets.TargetList = databaseTarget.CreatePrecursorTargets();

            //    foreach (var target in precursorExecutor.Targets.TargetList)
            //    {
            //        Console.WriteLine(target.ChargeState + "\t" + target.EmpiricalFormula + "\t" + target.Code +
            //                            "\t" +
            //                            target.MonoIsotopicMass + "\t" + target.MZ + "\t" +
            //                            target.NormalizedElutionTime);
            //    }

            //    // Execute the workflow with the new targets
            //    precursorExecutor.Execute();

            //    List<TargetedResultBase> precursorResultList =
            //        precursorExecutor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
            //    var precursorWorkflow =
            //        (UIMFTargetedMSMSWorkflowCollapseIMS)precursorExecutor.TargetedWorkflow;
            //    Dictionary<ChromPeak, XYData> precursorChromPeakToXYDataMap = precursorWorkflow.ChromPeakToXYDataMap;
            //    int maxChargeStateOfPrecursor = 0;

            //    foreach (TargetedResultBase precursorTargetedResultBase in precursorResultList)
            //    {
            //        //Console.WriteLine("Checking Target...");
            //        if (precursorTargetedResultBase.ErrorDescription != null &&
            //            precursorTargetedResultBase.ErrorDescription.Any())
            //        {
            //            //Console.WriteLine("Precursor error: " + precursorTargetedResultBase.ErrorDescription);
            //            continue;
            //        }
            //        if (precursorTargetedResultBase.FailedResult)
            //        {
            //            //Console.WriteLine("Precursor failed result: " + precursorTargetedResultBase.FailureType);
            //            continue;
            //        }

            //        var precursorTarget = precursorTargetedResultBase.Target as DatabaseSubTarget;

            //        foreach (var precursorPeakQualityData in precursorTargetedResultBase.ChromPeakQualityList)
            //        {
            //            if (precursorPeakQualityData.FitScore > 0.4) continue;

            //            ChromPeak precursorChromPeak = precursorPeakQualityData.Peak;
            //            XYData precursorProfileData = precursorChromPeakToXYDataMap[precursorChromPeak];
            //            precursorProfileData.NormalizeYData();
            //            var precursorXicProfile = new XICProfile(precursorPeakQualityData, precursorPeakQualityData.Peak);
            //            int chargeState = precursorPeakQualityData.IsotopicProfile.ChargeState;
            //            //Console.WriteLine("NET = " + precursorChromPeak.NETValue + "\tScan = " + precursorChromPeak.XValue + "\tPeakWidth = " + precursorChromPeak.Width + "\tFit = " + precursorPeakQualityData.FitScore + "\tMass = " + precursorPeakQualityData.IsotopicProfile.MonoIsotopicMass);

            //            // I don't believe any peaks that are not wide enough to be valid peptides
            //            if (precursorChromPeak.Width < 5 || precursorChromPeak.Width > 20) continue;

            //            // TODO: Remove this hard-coded filter
            //            //if (precursorPeakQualityData.IsotopicProfile.ChargeState != 3 || precursorChromPeak.NETValue < 0.58 || precursorChromPeak.NETValue > 0.59) continue;

            //            var precursorResult = new DatabaseSubTargetResult(precursorTarget,
            //                                                                databaseTarget,
            //                                                                precursorProfileData,
            //                                                                precursorXicProfile,
            //                                                                precursorPeakQualityData
            //                                                                    .FitScore);

            //            if (chargeState > maxChargeStateOfPrecursor) maxChargeStateOfPrecursor = chargeState;

            //            // Find out if this result relates to a previous result
            //            bool foundMatch = false;
            //            foreach (DatabaseMultipleSubTargetResult result in resultList)
            //            {
            //                if (result.DoesNewResultBelong(precursorResult))
            //                {
            //                    result.AddNewResult(precursorResult);
            //                    foundMatch = true;
            //                    break;
            //                }
            //            }

            //            if (!foundMatch)
            //            {
            //                var newResult =
            //                    new DatabaseMultipleSubTargetResult(precursorResult);
            //                resultList.Add(newResult);
            //            }
            //        }
            //    }

            //    foreach (DatabaseMultipleSubTargetResult result in resultList)
            //    {
            //        XYData precursorProfileData = result.PrecursorResultRep.XYData;

            //        int chargeStateToSearch = result.ChargeStateList.Max();
            //        double elutionTime = result.ElutionTime;

            //        //Console.WriteLine("Targeting Sequence = " + result.PrecursorResultRep.DatabaseSubTarget.Code + " NET = " + (float)elutionTime + "\tCS = " + chargeStateToSearch + "\tFit = " + result.PrecursorResultRep.XICProfile.DeconToolsFitScore);
            //        //precursorProfileData.Display();

            //        // Create fragment targets
            //        fragmentExecutor.Targets.TargetList = databaseTarget.CreateFragmentTargets((float)elutionTime,
            //                                                                                    chargeStateToSearch);
            //        fragmentExecutor.Execute();

            //        List<TargetedResultBase> fragmentResultList =
            //            fragmentExecutor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
            //        var fragmentWorkflow =
            //            (UIMFTargetedMSMSWorkflowCollapseIMS)fragmentExecutor.TargetedWorkflow;
            //        Dictionary<ChromPeak, XYData> fragmentChromPeakToXYDataMap =
            //            fragmentWorkflow.ChromPeakToXYDataMap;

            //        foreach (TargetedResultBase fragmentTargetedResultBase in fragmentResultList)
            //        {
            //            var fragmentTarget =
            //                fragmentTargetedResultBase.Target as DatabaseFragmentTarget;
            //            //Console.WriteLine(fragmentTarget.Fragment.ToString());

            //            if (fragmentTargetedResultBase.ErrorDescription != null &&
            //                fragmentTargetedResultBase.ErrorDescription.Any())
            //            {
            //                //Console.WriteLine(fragmentTargetedResultBase.ErrorDescription);
            //                continue;
            //            }
            //            if (fragmentTargetedResultBase.FailedResult)
            //            {
            //                //Console.WriteLine("Failed Result.");
            //                continue;
            //            }

            //            foreach (
            //                ChromPeakQualityData fragmentPeakQualityData in
            //                    fragmentTargetedResultBase.ChromPeakQualityList)
            //            {
            //                //Console.WriteLine("Fit Score = " + fragmentPeakQualityData.FitScore);
            //                if (fragmentPeakQualityData.FitScore > 0.5) continue;
            //                //if (fragmentPeakQualityData.FitScore >= 1) continue;

            //                ChromPeak fragmentChromPeak = fragmentPeakQualityData.Peak;
            //                XYData fragmentProfileData = fragmentChromPeakToXYDataMap[fragmentChromPeak];
            //                fragmentProfileData.NormalizeYData();
            //                var fragmentXicProfile = new XICProfile(fragmentPeakQualityData,
            //                                                                fragmentPeakQualityData.Peak);

            //                //Console.WriteLine(fragmentTarget.Fragment.ToString());
            //                //Console.WriteLine("NET = " + fragmentChromPeak.NETValue + "\tScan = " + fragmentChromPeak.XValue + "\tPeakWidth = " + fragmentChromPeak.Width + "\tFit = " + fragmentPeakQualityData.FitScore + "\tMass = " + fragmentPeakQualityData.IsotopicProfile.MonoIsotopicMass);
            //                double slope;
            //                double intercept;
            //                double rSquared;
            //                DataUtil.CorrelateXYData(precursorProfileData, fragmentProfileData, 1, out slope,
            //                                            out intercept, out rSquared);
            //                //Console.WriteLine("Slope = " + slope + "\tIntercept = " + intercept + "\tRSquared = " + rSquared);
            //                //fragmentProfileData.Display();

            //                //String label = fragmentTarget.Fragment.IonSymbol;
            //                //for (int i = 0; i < fragmentTarget.Fragment.ChargeState; i++)
            //                //{
            //                //    label += "+";
            //                //}

            //                //WriteFragmentProfileToFile(fragmentProfileData, label, textWriter, 509, 527);

            //                if (rSquared >= 0.4)
            //                {
            //                    if (fragmentTarget != null)
            //                    {
            //                        textWriter.WriteLine(result.PrecursorResultRep.DatabaseSubTarget.Code + "," +
            //                                                result.PrecursorResultRep.DatabaseSubTarget.EmpiricalFormula +
            //                                                "," + elutionTime + "," + chargeStateToSearch + "," +
            //                                                result.IsotopicFitScore + "," +
            //                                                fragmentTarget.Fragment.IonSymbol + "," +
            //                                                fragmentTarget.Fragment.ChargeState + "," + rSquared + "," +
            //                                                fragmentPeakQualityData.FitScore);
            //                        //Console.WriteLine(fragmentTarget.Fragment.ToString());
            //                        var fragmentResult =
            //                            new DatabaseFragmentTargetResult(fragmentTarget, fragmentProfileData,
            //                                                                fragmentXicProfile, fragmentPeakQualityData);
            //                        result.FragmentResultList.Add(fragmentResult);
            //                    }
            //                }
            //            }
            //        }

            //        DtaUtil.AppendToDtaFile(result, dtaWriter);
            //    }
            //}

            //dtaWriter.Close();
            //textWriter.Close();
        }

        [Test]
        public void TestNominalMassXICGeneration()
        {
            TextWriter textWriter = new StreamWriter("fragmentDataNM.csv");
            textWriter.WriteLine("Sequence,EmpFormula,NET,Charge,Fit,Fragment,FragmentCharge,RSquared,FragmentFit");

            const short minChargeState = 1;
            const short maxChargeState = 5;
            const double minMz = 200;
            const double maxMz = 2500;

            //string datasetPath = @"\\protoapps\UserData\Slysz\Standard_Testing\Targeted_FeatureFinding\UIMF_Targeted_MSMS_Testing\RawData\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo.UIMF";
            //string datasetPath = @"\\protoapps\UserData\Slysz\Standard_Testing\Targeted_FeatureFinding\UIMF_Targeted_MSMS_Testing\RawData\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo_Collision_Energy_Collapsed.UIMF";
            const string datasetPath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";

            var executorParameters = CreateWorkflowExecutorParameters();
            var workflowParameters = CreateWorkflowParameters();

            // Create workflow for precursors
            var precursorExecutor = new UIMFTargetedWorkflowExecutor(executorParameters, workflowParameters, datasetPath)
                                        {
                                            Targets = new TargetCollection ()
                                        };

            // Create workflow for fragments
            var fragmentWorkflowParameters = CreateWorkflowParameters();
            fragmentWorkflowParameters.ChromNETTolerance = 0.005;
            var fragmentExecutor = new UIMFTargetedWorkflowExecutor(executorParameters, fragmentWorkflowParameters, datasetPath);
            fragmentExecutor.Targets = new TargetCollection();

            double slope;
            double intercept;
            double rSquared;

            TextWriter dtaWriter = new StreamWriter("test_NM_dta.txt");

            var aminoAcidSet = new AminoAcidSet(Modification.Carbamidomethylation);

            //IEnumerable<string> peptideSequenceList = ReadPeptideList(@"../../../TestFiles/BSA_ST.txt");
            List<string> peptideSequenceList = new List<string> { "CCAADDKEACFAVEGPK" };

            foreach (string peptideSequence in peptideSequenceList)
            {
                IEnumerable<Composition> compositions = aminoAcidSet.GetCompositions(peptideSequence);

                foreach (Composition composition in compositions)
                {
                    var resultList = new List<DatabaseMultipleSubTargetResult>();
                    Composition peptideComposition = composition + Composition.H2O;

                    // Create the database target
                    var sequence = new Sequence(peptideComposition, peptideSequence, aminoAcidSet);

                    var databaseTarget = new DatabaseTarget(sequence, minMz, maxMz, minChargeState, maxChargeState);

                    // Get the list of precursor targets from the DatabaseTarget object and pass them in
                    precursorExecutor.Targets.TargetList = databaseTarget.CreatePrecursorTargets();

                    foreach (var target in precursorExecutor.Targets.TargetList)
                    {
                        Console.WriteLine(target.ChargeState + "\t" + target.EmpiricalFormula + "\t" + target.Code +
                                          "\t" +
                                          target.MonoIsotopicMass + "\t" + target.MZ + "\t" +
                                          target.NormalizedElutionTime);
                    }

                    // Execute the workflow with the new targets
                    precursorExecutor.Execute();

                    List<TargetedResultBase> precursorResultList =
                        precursorExecutor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
                    var precursorWorkflow =
                        (UIMFTargetedMSMSWorkflowCollapseIMS) precursorExecutor.TargetedWorkflow;
                    Dictionary<ChromPeak, XYData> precursorChromPeakToXYDataMap = precursorWorkflow.ChromPeakToXYDataMap;
                    int maxChargeStateOfPrecursor = 0;

                    foreach (TargetedResultBase precursorTargetedResultBase in precursorResultList)
                    {
                        //Console.WriteLine("Checking Target...");
                        if (precursorTargetedResultBase.ErrorDescription != null &&
                            precursorTargetedResultBase.ErrorDescription.Any())
                        {
                            //Console.WriteLine("Precursor error: " + precursorTargetedResultBase.ErrorDescription);
                            continue;
                        }
                        if (precursorTargetedResultBase.FailedResult)
                        {
                            //Console.WriteLine("Precursor failed result: " + precursorTargetedResultBase.FailureType);
                            continue;
                        }

                        var precursorTarget = precursorTargetedResultBase.Target as DatabaseSubTarget;

                        foreach (var precursorPeakQualityData in precursorTargetedResultBase.ChromPeakQualityList)
                        {
                            if (precursorPeakQualityData.FitScore > 0.4) continue;

                            ChromPeak precursorChromPeak = precursorPeakQualityData.Peak;
                            XYData precursorProfileData = precursorChromPeakToXYDataMap[precursorChromPeak];
                            precursorProfileData.NormalizeYData();
                            var precursorXicProfile = new XICProfile(precursorPeakQualityData, precursorPeakQualityData.Peak);
                            int chargeState = precursorPeakQualityData.IsotopicProfile.ChargeState;
                            //Console.WriteLine("NET = " + precursorChromPeak.NETValue + "\tScan = " + precursorChromPeak.XValue + "\tPeakWidth = " + precursorChromPeak.Width + "\tFit = " + precursorPeakQualityData.FitScore + "\tMass = " + precursorPeakQualityData.IsotopicProfile.MonoIsotopicMass);

                            // I don't believe any peaks that are not wide enough to be valid peptides
                            if (precursorChromPeak.Width < 5 || precursorChromPeak.Width > 20) continue;

                            // TODO: Remove this hard-coded filter
                            //if (precursorPeakQualityData.IsotopicProfile.ChargeState != 3 || precursorChromPeak.NETValue < 0.58 || precursorChromPeak.NETValue > 0.59) continue;

                            var precursorResult = new DatabaseSubTargetResult(precursorTarget,
                                                                              databaseTarget,
                                                                              precursorProfileData,
                                                                              precursorXicProfile,
                                                                              precursorPeakQualityData
                                                                                  .FitScore);

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
                                var newResult =
                                    new DatabaseMultipleSubTargetResult(precursorResult);
                                resultList.Add(newResult);
                            }
                        }
                    }

                    foreach (DatabaseMultipleSubTargetResult result in resultList)
                    {
                        XYData precursorProfileData = result.PrecursorResultRep.XYData;

                        int chargeStateToSearch = result.ChargeStateList.Max();
                        double elutionTime = result.ElutionTime;

                        //Console.WriteLine("Targeting Sequence = " + result.PrecursorResultRep.DatabaseSubTarget.Code + " NET = " + (float)elutionTime + "\tCS = " + chargeStateToSearch + "\tFit = " + result.PrecursorResultRep.XICProfile.DeconToolsFitScore);
                        //precursorProfileData.Display();

                        // Create fragment targets
                        fragmentExecutor.Targets.TargetList = databaseTarget.CreateFragmentTargets((float) elutionTime,
                                                                                                   chargeStateToSearch);
                        fragmentExecutor.Execute();

                        List<TargetedResultBase> fragmentResultList =
                            fragmentExecutor.TargetedWorkflow.Run.ResultCollection.GetMassTagResults();
                        var fragmentWorkflow =
                            (UIMFTargetedMSMSWorkflowCollapseIMS) fragmentExecutor.TargetedWorkflow;
                        Dictionary<ChromPeak, XYData> fragmentChromPeakToXYDataMap =
                            fragmentWorkflow.ChromPeakToXYDataMap;

                        foreach (TargetedResultBase fragmentTargetedResultBase in fragmentResultList)
                        {
                            var fragmentTarget =
                                fragmentTargetedResultBase.Target as DatabaseFragmentTarget;
                            //Console.WriteLine(fragmentTarget.Fragment.ToString());

                            if (fragmentTargetedResultBase.ErrorDescription != null &&
                                fragmentTargetedResultBase.ErrorDescription.Any())
                            {
                                //Console.WriteLine(fragmentTargetedResultBase.ErrorDescription);
                                continue;
                            }
                            if (fragmentTargetedResultBase.FailedResult)
                            {
                                //Console.WriteLine("Failed Result.");
                                continue;
                            }

                            foreach (
                                ChromPeakQualityData fragmentPeakQualityData in
                                    fragmentTargetedResultBase.ChromPeakQualityList)
                            {
                                //Console.WriteLine("Fit Score = " + fragmentPeakQualityData.FitScore);
                                if (fragmentPeakQualityData.FitScore > 0.5) continue;
                                //if (fragmentPeakQualityData.FitScore >= 1) continue;

                                ChromPeak fragmentChromPeak = fragmentPeakQualityData.Peak;
                                XYData fragmentProfileData = fragmentChromPeakToXYDataMap[fragmentChromPeak];
                                fragmentProfileData.NormalizeYData();
                                var fragmentXicProfile = new XICProfile(fragmentPeakQualityData,
                                                                               fragmentPeakQualityData.Peak);

                                //Console.WriteLine(fragmentTarget.Fragment.ToString());
                                //Console.WriteLine("NET = " + fragmentChromPeak.NETValue + "\tScan = " + fragmentChromPeak.XValue + "\tPeakWidth = " + fragmentChromPeak.Width + "\tFit = " + fragmentPeakQualityData.FitScore + "\tMass = " + fragmentPeakQualityData.IsotopicProfile.MonoIsotopicMass);
                                DataUtil.CorrelateXYData(precursorProfileData, fragmentProfileData, 1, out slope,
                                                         out intercept, out rSquared);
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
                                    if (fragmentTarget != null)
                                    {
                                        textWriter.WriteLine(result.PrecursorResultRep.DatabaseSubTarget.Code + "," +
                                                             result.PrecursorResultRep.DatabaseSubTarget.EmpiricalFormula +
                                                             "," + elutionTime + "," + chargeStateToSearch + "," +
                                                             result.IsotopicFitScore + "," +
                                                             fragmentTarget.Fragment.IonSymbol + "," +
                                                             fragmentTarget.Fragment.ChargeState + "," + rSquared + "," +
                                                             fragmentPeakQualityData.FitScore);
                                        //Console.WriteLine(fragmentTarget.Fragment.ToString());
                                        var fragmentResult =
                                            new DatabaseFragmentTargetResult(fragmentTarget, fragmentProfileData,
                                                                             fragmentXicProfile, fragmentPeakQualityData);
                                        result.FragmentResultList.Add(fragmentResult);
                                    }
                                }
                            }
                        }

                        DtaUtil.AppendToDtaFile(result, dtaWriter);
                    }
                }
            }

            dtaWriter.Close();
            textWriter.Close();
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
            //workflowParameters.ChromGenTolerance = 1.0005/2;
            workflowParameters.ChromGenTolerance = 0.5;
            workflowParameters.ChromGenToleranceUnit = Globals.ToleranceUnit.MZ;
            
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

        private IEnumerable<string> ReadPeptideList(string fileLocation)
        {
            var peptideList = new List<string>();

            TextReader textReader = new StreamReader(fileLocation);

            string line = "";
            while ((line = textReader.ReadLine()) != null)
            {
                // If peptide in form of X.PEPTIDE.X, then grab the PEPTIDE string
                if (line.Contains("."))
                {
                    line = line.Split('.')[1];
                }

                // Cannot yet handle this mod
                if (line.Contains("@")) continue;

                // Remove the "!" since C is a static mod that is already expected
                line = line.Replace("!", "");

                peptideList.Add(line);
            }

            textReader.Close();

            return peptideList;
        }
    }
}
