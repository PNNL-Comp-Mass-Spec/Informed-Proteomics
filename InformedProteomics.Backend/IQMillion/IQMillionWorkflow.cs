﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using DeconTools.Backend;
using DeconTools.Backend.Core;
using DeconTools.Backend.ProcessingTasks;
using DeconTools.Backend.ProcessingTasks.ChromatogramProcessing;
using DeconTools.Backend.ProcessingTasks.FitScoreCalculators;
using DeconTools.Backend.ProcessingTasks.PeakDetectors;
using DeconTools.Backend.ProcessingTasks.ResultValidators;
using DeconTools.Backend.ProcessingTasks.Smoothers;
using DeconTools.Backend.ProcessingTasks.TargetedFeatureFinders;
using DeconTools.Backend.ProcessingTasks.TheorFeatureGenerator;
using DeconTools.Workflows.Backend.Core;

namespace InformedProteomics.Backend.IQMillion
{
    internal class IQMillionWorkflow : TargetedWorkflow
    {
        #region Constructors

        public IQMillionWorkflow(Run run, TargetedWorkflowParameters parameters)
            : base(run, parameters)
        {

        }

        public IQMillionWorkflow(TargetedWorkflowParameters parameters)
            : base(parameters)
        {
        }

        #endregion


        #region IWorkflow Members

        protected override void DoMainInitialization()
        {
            ValidateParameters();

            _theorFeatureGen = new JoshTheorFeatureGenerator(DeconTools.Backend.Globals.LabellingType.NONE, 0.005);
            _chromGen = new PeakChromatogramGenerator(_workflowParameters.ChromGenTolerance, _workflowParameters.ChromGeneratorMode,
                                                      DeconTools.Backend.Globals.IsotopicProfileType.UNLABELLED,
                                                      _workflowParameters.ChromGenToleranceUnit)
            {
                TopNPeaksLowerCutOff = 0.333,
                NETWindowWidthForAlignedData = (float)_workflowParameters.ChromNETTolerance * 2,
                NETWindowWidthForNonAlignedData = (float)_workflowParameters.ChromNETTolerance * 2
            };

            //only

            const bool allowNegativeValues = false;
            _chromSmoother = new SavitzkyGolaySmoother(_workflowParameters.ChromSmootherNumPointsInSmooth, 2, allowNegativeValues);
            _chromPeakDetector = new ChromPeakDetector(_workflowParameters.ChromPeakDetectorPeakBR, _workflowParameters.ChromPeakDetectorSigNoise);

            _chromPeakSelector = CreateChromPeakSelector(_workflowParameters);

            _iterativeTFFParameters = new IterativeTFFParameters();
            _iterativeTFFParameters.ToleranceInPPM = _workflowParameters.MSToleranceInPPM;

            _msfeatureFinder = new IterativeTFF(_iterativeTFFParameters);
            _fitScoreCalc = new MassTagFitScoreCalculator();
            _resultValidator = new ResultValidatorTask();
            _chromatogramCorrelator = new ChromatogramCorrelator(_workflowParameters.ChromSmootherNumPointsInSmooth, 0.01, _workflowParameters.ChromGenTolerance);

            ChromatogramXYData = new XYData();
            MassSpectrumXYData = new XYData();
            ChromPeaksDetected = new List<ChromPeak>();
        }

        public override void DoWorkflow()
        {
            Result = Run.ResultCollection.GetTargetedResult(Run.CurrentMassTag);
            Result.ResetResult();

            ExecuteTask(_theorFeatureGen);
            ExecuteTask(_chromGen);
            ExecuteTask(_chromSmoother);
            updateChromDataXYValues(Run.XYData);

            ExecuteTask(_chromPeakDetector);
            UpdateChromDetectedPeaks(Run.PeakList);

            //UpdateChromPeaksWithXYData();

            ExecuteTask(_chromPeakSelector);
            ChromPeakSelected = Result.ChromPeakSelected;

            Result.ResetMassSpectrumRelatedInfo();

            ExecuteTask(MSGenerator);
            updateMassSpectrumXYValues(Run.XYData);

            ExecuteTask(_msfeatureFinder);

            ExecuteTask(_fitScoreCalc);
            ExecuteTask(_resultValidator);

            if (_workflowParameters.ChromatogramCorrelationIsPerformed)
            {
                ExecuteTask(_chromatogramCorrelator);
            }
        }

        protected override DeconTools.Backend.Globals.ResultType GetResultType()
        {
            return DeconTools.Backend.Globals.ResultType.BASIC_TARGETED_RESULT;
        }

        #endregion

    }
}

