using System;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Add functionality to progress reporting using <see cref="IProgress{T}"/>
    /// </summary>
    public class ProgressData
    {
        /// <summary>
        /// Status string - for reporting textual information about the current task
        /// </summary>
        public string Status { get; set; }

        /// <summary>
        /// Internal status string - for tracking nested progress status
        /// </summary>
        public string StatusInternal { get; set; }

        /// <summary>
        /// Referenced <see cref="IProgress{T}"/> object, that all updates are pushed out to.
        /// </summary>
        public IProgress<ProgressData> ProgressObj { get; set; }

        /// <summary>
        /// When true, anything that will cause progress to go backwards will cause an exception; otherwise, such changes are silently handled to prevent backwards progress. Should not be true in general release code.
        /// </summary>
        public bool ThrowExceptionOnBackwardsProgress
        {
            get { return _throwExceptionOnBackwardsProgress; }
            set
            {
                if (_throwExceptionOnBackwardsProgress != value)
                {
                    _throwExceptionOnBackwardsProgress = value;
                    CheckForwardOnlyLogic();
                }
            }
        }

        /// <summary>
        /// When set to true, logic is used that will prevent progress from jumping backwards (errors are silently ignored; see <see cref="ThrowExceptionOnBackwardsProgress"/> to trigger exceptions instead)
        /// </summary>
        public bool PreventBackwardsProgress
        {
            get { return _preventBackwardsProgress; }
            set
            {
                if (_preventBackwardsProgress != value)
                {
                    _preventBackwardsProgress = value;
                    CheckForwardOnlyLogic();
                }
            }
        }

        /// <summary>
        /// The current percent progress of the task. Updated using <see cref="Report(double,string)"/> or variants
        /// </summary>
        /// <remarks>Value between 0 and 100</remarks>
        public double Percent
        {
            get
            {
                if (IsPartialRange)
                {
                    return _percent * ((MaxPercentage - MinPercentage) / 100.0) + MinPercentage;
                }
                return _percent;
            }
            private set
            {
                if (_useForwardOnlyLogic)
                {
                    if (ThrowExceptionOnBackwardsProgress && _percent > value && value < 3 * _lastPercentChange)
                    {
                        throw new Exception("Progress percent change will cause anomaly in progress reporting!");
                    }

                    if (value > _percent)
                    {
                        _lastPercentChange = value - _percent;
                        _percent = value;
                    }
                }
                else
                {
                    _percent = value;
                }
            }
        }

        /// <summary>
        /// If the progress reporting will be blocked into ranges
        /// Setting this to "true" will reset MinPercentage and MaxPercentage to 0.
        /// </summary>
        public bool IsPartialRange {
            get { return _isPartialRange; }
            set
            {
                if (_isPartialRange == value)
                {
                    return;
                }
                // if setting to false, change percent accordingly to prevent sudden, odd jumps
                if (!value)
                {
                    _percent = Percent;
                    _isPartialRange = value;
                }
                else
                {
                    _isPartialRange = value;
                    // if setting to true, change the min percentage to prevent sudden, odd jumps
                    if (ThrowExceptionOnBackwardsProgress && _percent > _maxPercentage)
                    {
                        throw new Exception("Progress step change from non-partial to partial will cause anomaly in progress reporting!");
                    }

                    MinPercentage = _percent;
                }
            } }

        /// <summary>
        /// Must be less than current MaxPercentage
        /// </summary>
        /// <remarks>Will set IsPartialRange to true</remarks>
        public double MinPercentage
        {
            get
            {
                return _minPercentage;
            }
            set
            {
                CheckSetMinMaxRange(value, _maxPercentage);
            }
        }

        /// <summary>
        /// Must be greater than current MinPercentage
        /// </summary>
        /// <remarks>Will set IsPartialRange to true</remarks>
        public double MaxPercentage
        {
            get
            {
                return _maxPercentage;
            }
            set
            {
                CheckSetMinMaxRange(_minPercentage, value);
            }
        }

        /// <summary>
        /// Throttling for console output - used with ShouldUpdate() to provide a simple throttle to reduce the console output
        /// </summary>
        public double UpdateFrequencySeconds { get; set; }

        /// <summary>
        /// Last output time, for throttling updates for console output
        /// </summary>
        /// <remarks>static for the case of multiple ProgressData objects being fed to "Progress.Report()"</remarks>
        public static DateTime LastUpdated { get; private set; }
        private static readonly string UpdateLock = string.Empty; // Only because we need a reference type for a lock

        private double _percent = 0;
        private double _lastPercentChange = 0;
        private bool _throwExceptionOnBackwardsProgress = false;
        private bool _preventBackwardsProgress = false;
        private bool _useForwardOnlyLogic = false;
        private bool _isPartialRange = false;

        private void CheckForwardOnlyLogic()
        {
            if (ThrowExceptionOnBackwardsProgress)
            {
                PreventBackwardsProgress = true;
            }
            _useForwardOnlyLogic = ThrowExceptionOnBackwardsProgress || PreventBackwardsProgress;
        }

        /// <summary>
        /// Track if a partial range (not 0-100%) has been set previously. This should never be set to false outside of object construction.
        /// </summary>
        private bool HasUsedPartialRange
        {
            get { return _hasUsedPartialRangeWithAReallyLongAndNastyNameSoThatNoOneEverWantsToUseUtBesidesWhereItIsSupposedToBeUsed; }
            set
            {
                if (!_hasUsedPartialRangeWithAReallyLongAndNastyNameSoThatNoOneEverWantsToUseUtBesidesWhereItIsSupposedToBeUsed && value)
                {
                    _hasUsedPartialRangeWithAReallyLongAndNastyNameSoThatNoOneEverWantsToUseUtBesidesWhereItIsSupposedToBeUsed = true;
                }
            }
        }

        /// <summary>
        /// Backing variable for HasUsedPartialRange. ONLY USE INSIDE OF HasUsedPartialRange GETTER/SETTER.
        /// </summary>
        private bool _hasUsedPartialRangeWithAReallyLongAndNastyNameSoThatNoOneEverWantsToUseUtBesidesWhereItIsSupposedToBeUsed = false;

        private double _minPercentage = 0;
        private double _maxPercentage = 100;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="progress">The progress object that "ProgressData.Report" should call "Report" on</param>
        /// <param name="preventBackwardsProgress">Set to false to disable the logic preventing reverse progress</param>
        public ProgressData(IProgress<ProgressData> progress = null, bool preventBackwardsProgress = true)
        {
            LastUpdated = DateTime.MinValue;
            UpdateFrequencySeconds = 0.0001;
            IsPartialRange = false;
            ProgressObj = progress;
            PreventBackwardsProgress = preventBackwardsProgress;
            if (progress == null)
            {
                ProgressObj = new Progress<ProgressData>();
            }
        }

        /// <summary>
        /// Change to a new range block
        /// </summary>
        /// <param name="newMaxPercentage">New max percent for range, must be greater than current max percent.</param>
        /// <param name="newStatus">Updated status string, null for no update</param>
        /// <remarks>Will set IsPartialRange to true</remarks>
        /// <remarks>If current max percent is 100, the new max percent can be any value between 0 and 100</remarks>
        public void StepRange(double newMaxPercentage, string newStatus = null)
        {
            if (!IsPartialRange)
            {
                if (ThrowExceptionOnBackwardsProgress && _percent > _maxPercentage)
                {
                    throw new Exception("Progress step change from non-partial to partial will cause anomaly in progress reporting!");
                }

                _minPercentage = _percent;
                // Set the current max range to zero, so that the call to CheckSetMinMaxRange works as we need it to.
                // Only do it when the max percentage is 100+, to allow potential removal then re-add of percent stepping.
                if (_maxPercentage >= 100.0)
                {
                    _maxPercentage = 0;
                }
            }
            if (newStatus != null)
            {
                Status = newStatus;
            }
            CheckSetMinMaxRange(_maxPercentage, newMaxPercentage);
        }

        /// <summary>
        /// Perform validity checks on new min/max percent values, and then set them accordingly
        /// </summary>
        /// <param name="newMin"></param>
        /// <param name="newMax"></param>
        private void CheckSetMinMaxRange(double newMin, double newMax)
        {
            if (_useForwardOnlyLogic)
            {
                if (ThrowExceptionOnBackwardsProgress)
                {
                    if (_percent > _maxPercentage && !IsPartialRange)
                    {
                        throw new Exception("Progress step change from non-partial to partial will cause anomaly in progress reporting!");
                    }
                    if (newMin < _minPercentage)
                    {
                        throw new Exception("Progress min change will cause anomaly in progress reporting!");
                    }
                    if (newMax < _maxPercentage && HasUsedPartialRange)
                    {
                        throw new Exception("Progress max change will cause anomaly in progress reporting!");
                    }
                }

                // Prevent progress range from decreasing
                if (newMin < _minPercentage)
                {
                    newMin = _minPercentage;
                }
                // newMax less than old max is a problem, unless this is the first setting of the max for partial ranges
                if (newMax < _maxPercentage && HasUsedPartialRange)
                {
                    newMax = _maxPercentage;
                }
            }
            if (newMin > 0.0 || newMax < 100.0)
            {
                HasUsedPartialRange = true;
            }
            if (newMax > newMin)
            {
                _minPercentage = newMin;
                _maxPercentage = newMax;
            }
            if (_minPercentage < 0)
            {
                _minPercentage = 0;
            }
            if (_maxPercentage > 100.0)
            {
                _maxPercentage = 100;
            }

            _isPartialRange = true;

            // Trigger an update, with the proper minimum value for the range
            _lastPercentChange = 0.0;
            _percent = 0; // Proper minimum value for range, changed outside of the control checks
            Report(0.0);
        }

        /// <summary>
        /// Update percent, and return object. For single-lining a progress update and report with <see cref="IProgress{T}.Report"/>
        /// </summary>
        /// <param name="pct"></param>
        /// <returns></returns>
        [Obsolete("Use Report() instead, with ProgressObj set.")]
        public ProgressData UpdatePercent(double pct)
        {
            Percent = pct;
            return this;
        }

        /// <summary>
        /// Check function to limit output frequency, when outputting to console.
        /// </summary>
        /// <returns></returns>
        public bool ShouldUpdate()
        {
            var update = false;
            lock (UpdateLock)
            {
                if (DateTime.UtcNow >= LastUpdated.AddSeconds(UpdateFrequencySeconds))
                {
                    LastUpdated = DateTime.UtcNow;
                    update = true;
                }
            }
            return update;
        }

        /// <summary>
        /// Updates the status, then calls the stored progress object's "Report"
        /// </summary>
        /// <param name="newStatus">Updated status string</param>
        public void Report(string newStatus)
        {
            Status = newStatus;
            ProgressObj.Report(this);
        }

        /// <summary>
        /// Updates the percent, then calls the stored progress object's "Report"
        /// </summary>
        /// <param name="pct">percent progress, 0 to 100</param>
        /// <param name="newStatus">Updated status string, null for no update</param>
        public void Report(double pct, string newStatus = null)
        {
            if (pct > 100)
            {
                pct = 100;
            }
            // drop large spikes that may have been triggered by multithreading
            if (_useForwardOnlyLogic && pct - _percent > 70 && !(pct.Equals(0.0) && _percent.Equals(100.0)))
            {
                pct = _percent;
            }

            Percent = pct;
            if (newStatus != null)
            {
                Status = newStatus;
            }
            ProgressObj.Report(this);
        }

        /// <summary>
        /// Updates the percent, then calls the stored progress object's "Report"
        /// </summary>
        /// <param name="pct">percent progress, 0 to 1</param>
        /// <param name="newStatus">Updated status string, null for no update</param>
        public void ReportDecimal(double pct, string newStatus = null)
        {
            Report(pct * 100.0, newStatus);
        }
        /// <summary>
        /// Updates the percent, then calls the stored progress object's "Report"
        /// </summary>
        /// <param name="count">The count progress, or numerator</param>
        /// <param name="total">The total number of objects to be counted, or denominator</param>
        /// <param name="newStatus">Updated status string, null for no update</param>
        public void Report(double count, double total, string newStatus = null)
        {
            ReportDecimal(count / total, newStatus);
        }
    }
}
