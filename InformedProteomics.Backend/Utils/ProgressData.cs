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
        /// The current percent progress of the task. Updated using <see cref="Report(double,string)"/> or variants
        /// </summary>
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
            private set { _percent = value; }
        }

        /// <summary>
        /// If the progress reporting will be blocked into ranges
        /// Setting this to "true" will reset MinPercentage and MaxPercentage to 0.
        /// </summary>
        public bool IsPartialRange { get; set; }

        /// <summary>
        /// Must be less than current MaxPercentage
        /// </summary>
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
        private double _minPercentage = 0;
        private double _maxPercentage = 100;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="progress">The progress object that "ProgressData.Report" should call "Report" on</param>
        public ProgressData(IProgress<ProgressData> progress = null)
        {
            LastUpdated = DateTime.MinValue;
            UpdateFrequencySeconds = 0.0001;
            IsPartialRange = false;
            ProgressObj = progress;
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
                IsPartialRange = true;

                _minPercentage = 0;
                if (_maxPercentage >= 100.0)
                {
                    _maxPercentage = 100;
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

            // Trigger an update, with the proper minimum value for the range
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
