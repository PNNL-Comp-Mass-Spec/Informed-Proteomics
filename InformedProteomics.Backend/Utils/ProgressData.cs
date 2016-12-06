using System;

namespace InformedProteomics.Backend.Utils
{
    public class ProgressData
    {
        public string Status { get; set; }
        public string StatusInternal { get; set; }
        public IProgress<ProgressData> ProgressObj { get; set; }

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

        public double UpdateFrequencySeconds { get; set; }

        // static for the case of multiple ProgressData objects being fed to "Progress.Report()"
        public static DateTime LastUpdated { get; private set; }
        private static string _updateLock = String.Empty; // Only because we need a reference type for a lock

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
                    _maxPercentage = 0;
                }
            }
            if (newStatus != null)
            {
                Status = newStatus;
            }
            CheckSetMinMaxRange(_maxPercentage, newMaxPercentage);
        }

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

        public ProgressData UpdatePercent(double pct)
        {
            Percent = pct;
            return this;
        }

        public bool ShouldUpdate()
        {
            var update = false;
            lock (_updateLock)
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
