using System;

namespace InformedProteomics.Backend.Utils
{
    public class ProgressData
    {
        public string Status { get; set; }

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
            set { _percent = value; }
        }

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

        public ProgressData()
        {
            LastUpdated = DateTime.MinValue;
            UpdateFrequencySeconds = 0.0001;
            IsPartialRange = false;
        }

        /// <summary>
        /// Change to a new range block
        /// </summary>
        /// <param name="newMaxPercentage">New max percent for range, must be greater than current max percent.</param>
        public void StepRange(double newMaxPercentage)
        {
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
    }
}
