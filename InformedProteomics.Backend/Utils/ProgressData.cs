using System;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Add functionality to progress reporting using <see cref="IProgress{T}"/>
    /// </summary>
    [Obsolete("Use PRISM.ProgressData", true)]
    public class ProgressData : PRISM.ProgressData
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="progress">The progress object that "ProgressData.Report" should call "Report" on</param>
        /// <param name="preventBackwardsProgress">Set to false to disable the logic preventing reverse progress</param>
        public ProgressData(IProgress<ProgressData> progress = null, bool preventBackwardsProgress = true) : base(progress as IProgress<PRISM.ProgressData>, preventBackwardsProgress)
        {
        }
    }
}
