using System;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Scoring.Interfaces
{
    /// <summary>
    /// This is an interface for scorers that creates user-consumable scores.
    /// </summary>
    public interface IInformedScorer
    {
        /// <summary>
        /// Gets the score cut off for filtering out PSMs before they are passed into the generating function.
        /// </summary>
        double ScoreCutOff { get; }

        /// <summary>
        /// Gets the number of sequence fragments found in the spectrum.
        /// </summary>
        /// <param name="sequence">The sequence to count fragments for.</param>
        /// <returns>The number of matched fragments.</returns>
        int GetNumMatchedFragments(Sequence sequence);

        /// <summary>
        /// Gets the user-presentable score that will be displayed on results file.
        /// </summary>
        /// <param name="sequence">The sequence to calculate a score for.</param>
        /// <returns>The score.</returns>
        double GetUserVisibleScore(Sequence sequence);
    }
}
