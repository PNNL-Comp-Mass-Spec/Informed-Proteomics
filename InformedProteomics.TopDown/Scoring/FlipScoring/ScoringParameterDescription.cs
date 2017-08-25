using System.IO;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.TopDown.Scoring.FlipScoring
{
    /// <summary>
    /// Class containing the description for a scoring parameter set.
    /// This determines which scoring parameters to select.
    /// </summary>
    public class ScoringParameterDescription
    {
        /// <summary>
        /// Gets or sets the MS/MS activation/fragmentation method.
        /// </summary>
        public ActivationMethod ActivationMethod { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether.
        /// </summary>
        public bool IsTopDown { get; set; }

        /// <summary>
        /// Gets the path to the correct scoring parameter set based on the description.
        /// </summary>
        /// <param name="paramPath">
        /// The base path for the scoring parameters.
        /// The default (if left null) is bin/scoringParams.
        /// </param>
        /// <returns>The path as a string.</returns>
        public string GetPath(string paramPath = null)
        {
            paramPath = paramPath ?? "scoringParams";
            string topdown = this.IsTopDown ? "topDown" : "bottomUp";
            return Path.Combine(paramPath, string.Format("{0}_{1}", this.ActivationMethod, topdown));
        }

        /// <summary>
        /// Overloaded equals. Compares based on <see cref="ActivationMethod" /> and <see cref="IsTopDown" />.
        /// </summary>
        /// <param name="other">The scoringparameter description to compare to.</param>
        /// <returns>A value indicating whether the scoring parameter is equal to the other parameter.</returns>
        protected bool Equals(ScoringParameterDescription other)
        {
            return this.ActivationMethod == other.ActivationMethod && this.IsTopDown == other.IsTopDown;
        }

        /// <summary>
        /// Overloaded object-based equals.
        /// </summary>
        /// <param name="obj">The object to compare to.</param>
        /// <returns>A value indicating if this scoring parameter object is equal to the other object.</returns>
        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return this.Equals((ScoringParameterDescription)obj);
        }

        /// <summary>
        /// Overloaded hash code based on <see cref="ActivationMethod" /> and <see cref="IsTopDown" />.
        /// </summary>
        /// <returns>Hashed scoring parameters.</returns>
        public override int GetHashCode()
        {
            unchecked
            {
                return ((int)this.ActivationMethod * 397) ^ this.IsTopDown.GetHashCode();
            }
        }
    }
}
