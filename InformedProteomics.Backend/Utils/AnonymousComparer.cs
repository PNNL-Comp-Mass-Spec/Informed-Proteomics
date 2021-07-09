using System;
using System.Collections.Generic;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Class used for helping create an IComparer class for binary search
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class AnonymousComparer<T> : IComparer<T>
    {
        private readonly Comparison<T> _mComparison;

        /// <summary>
        /// Constructor that requires the Comparison be passed in
        /// </summary>
        /// <param name="comparison">The Comparison to be used for the binary search</param>
        public AnonymousComparer(Comparison<T> comparison)
        {
            _mComparison = comparison ?? throw new ArgumentNullException(nameof(comparison));
        }

        /// <summary>
        /// Compares 2 objects using the Comparison passed in when creating the AnonymousComparer class
        /// </summary>
        /// <param name="x">The first object</param>
        /// <param name="y">The second object</param>
        /// <returns>
        /// Less than zero if the first object precedes the second.
        /// Zero if the objects occur in the same position.
        /// Greater than zero if the first object follows the second.
        /// </returns>
        public int Compare(T x, T y)
        {
            return _mComparison(x, y);
        }
    }
}
