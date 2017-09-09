namespace InformedProteomics.Backend.Database
{
    /// <summary>
    /// Type of peptide cleavage to use (for Top-Down searches)
    /// </summary>
    public enum InternalCleavageType
    {
        /// <summary>
        /// Only intact proteins
        /// </summary>
        /// <remarks>No cleavages will be performed on the protein sequences
        /// NOTE: Allows C-term cleavages up to the specified maximum</remarks>
        NoInternalCleavage = 0,

        /// <summary>
        /// Only sequences that match either the N terminus cleavages parameter or the C terminus cleavages parameter
        /// </summary>
        /// <remarks>If the C terminus has no cleavage, the N terminus can have cleavages up to the specified N terminus cleavages max;
        /// if the N terminus has no cleavages, the C terminus can have cleavages up to the specified C terminus cleavages max.
        /// NOTE: Any sequence minLength &lt;= length &lt;= maxLength, with C or N terminus cleavages</remarks>
        SingleInternalCleavage = 1,

        /// <summary>
        /// Any sequence in a protein that matches the parameters for N and C terminus cleavages
        /// </summary>
        /// <remarks>All cleavage combinations that meet the parameters for N and C terminus cleavage maximums will be checked.
        /// NOTE: Any sequence minLength &lt;= length &lt;= maxLength</remarks>
        MultipleInternalCleavages = 2
    }
}