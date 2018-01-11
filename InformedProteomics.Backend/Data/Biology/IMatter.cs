namespace InformedProteomics.Backend.Data.Biology
{
    /// <summary>
    /// Interface for biological components classified as Matter (having a mass)
    /// </summary>
    public interface IMatter
    {
        //double GetMass();

        /// <summary>
        /// Mass of the biological component
        /// </summary>
        double Mass { get; }
    }
}
