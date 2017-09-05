namespace InformedProteomics.Backend.Data.Biology
{
    /// <summary>
    /// Interface for biological components that have a elemental composition
    /// </summary>
    public interface IMolecule : IMatter
    {
        /// <summary>
        /// Elemental composition of the biological component
        /// </summary>
        Composition.Composition Composition { get; }
    }
}
