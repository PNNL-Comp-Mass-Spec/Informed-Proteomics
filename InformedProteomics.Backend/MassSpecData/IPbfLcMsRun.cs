namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Interface for PbfLcMsRun, primarily for cases where PbfLcMsRun properties/functions are needed within the backend
    /// </summary>
    public interface IPbfLcMsRun : ILcMsRun
    {
        /// <summary>
        /// SHA-1 Checksum of the pbf file, calculated on first access to this property - lowercase, hex only
        /// </summary>
        string PbfFileChecksum { get; }
    }
}
