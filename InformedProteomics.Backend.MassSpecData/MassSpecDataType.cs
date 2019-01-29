namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Enum of distinctly handled file types
    /// </summary>
    public enum MassSpecDataType
    {
        /// <summary>
        /// Thermo .raw data, read with XcaliburReader if DLLs are available
        /// </summary>
        // ReSharper disable once IdentifierTypo
        XCaliburRun,

        /// <summary>
        /// mzML file, read with MzMLReader usually
        /// </summary>
        MzMLFile,

        /// <summary>
        /// PBF file, read with PbfLcMsRun
        /// </summary>
        PbfFile,

        /// <summary>
        /// Deconvoluted PBF file, read with DPbfLcMsRun
        /// </summary>
        DeconvolutedPbfFile,

        /// <summary>
        /// Other source file, read with ProteoWizard if DLLs are available
        /// </summary>
        Unknown,
    }
}
