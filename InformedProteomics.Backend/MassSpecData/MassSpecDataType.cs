namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Enum of distinctly handled file types
    /// </summary>
    public enum MassSpecDataType
    {
        /// <summary>
        /// Thermo Finnigan .RAW data, read with XCaliburReader if DLLs are available
        /// </summary>
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
