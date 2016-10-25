using System;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.MassSpecData
{
    /// <summary>
    /// Reduced version of PbfLcMsRun for holding deconvoluted spectra
    /// </summary>
    public class DPbfLcMsRun : PbfLcMsRun
    {
        /// <summary>
        /// File extension
        /// </summary>
        public new const string FileExtensionConst = ".dpbf";
        protected override string FileExtensionVirtual { get { return FileExtensionConst; } }

        /// <summary>
        /// File extension - overridden. See <see cref="FileExtensionConst"/> for static access.
        /// </summary>
        public override bool ContainsChromatograms { get { return false; } }

        /// <summary>
        /// Function to convert a spectra file name/path to a *.pbf name, even when it has multiple extensions (i.e., .mzML.gz)
        /// </summary>
        /// <param name="specFileName"></param>
        /// <returns></returns>
        /// <remarks>It is recommended that "MassSpecDataReaderFactory.NormalizeDatasetPath" be called prior to calling this function, and that the returned string be used instead of the original path</remarks>
        public new static string GetPbfFileName(string specFileName)
        {
            return MassSpecDataReaderFactory.ChangeExtension(specFileName, FileExtensionConst);
        }

        /// <summary>
        /// Gets valid possible pbf file paths
        /// </summary>
        /// <param name="specFilePath">Path to the spectra file</param>
        /// <param name="pbfPath">Path to the default pbf file (in the same folder as the spectra file dataset)</param>
        /// <param name="fileName"></param>
        /// <param name="tempPath"></param>
        /// <returns>The default path to the pbf file, unless a valid pbf file exists at the temp path</returns>
        public new static string GetCheckPbfFilePath(string specFilePath, out string pbfPath, out string fileName, out string tempPath)
        {
            return GetCheckPbfFilePath(specFilePath, out pbfPath, out fileName, out tempPath, FileExtensionConst);
        }

        /// <summary>
        /// Constructor for opening a DPBF file
        /// </summary>
        /// <param name="specFileName"></param>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        public DPbfLcMsRun(string specFileName, double precursorSignalToNoiseRatioThreshold = 0.0,
            double productSignalToNoiseRatioThreshold = 0.0)
            : base(precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold)
        {
            OpenPbfFile(specFileName);
        }

        /// <summary>
        /// Constructor for creating and/or opening a DPBF file
        /// </summary>
        /// <param name="specFileName"></param>
        /// <param name="msdr"></param>
        /// <param name="pbfFileName"></param>
        /// <param name="precursorSignalToNoiseRatioThreshold"></param>
        /// <param name="productSignalToNoiseRatioThreshold"></param>
        /// <param name="progress"></param>
        public DPbfLcMsRun(string specFileName, IMassSpecDataReader msdr, string pbfFileName = null,
            double precursorSignalToNoiseRatioThreshold = 0.0, double productSignalToNoiseRatioThreshold = 0.0,
            IProgress<ProgressData> progress = null)
            : base(precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold)
        {
            GetPbfFile(specFileName, msdr, pbfFileName, progress);
        }
    }
}
