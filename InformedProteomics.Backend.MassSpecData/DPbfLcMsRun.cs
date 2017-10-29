using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Spectrometry;
using PRISM;

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

        /// <summary>
        /// File extension used for this type
        /// </summary>
        protected override string FileExtensionVirtual => FileExtensionConst;

        /// <summary>
        /// File extension - overridden. See <see cref="FileExtensionConst"/> for static access.
        /// </summary>
        public override bool ContainsChromatograms => false;

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
        public override string GetCheckPbfFilePath(string specFilePath, out string pbfPath, out string fileName, out string tempPath)
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
        /// <param name="keepDataReaderOpen"></param>
        public DPbfLcMsRun(string specFileName, IMassSpecDataReader msdr, string pbfFileName = null,
            double precursorSignalToNoiseRatioThreshold = 0.0, double productSignalToNoiseRatioThreshold = 0.0,
            IProgress<ProgressData> progress = null, bool keepDataReaderOpen = false)
            : base(precursorSignalToNoiseRatioThreshold, productSignalToNoiseRatioThreshold)
        {
            GetPbfFile(specFileName, msdr, pbfFileName, progress, keepDataReaderOpen);
        }

        /// <summary>
        /// Reads a Spectrum from the DPBF file with isotope peaks populated
        /// </summary>
        /// <param name="fullData">The file used to create the DPBF file (throws <see cref="System.ArgumentException"/> if it is not) - PBF file preferred.</param>
        /// <param name="scanNum">The scan to read</param>
        /// <returns></returns>
        /// <exception cref="System.ArgumentException">If the checksum of the source file does not match the checksum stored in the DPBF file</exception>
        public DeconvolutedSpectrum GetSpectrumWithIsotopePeaks(IMassSpecDataReader fullData, int scanNum)
        {
            if (SrcFileChecksum != fullData.SrcFileChecksum || (fullData is PbfLcMsRun && SrcFileChecksum != ((PbfLcMsRun)fullData).PbfFileChecksum))
            {
                throw new ArgumentException("Supplied file was not used to create this DPBF file!", nameof(fullData));
            }

            if (!(GetSpectrum(scanNum, true) is DeconvolutedSpectrum spec))
            {
                return null;
            }

            var pbfSpec = fullData.ReadMassSpectrum(scanNum, true);

            foreach (var peakBase in spec.Peaks)
            {
                var peak = peakBase as DeconvolutedPeak;
                if (peak == null)
                {
                    continue;
                }

                peak.SetObservedPeaksFromSpectrum(pbfSpec);
            }

            return spec;
        }

        /// <summary>
        /// Read a spectrum from the current position in <paramref name="reader"/>, with the option to only read the metadata.
        /// </summary>
        /// <param name="reader"></param>
        /// <param name="includePeaks"></param>
        /// <returns></returns>
        protected internal override Spectrum ReadSpectrum(BinaryReader reader, bool includePeaks = true)
        {
            // Must reflect all changes to WriteSpectrum

            var scanNum = reader.ReadInt32();
            var c = new char[NativeIdLength];
            reader.Read(c, 0, NativeIdLength);
            var nativeId = (new string(c)).Trim();
            var msLevel = reader.ReadByte();
            var elutionTime = reader.ReadDouble();
            var tic = reader.ReadSingle();

            double? precursorMass = reader.ReadDouble();
            if (Math.Abs(precursorMass.Value) < float.Epsilon) precursorMass = null;
            int? precursorCharge = reader.ReadInt32();
            if (precursorCharge == 0) precursorCharge = null;
            var activationMethod = (ActivationMethod)reader.ReadByte();
            var isolationWindowTargetMz = reader.ReadDouble();
            var isolationWindowLowerOffset = reader.ReadDouble();
            var isolationWindowUpperOffset = reader.ReadDouble();

            var peakList = new List<DeconvolutedPeak>();
            var numPeaks = reader.ReadInt32();

            for (var i = 0; i < numPeaks; i++)
            {
                var mz = reader.ReadDouble();
                var intensity = reader.ReadSingle();
                var charge = reader.ReadInt32();
                var corr = reader.ReadSingle();
                var dist = reader.ReadSingle();
                var peak = new DeconvolutedPeak(mz, intensity, charge, corr, dist);
                var isoPeakCount = (int)reader.ReadByte();
                peak.ObservedPeakIndices.Capacity = isoPeakCount;
                for (var j = 0; j < isoPeakCount; j++)
                {
                    peak.ObservedPeakIndices.Add(reader.ReadUInt16());
                }
                peakList.Add(peak);
            }

            var spec = new DeconvolutedSpectrum(peakList, scanNum)
            {
                ActivationMethod = activationMethod,
                IsolationWindow = new IsolationWindow(
                    isolationWindowTargetMz,
                    isolationWindowLowerOffset,
                    isolationWindowUpperOffset,
                    precursorMass,
                    precursorCharge
                ),
                MsLevel = msLevel,
                ElutionTime = elutionTime,
                NativeId = nativeId,
                TotalIonCurrent = tic
            };
            return spec;
        }

        /// <summary>
        /// Write the supplied spectrum to the current position in <paramref name="writer"/>
        /// </summary>
        /// <param name="specIn"></param>
        /// <param name="writer"></param>
        protected internal override void WriteSpectrum(Spectrum specIn, BinaryWriter writer)
        {
            // All changes made here must be duplicated to ReadSpectrum() and GetPeakMetadataForSpectrum()
            if (!(specIn is DeconvolutedSpectrum spec))
            {
                throw new ArgumentException("Input spectrum must be DeconvolutedSpectrum!", nameof(specIn));
            }

            // scan number: 4
            writer.Write(spec.ScanNum);

            // NativeID: 50
            // pad or truncate to keep in limit (may have to change in future...)
            writer.Write(spec.NativeId.PadRight(NativeIdLength).ToCharArray(0, NativeIdLength), 0, NativeIdLength);

            // ms level: 1
            writer.Write(Convert.ToByte(spec.MsLevel));

            // elution time: 4
            writer.Write(spec.ElutionTime);

            // Total Ion Current: 4
            writer.Write(Convert.ToSingle(spec.TotalIonCurrent));

            var isolationWindow = spec.IsolationWindow;
            // precursor mass: 8
            writer.Write(isolationWindow.MonoisotopicMz ?? 0.0);
            // precursor charge: 4
            writer.Write(isolationWindow.Charge ?? 0);
            // Activation method: 1
            writer.Write((byte)spec.ActivationMethod);
            // Isolation window target m/z: 8
            writer.Write(isolationWindow.IsolationWindowTargetMz);
            // Isolation window lower offset: 8
            writer.Write(isolationWindow.IsolationWindowLowerOffset);
            // Isolation window upper offset: 8
            writer.Write(isolationWindow.IsolationWindowUpperOffset);

            // Guarantee sorted peaks.
            //var peaks = spec.Peaks.ToList();
            Array.Sort(spec.Peaks);
            //peaks.Sort();
            // Number of peaks: 4
            writer.Write(spec.Peaks.Length);
            //writer.Write(peaks.Count);
            foreach (var peakIn in spec.Peaks)
            //foreach (var peak in peaks)
            {
                var peak = peakIn as DeconvolutedPeak;
                if (peak == null)
                {
                    throw new ArgumentException("Input spectrum peaks array must contain only DeconvolutedPeaks!");
                }
                // m/z: 8
                writer.Write(peak.Mz);
                // intensity: 4
                writer.Write(Convert.ToSingle(peak.Intensity));
                // charge: 4
                writer.Write(peak.Charge);
                // correlation: 4
                writer.Write(Convert.ToSingle(peak.Corr));
                // distance: 4
                writer.Write(Convert.ToSingle(peak.Dist));

                // Output indices of isotope peaks where index (in PbfLcMsRun) is < 65535
                var isotopePeaksInRange = peak.ObservedPeakIndices.Where(x => x >= ushort.MinValue && x <= ushort.MaxValue).ToList();
                if (isotopePeaksInRange.Count != peak.ObservedPeakIndices.Count)
                {
                    WarnOnce(WriteWarnings.PeakIndexCannotBeStored, "Cannot output all observed peaks for some peaks - more than 65535 peaks in observed spectrum");
                }
                // Count of isotope peaks, as byte: 1
                writer.Write((byte)isotopePeaksInRange.Count);
                // Output indices, as unsigned shorts
                foreach (var index in isotopePeaksInRange)
                {
                    writer.Write((ushort) index);
                }
            }
        }

        private enum WriteWarnings
        {
            PeakIndexCannotBeStored
        }
        private readonly List<WriteWarnings> warningTracker = new List<WriteWarnings>();

        private void WarnOnce(WriteWarnings code, string message)
        {
            if (warningTracker.Contains(code))
            {
                return;
            }
            warningTracker.Add(code);
            System.Console.WriteLine("WARNING: {0}!");
        }
    }
}
