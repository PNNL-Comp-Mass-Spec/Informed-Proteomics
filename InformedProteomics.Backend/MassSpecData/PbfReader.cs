using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public class PbfReader : IMassSpecDataReader
    {
        public const int FileFormatId = 150601;

        public PbfReader(string specFileName)
        {
            _specFileName = specFileName;
        }

        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            using (var fileStream = File.Open(_specFileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            using(var bufferedStream = new BufferedStream(fileStream))
            using (var reader = new BinaryReader(bufferedStream))
            {
                while (reader.BaseStream.Position != (reader.BaseStream.Length - sizeof(int)))
                {
                    var spec = ReadSpectrumFromPbf(reader);
                    //Console.WriteLine("{0}\t{1}\t{2}", spec.ScanNum, spec.MsLevel, spec.Peaks.Length);
                    yield return spec;
                }
            }
        }

        public static bool CheckFileFormatVersion(string filePath)
        {
            var fs = File.OpenRead(filePath);
            using (var reader = new BinaryReader(fs))
            {
                fs.Seek(-1 * sizeof(int), SeekOrigin.End);

                var fileFormatId = reader.ReadInt32();
                if (fileFormatId != FileFormatId) return false;
            }
            return true;
        }


        public void Close()
        {
        }

        public static void WriteSpectrumAsPbf(Spectrum spec, BinaryWriter writer)
        {
            // scan number: 4
            writer.Write(spec.ScanNum);

            // ms level: 1
            writer.Write(Convert.ToByte(spec.MsLevel));

            // elution time: 4
            writer.Write(spec.ElutionTime);

            var productSpec = spec as ProductSpectrum;
            if (productSpec != null)    // product spectrum
            {
                var isolationWindow = productSpec.IsolationWindow;
                // precursor mass: 8
                writer.Write(isolationWindow.MonoisotopicMz ?? 0.0);
                // precursor charge: 4
                writer.Write(isolationWindow.Charge ?? 0);
                // Activation method: 1
                writer.Write((byte)productSpec.ActivationMethod);
                // Isolation window target m/z: 8
                writer.Write(productSpec.IsolationWindow.IsolationWindowTargetMz);
                // Isolation window lower offset: 8
                writer.Write(productSpec.IsolationWindow.IsolationWindowLowerOffset);
                // Isolation window uppoer offset: 8
                writer.Write(productSpec.IsolationWindow.IsolationWindowUpperOffset);
            }
            // Number of peaks: 4
            writer.Write(spec.Peaks.Length);
            foreach (var peak in spec.Peaks)
            {
                // m/z: 8
                writer.Write(peak.Mz);
                // intensity: 4
                writer.Write(Convert.ToSingle(peak.Intensity));
            }
        }

        private static Spectrum ReadSpectrumFromPbf(BinaryReader reader)
        {
            var scanNum = reader.ReadInt32();
            var msLevel = reader.ReadByte();
            var elutionTime = reader.ReadDouble();

            if (msLevel > 1)
            {
                double? precursorMass = reader.ReadDouble();
                if (precursorMass == 0.0) precursorMass = null;
                int? precursorCharge = reader.ReadInt32();
                if (precursorCharge == 0) precursorCharge = null;
                var activationMethod = (ActivationMethod)reader.ReadByte();
                var isolationWindowTargetMz = reader.ReadDouble();
                var isolationWindowLowerOffset = reader.ReadDouble();
                var isolationWindowUpperOffset = reader.ReadDouble();
                var peakList = ReadPeakList(reader);
                return new ProductSpectrum(peakList, scanNum)
                    {
                        MsLevel = msLevel,
                        ElutionTime = elutionTime,
                        ActivationMethod = activationMethod,
                        IsolationWindow = new IsolationWindow(
                            isolationWindowTargetMz,
                            isolationWindowLowerOffset,
                            isolationWindowUpperOffset,
                            precursorMass,
                            precursorCharge
                            )
                    };
            }
            else
            {
                var peakList = ReadPeakList(reader);
                return new Spectrum(peakList, scanNum)
                    {
                        ElutionTime = elutionTime
                    };
            }
        }

        private static List<Peak> ReadPeakList(BinaryReader reader)
        {
            var peakList = new List<Peak>();
            var numPeaks = reader.ReadInt32();
            for (var i = 0; i < numPeaks; i++)
            {
                var mz = reader.ReadDouble();
                var intensity = reader.ReadSingle();
                peakList.Add(new Peak(mz, intensity));
            }
            return peakList;
        }

        private readonly string _specFileName;
    }
}
