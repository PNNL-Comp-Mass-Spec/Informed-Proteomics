using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public class PbfReader : IMassSpecDataReader
    {
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
                while (reader.BaseStream.Position != reader.BaseStream.Length)
                {
                    var spec = ReadSpectrumFromPbf(reader);
                    //Console.WriteLine("{0}\t{1}\t{2}", spec.ScanNum, spec.MsLevel, spec.Peaks.Length);
                    yield return spec;
                }
            }
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

            var productSpec = spec as ProductSpectrum;
            if (productSpec != null)    // product spectrum
            {
                // Activation method: 1
                writer.Write((byte)productSpec.ActivationMethod);
                // Isolation window target m/z: 8
                writer.Write(productSpec.IsolationInfo.IsolationWindowTargetMz);
                // Isolation window lower offset: 8
                writer.Write(productSpec.IsolationInfo.IsolationWindowLowerOffset);
                // Isolation window uppoer offset: 8
                writer.Write(productSpec.IsolationInfo.IsolationWindowUpperOffset);
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
            if (msLevel > 1)
            {
                var activationMethod = (ActivationMethod)reader.ReadByte();
                var isolationWindowTargetMz = reader.ReadDouble();
                var isolationWindowLowerOffset = reader.ReadDouble();
                var isolationWindowUpperOffset = reader.ReadDouble();
                var peakList = ReadPeakList(reader);
                return new ProductSpectrum(peakList, scanNum)
                    {
                        MsLevel = msLevel,
                        ActivationMethod = activationMethod,
                        IsolationInfo = new IsolationInfo(
                            isolationWindowTargetMz,
                            isolationWindowLowerOffset,
                            isolationWindowUpperOffset)
                    };
            }
            else
            {
                var peakList = ReadPeakList(reader);
                return new Spectrum(peakList, scanNum);
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
