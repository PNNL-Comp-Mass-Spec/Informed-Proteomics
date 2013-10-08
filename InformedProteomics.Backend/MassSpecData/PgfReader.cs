using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.MassSpecData
{
    public class PgfReader: IMassSpecDataReader
    {
        private readonly string _specFileName;

        public PgfReader(string specFileName)
        {
            _specFileName = specFileName;
        }

        public IEnumerable<Spectrum> ReadAllSpectra()
        {
            var scanNum = -1;
            var msLevel = -1;
            var activation = ActivationMethod.Unknown;
            var isolationWindowTargetMz = 0.0;
            var isolationWindowLowerOffset = 0.0;
            var isolationWindowUpperOffset = 0.0;
            List<Peak> peakList = null;

            using (var fileStream = File.Open(_specFileName, FileMode.Open, FileAccess.Read, FileShare.Read))
            using(var bufferedStream = new BufferedStream(fileStream))
            using (var reader = new StreamReader(bufferedStream))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    if (line.Length == 0) continue;
                    if (line.Equals("BEGIN IONS"))
                    {
                        peakList = new List<Peak>();
                    }
                    else if (line.StartsWith("SCANS="))
                    {
                        scanNum = Convert.ToInt32(line.Substring("SCANS=".Length));
                    }
                    else if (line.StartsWith("_MSLEVEL="))
                    {
                        msLevel = Convert.ToInt32(line.Substring("_MSLEVEL=".Length));
                    }
                    else if (line.StartsWith("_ACTIVATION="))
                    {
                        var activationStr = line.Substring("_ACTIVATION=".Length);
                        if (activationStr.Equals("CID")) activation = ActivationMethod.CID;
                        else if (activationStr.Equals("ETD")) activation = ActivationMethod.ETD;
                        else if (activationStr.Equals("HCD")) activation = ActivationMethod.HCD;
                        else if (activationStr.Equals("ECD")) activation = ActivationMethod.ECD;
                        else if (activationStr.Equals("PQD")) activation = ActivationMethod.PQD;
                        else activation = ActivationMethod.Unknown;
                    }
                    else if (line.StartsWith("_ISOTARGET="))
                    {
                        isolationWindowTargetMz = Convert.ToDouble(line.Substring("_ISOTARGET=".Length));
                    }
                    else if (line.StartsWith("_ISOLOWER="))
                    {
                        isolationWindowLowerOffset = Convert.ToDouble(line.Substring("_ISOLOWER=".Length));
                    }
                    else if (line.StartsWith("_ISOUPPER="))
                    {
                        isolationWindowUpperOffset = Convert.ToDouble(line.Substring("_ISOUPPER=".Length));
                    }
                    else if (line.Equals("END IONS"))
                    {
                        if (msLevel == 1)
                        {
                            yield return new Spectrum(peakList, scanNum);
                        }
                        else
                        {
                            yield return new ProductSpectrum(peakList, scanNum)
                            {
                                MsLevel = msLevel,
                                ActivationMethod = activation,
                                IsolationInfo = new IsolationInfo(
                                    isolationWindowTargetMz,
                                    isolationWindowLowerOffset,
                                    isolationWindowUpperOffset)
                            };
                        }
                    }
                    else
                    {
                        var token = line.Split('\t');
                        if (token.Length != 2) continue;
                        var mz = Convert.ToSingle(token[0]);
                        var intensity = Convert.ToDouble(token[1]);
                        var peak = new Peak(mz, intensity);
                        if (peakList == null) yield break; // Wrong format! Should throw an error.
                        peakList.Add(peak);
                    }
                }
            }
        }

        public void Close()
        {
            //if (_reader != null) _reader.Close();
        }

        public static void WriteSpectrumAsPgf(Spectrum spec, StreamWriter writer)
        {
            writer.WriteLine("BEGIN IONS");
            writer.WriteLine("SCANS={0}", spec.ScanNum);
            writer.WriteLine("_MSLEVEL={0}", spec.MsLevel);
            var productSpec = spec as ProductSpectrum;
            if (productSpec != null)
            {
                writer.WriteLine("_ACTIVATION={0}", productSpec.ActivationMethod);
                writer.WriteLine("_ISOTARGET={0}", productSpec.IsolationInfo.IsolationWindowTargetMz);
                writer.WriteLine("_ISOLOWER={0}", productSpec.IsolationInfo.IsolationWindowLowerOffset);
                writer.WriteLine("_ISOUPPER={0}", productSpec.IsolationInfo.IsolationWindowUpperOffset);
            }
            foreach (var peak in spec.Peaks)
            {
                writer.WriteLine("{0}\t{1}", Convert.ToSingle(peak.Mz), Convert.ToSingle(peak.Intensity));
            }
            writer.WriteLine("END IONS");
        }
    }
}
