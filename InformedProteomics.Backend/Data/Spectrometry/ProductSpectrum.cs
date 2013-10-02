using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class ProductSpectrum: Spectrum
    {
        public ProductSpectrum(double[] mzArr, double[] intensityArr, int scanNum) : base(mzArr, intensityArr, scanNum)
        {
        }

        public ProductSpectrum(ICollection<Peak> peaks, int scanNum): base(peaks, scanNum)
        {
        }

        public ActivationMethod ActivationMethod { get; set; }
        public IsolationInfo IsolationInfo { get; set; }

        public void SetMsLevel(int msLevel)
        {
            MsLevel = msLevel;
        }

        public override void WriteTo(StreamWriter writer)
        {
            writer.WriteLine("BEGIN IONS");
            writer.WriteLine("SCANS={0}", ScanNum);
            writer.WriteLine("_MSLEVEL={0}", MsLevel);
            writer.WriteLine("_ACTIVATION={0}", ActivationMethod.ToString());
            writer.WriteLine("_ISOTARGET={0}", IsolationInfo.IsolationWindowTargetMz);
            writer.WriteLine("_ISOLOWER={0}", IsolationInfo.IsolationWindowLowerOffset);
            writer.WriteLine("_ISOUPPER={0}", IsolationInfo.IsolationWindowUpperOffset);
            foreach (var peak in Peaks)
            {
                writer.WriteLine("{0}\t{1}", Convert.ToSingle(peak.Mz), Convert.ToSingle(peak.Intensity));
            }
            writer.WriteLine("END IONS");
        }


    }
}
