using System;
using System.Text;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Spectrum
    {
        public Spectrum(double[] mzArr, double[] intensityArr)
        {
            MzArr = mzArr;
            IntensityArr = intensityArr;
        }

        public int MsLevel
        {
            get { return _msLevel; }
            set { _msLevel = value; }
        }

        public double[] MzArr { get; private set; }
        public double[] IntensityArr { get; private set; }

        public void Display()
        {
            var sb = new StringBuilder();
            sb.Append("--------- Spectrum -----------------\n");
            for (int i = 0; i < MzArr.Length; i++)
            {
                sb.Append(MzArr[i]);
                sb.Append("\t");
                sb.Append(IntensityArr[i]);
                sb.Append("\n");
            }
            sb.Append("--------------------------- end ---------------------------------------\n");

            Console.Write(sb.ToString());
        }

        private int _msLevel = 1;
    }
}
