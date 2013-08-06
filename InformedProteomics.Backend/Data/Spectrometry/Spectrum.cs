using System;
using System.Collections.Generic;
using System.Text;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Spectrum
    {
        public Spectrum(IList<double> mzArr, IList<double> intensityArr)
        {
            Peaks = new Peak[mzArr.Count];
            for(var i=0; i<mzArr.Count; i++) Peaks[i] = new Peak(mzArr[i], intensityArr[i]);
        }

        public int MsLevel
        {
            get { return _msLevel; }
            set { _msLevel = value; }
        }

        public Peak[] Peaks { get; private set; }

        /// <summary>
        /// Finds the maximum intensity peak within the specified range
        /// </summary>
        /// <param name="mz">m/z</param>
        /// <param name="tolerance">tolerance</param>
        /// <returns>maximum intensity peak</returns>
        public Peak FindPeak(double mz, Tolerance tolerance)
        {
            var tolTh = tolerance.GetToleranceAsTh(mz);
            var minMz = mz - tolTh;
            var maxMz = mz + tolTh;

            return FindPeak(minMz, maxMz);
        }

        /// <summary>
        /// Finds the maximum intensity peak within the specified range
        /// </summary>
        /// <param name="minMz">minimum m/z</param>
        /// <param name="maxMz">maximum m/z</param>
        /// <returns>maximum intensity peak within the range</returns>
        public Peak FindPeak(double minMz, double maxMz)
        {
            var index = Array.BinarySearch(Peaks, new Peak((minMz+maxMz)/2, 0));
            if (index < 0) index = ~index;

            Peak bestPeak = null;
            var bestIntensity = 0.0;

            // go down
            var i = index;
            while (i >= 0 && i < Peaks.Length)
            {
                if (Peaks[i].Mz <= minMz) break;
                if (Peaks[i].Intensity > bestIntensity)
                {
                    bestIntensity = Peaks[i].Intensity;
                    bestPeak = Peaks[i];
                }
                --i;
            }

            // go up
            i = index + 1;
            while (i >= 0 && i < Peaks.Length)
            {
                if (Peaks[i].Mz >= maxMz) break;
                if (Peaks[i].Intensity > bestIntensity)
                {
                    bestIntensity = Peaks[i].Intensity;
                    bestPeak = Peaks[i];
                }
                ++i;
            }
            return bestPeak;
        }

        public void Display()
        {
            var sb = new StringBuilder();
            sb.Append("--------- Spectrum -----------------\n");
            foreach (var peak in Peaks)
            {
                sb.Append(peak.Mz);
                sb.Append("\t");
                sb.Append(peak.Intensity);
                sb.Append("\n");
            }
            sb.Append("--------------------------- end ---------------------------------------\n");

            Console.Write(sb.ToString());
        }

        private int _msLevel = 1;
    }
}
