using System.Collections.Generic;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Centroid spectra (copied from Skyline)
    /// </summary>
    public class Centroider
    {
        /// <summary>
        /// List of m/zs to be centroided
        /// </summary>
        public IList<double> Mzs { get; }

        /// <summary>
        /// List of intensities to be centroided
        /// </summary>
        public IList<double> Intensities { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="mzs"></param>
        /// <param name="intensities"></param>
        public Centroider(IList<double> mzs, IList<double> intensities)
        {
            Mzs = mzs;
            Intensities = intensities;
        }

        /// <summary>
        /// Get the centroided data
        /// </summary>
        /// <param name="centroidedMzs"></param>
        /// <param name="centroidedIntensities"></param>
        public void GetCentroidedData(out double[] centroidedMzs, out double[] centroidedIntensities)
        {
            var list1 = new List<double>();
            var list2 = new List<double>();
            var flag = true;
            var num1 = 0.0;
            var num2 = 0.0;
            var num3 = 0.0;
            for (var index = 0; index < Mzs.Count; ++index)
            {
                var num4 = Intensities[index];
                if (num4 < num3)
                {
                    flag = false;
                }
                else
                {
                    if (!flag)
                    {
                        if (num2 > 0.0)
                        {
                            list1.Add(num1);
                            list2.Add(num2);
                        }
                        num2 = 0.0;
                    }
                    flag = true;
                }
                var num5 = num2 + num4;
                if (num5 > 0.0)
                {
                    var num6 = Mzs[index];
                    num1 = (num1 * num2 + num6 * num4) / num5;
                    num2 = num5;
                }
                num3 = num4;
            }
            if (num2 > 0.0)
            {
                list1.Add(num1);
                list2.Add(num2);
            }
            centroidedMzs = list1.ToArray();
            centroidedIntensities = list2.ToArray();
        }
    }
}
