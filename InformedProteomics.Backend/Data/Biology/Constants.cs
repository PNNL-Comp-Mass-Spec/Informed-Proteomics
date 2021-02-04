using System;

namespace InformedProteomics.Backend.Data.Biology
{
    /// <summary>
    /// Common constants
    /// </summary>
    public static class Constants
    {
        /// <summary>
        /// Convert between nominal mass and monoisotopic mass with minimal error
        /// </summary>
        public const double RescalingConstant = 0.9995;

        /// <summary>
        /// Used to convert between bin number and mass
        /// </summary>
        public const double RescalingConstantHighPrecision = 274.335215;
        //public static readonly double H2O = Annotation.Composition.H2O.GetMass();

        /// <summary>
        /// Carbon-13 isotopic mass
        /// </summary>
        public const double C13 = 13.00335483;

        /// <summary>
        /// Mass difference between Carbon-13 and Carbon-12
        /// </summary>
        public const double C13MinusC12 = C13 - 12.0;

        /// <summary>
        /// Mass of a proton
        /// </summary>
        public const double Proton = 1.00727649;

        /// <summary>
        /// Get the bin number of the supplied mass <paramref name="m"/>
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static int GetBinNum(double m)
        {
            return (int)Math.Round(m * RescalingConstant);
        }

        /// <summary>
        /// Get the high-precision bin number of the supplied mass <paramref name="m"/>
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static int GetBinNumHighPrecision(double m)
        {
            return (int)Math.Round(m * RescalingConstantHighPrecision);
        }
    }
}
