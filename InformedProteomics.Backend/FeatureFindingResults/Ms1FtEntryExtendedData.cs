using CsvHelper.Configuration;

namespace InformedProteomics.Backend.FeatureFindingResults
{
    /// <summary>
    /// Extended data that may be written to a ms1ft file
    /// </summary>
    public class Ms1FtEntryExtendedData
    {
        /// <summary>
        /// The best even charge for the feature
        /// </summary>
        public int BestEvenCharge { get; set; }

        /// <summary>
        /// The best odd charge for the feature
        /// </summary>
        public int BestOddCharge { get; set; }

        /// <summary>
        /// The best correlation score across the charge (even charge)
        /// </summary>
        public double CorrEvenCharge { get; set; }

        /// <summary>
        /// The best correlation score across the charge (odd charge)
        /// </summary>
        public double CorrOddCharge { get; set; }

        /// <summary>
        /// The best intensity score across the charge (even charge)
        /// </summary>
        public double IntensityEvenCharge { get; set; }

        /// <summary>
        /// The best intensity score across the charge (odd charge)
        /// </summary>
        public double IntensityOddCharge { get; set; }

        /// <summary>
        /// The envelope correlation score across the charge (even charge)
        /// </summary>
        public double SummedCorrEvenCharge { get; set; }

        /// <summary>
        /// The envelope correlation score across the charge (odd charge)
        /// </summary>
        public double SummedCorrOddCharge { get; set; }

        /// <summary>
        /// The envelope intensity score across the charge (even charge)
        /// </summary>
        public double SummedIntensityEvenCharge { get; set; }

        /// <summary>
        /// The envelope intensity score across the charge (odd charge)
        /// </summary>
        public double SummedIntensityOddCharge { get; set; }

        /// <summary>
        /// The XIC correlation between best charges (even charge)
        /// </summary>
        public double XicCorrBetCharges1 { get; set; }

        /// <summary>
        /// The XIC correlation between best charges (odd charge)
        /// </summary>
        public double XicCorrBetCharges2 { get; set; }

        /// <summary>
        /// The abundance distribution across the charge (even charge)
        /// </summary>
        public double AbundanceRatioEvenCharge { get; set; }

        /// <summary>
        /// The abundance distribution across the charge (odd charge)
        /// </summary>
        public double AbundanceRatioOddCharge { get; set; }

        /// <summary>
        /// Class mapping <see cref="Ms1FtEntryExtendedData"/> properties to text file columns
        /// </summary>
        public class Ms1FtEntryExtendedDataMap : ClassMap<Ms1FtEntryExtendedData>
        {
            /// <summary>
            /// Column count used to provide indices
            /// </summary>
            protected int ColumnCount;

            /// <summary>
            /// Constructor: Create the mapping
            /// </summary>
            public Ms1FtEntryExtendedDataMap()
            {
                ColumnCount = new Ms1FtEntry.Ms1FtEntryMap().GetColumnCount;

                Map(x => x.BestEvenCharge).Index(ColumnCount++).Name("BestEvenCharge");
                Map(x => x.BestOddCharge).Index(ColumnCount++).Name("BestOddCharge");
                Map(x => x.CorrEvenCharge).Index(ColumnCount++).Name("CorrEvenCharge").TypeConverterOption.Format("0.000");
                Map(x => x.CorrOddCharge).Index(ColumnCount++).Name("CorrOddCharge").TypeConverterOption.Format("0.000");
                Map(x => x.IntensityEvenCharge).Index(ColumnCount++).Name("IntensityEvenCharge").TypeConverterOption.Format("0.000");
                Map(x => x.IntensityOddCharge).Index(ColumnCount++).Name("IntensityOddCharge").TypeConverterOption.Format("0.000");
                Map(x => x.SummedCorrEvenCharge).Index(ColumnCount++).Name("SummedCorrEvenCharge").TypeConverterOption.Format("0.000");
                Map(x => x.SummedCorrOddCharge).Index(ColumnCount++).Name("SummedCorrOddCharge").TypeConverterOption.Format("0.000");
                Map(x => x.SummedIntensityEvenCharge).Index(ColumnCount++).Name("SummedIntensityEvenCharge").TypeConverterOption.Format("0.000");
                Map(x => x.SummedIntensityOddCharge).Index(ColumnCount++).Name("SummedIntensityOddCharge").TypeConverterOption.Format("0.000");
                Map(x => x.XicCorrBetCharges1).Index(ColumnCount++).Name("XicCorrBetCharges1").TypeConverterOption.Format("0.000");
                Map(x => x.XicCorrBetCharges2).Index(ColumnCount++).Name("XicCorrBetCharges2").TypeConverterOption.Format("0.000");
                Map(x => x.AbundanceRatioEvenCharge).Index(ColumnCount++).Name("AbundanceRatioEvenCharge").TypeConverterOption.Format("0.000");
                Map(x => x.AbundanceRatioOddCharge).Index(ColumnCount++).Name("AbundanceRatioOddCharge").TypeConverterOption.Format("0.000");
            }
        }
    }
}
