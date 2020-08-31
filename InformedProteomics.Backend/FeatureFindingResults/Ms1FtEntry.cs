using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using CsvHelper;
using CsvHelper.Configuration;

namespace InformedProteomics.Backend.FeatureFindingResults
{
    /// <summary>
    /// Data that is written to a ms1ft file
    /// </summary>
    public class Ms1FtEntry
    {
        // Ignore Spelling: Xic

        /// <summary>
        /// Feature ID
        /// </summary>
        public int FeatureId { get; set; }

        /// <summary>
        /// The first scan containing the feature
        /// </summary>
        public int MinScan { get; set; }

        /// <summary>
        /// The last scan containing the feature
        /// </summary>
        public int MaxScan { get; set; }

        /// <summary>
        /// The minimum charge seen for the feature
        /// </summary>
        public int MinCharge { get; set; }

        /// <summary>
        /// The maximum charge seen for the feature
        /// </summary>
        public int MaxCharge { get; set; }

        /// <summary>
        /// The feature monoisotopic mass
        /// </summary>
        public double MonoMass { get; set; }

        /// <summary>
        /// The representative scan number for the feature
        /// </summary>
        public int RepresentativeScan { get; set; }

        /// <summary>
        /// The representative charge number for the feature
        /// </summary>
        public int RepresentativeCharge { get; set; }

        /// <summary>
        /// The representative m/z number for the feature
        /// </summary>
        public double RepresentativeMz { get; set; }

        /// <summary>
        /// Feature abundance
        /// </summary>
        public double Abundance { get; set; }

        /// <summary>
        /// Scan number where the feature was observed with the highest intensity
        /// </summary>
        public int ApexScanNum { get; set; }

        /// <summary>
        /// The highest observed intensity for the feature
        /// </summary>
        public double ApexIntensity { get; set; }

        /// <summary>
        /// The minimum elution time for the feature
        /// </summary>
        public double MinElutionTime { get; set; }

        /// <summary>
        /// The maximum elution time for the feature
        /// </summary>
        public double MaxElutionTime { get; set; }

        /// <summary>
        /// The length of time the feature was eluting
        /// </summary>
        public double ElutionLength { get; set; }

        /// <summary>
        /// The theoretical envelope of the feature
        /// </summary>
        public string Envelope { get; set; } = "";

        /// <summary>
        /// The feature score
        /// </summary>
        public double LikelihoodRatio { get; set; }

        /// <summary>
        /// Extended data for the feature
        /// </summary>
        public Ms1FtEntryExtendedData ExtendedData { get; set; }

        /// <summary>
        /// Write the data to a ms1ft file
        /// </summary>
        /// <param name="filePath">path where the .ms1ft file should be written</param>
        /// <param name="features">features to output to the file</param>
        /// <param name="writeExtendedData">if true, the data in <see cref="ExtendedData"/> will also be output to the .ms1ft file</param>
        public static void WriteToFile(string filePath, IEnumerable<Ms1FtEntry> features, bool writeExtendedData = false)
        {
            using (var tsv = new CsvWriter(new StreamWriter(new FileStream(filePath, FileMode.Create, FileAccess.Write, FileShare.ReadWrite)), CultureInfo.InvariantCulture))
            {
                SetCsvWriterConfig(tsv.Configuration);
                if (writeExtendedData)
                {
                    tsv.Configuration.RegisterClassMap<Ms1FtEntryExtendedData.Ms1FtEntryExtendedDataMap>();
                    tsv.Configuration.RegisterClassMap<Ms1FtExtendedEntryMap>();
                }
                else
                {
                    tsv.Configuration.RegisterClassMap<Ms1FtEntryMap>();
                }
                tsv.WriteRecords(features);
            }
        }

        /// <summary>
        /// Read the data from a ms1ft file
        /// </summary>
        /// <param name="filePath">path to .ms1ft file</param>
        /// <param name="readExtendedData">if true, and file has extended data, <see cref="ExtendedData"/> will be populated from the file</param>
        /// <returns></returns>
        public static IEnumerable<Ms1FtEntry> ReadFromFile(string filePath, bool readExtendedData = false)
        {
            using (var stream = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
            {
                var hasExtended = false;
                var line = stream.ReadLine();
                if (!string.IsNullOrWhiteSpace(line) && line.IndexOf("BestEvenCharge", StringComparison.OrdinalIgnoreCase) >= 0)
                {
                    hasExtended = true;
                }
                stream.BaseStream.Seek(0, SeekOrigin.Begin);
                stream.DiscardBufferedData();
                using (var tsv = new CsvReader(stream, CultureInfo.InvariantCulture))
                {
                    SetCsvReaderConfig(tsv.Configuration);
                    if (readExtendedData && hasExtended)
                    {
                        tsv.Configuration.RegisterClassMap<Ms1FtEntryExtendedData.Ms1FtEntryExtendedDataMap>();
                        tsv.Configuration.RegisterClassMap<Ms1FtExtendedEntryMap>();
                    }
                    else
                    {
                        tsv.Configuration.RegisterClassMap<Ms1FtEntryMap>();
                    }
                    foreach (var record in tsv.GetRecords<Ms1FtEntry>())
                    {
                        if (readExtendedData && !hasExtended)
                        {
                            record.ExtendedData = new Ms1FtEntryExtendedData();
                        }
                        yield return record;
                    }
                }
            }
        }

        private static void SetCsvReaderConfig(IReaderConfiguration config)
        {
            config.Delimiter = "\t";

            //    delegate(string header, int i) { return header?.Trim().ToLower(); };
            config.PrepareHeaderForMatch = (header, i) => header?.Trim().ToLower();

            config.HeaderValidated = null;
            config.MissingFieldFound = null;
            //config.BadDataFound = null;
            config.Comment = '#';
            config.AllowComments = true;
        }

        private static void SetCsvWriterConfig(IWriterConfiguration config)
        {
            config.Delimiter = "\t";
            config.Comment = '#';
            config.AllowComments = true;
        }

        /// <summary>
        /// Class mapping <see cref="Ms1FtEntry"/> properties to text file columns, excluding extended data
        /// </summary>
        public class Ms1FtEntryMap : ClassMap<Ms1FtEntry>
        {
            /// <summary>
            /// Column count used to provide indices
            /// </summary>
            protected int ColumnCount;

            /// <summary>
            /// Column count used to provide indices
            /// </summary>
            public int GetColumnCount => ColumnCount;

            /// <summary>
            /// Constructor: Create the mapping
            /// </summary>
            [SuppressMessage("ReSharper", "VirtualMemberCallInConstructor")]
            // ReSharper disable once MemberCanBeProtected.Local
            public Ms1FtEntryMap()
            {
                Map(x => x.FeatureId).Index(ColumnCount++).Name("FeatureID");
                Map(x => x.MinScan).Index(ColumnCount++).Name("MinScan");
                Map(x => x.MaxScan).Index(ColumnCount++).Name("MaxScan");
                Map(x => x.MinCharge).Index(ColumnCount++).Name("MinCharge");
                Map(x => x.MaxCharge).Index(ColumnCount++).Name("MaxCharge");
                Map(x => x.MonoMass).Index(ColumnCount++).Name("MonoMass").TypeConverterOption.Format("0.0000");
                Map(x => x.RepresentativeScan).Index(ColumnCount++).Name("RepScan");
                Map(x => x.RepresentativeCharge).Index(ColumnCount++).Name("RepCharge");
                Map(x => x.RepresentativeMz).Index(ColumnCount++).Name("RepMz").TypeConverterOption.Format("0.0000");
                Map(x => x.Abundance).Index(ColumnCount++).Name("Abundance").TypeConverterOption.Format("0.00");
                Map(x => x.ApexScanNum).Index(ColumnCount++).Name("ApexScanNum");
                Map(x => x.ApexIntensity).Index(ColumnCount++).Name("ApexIntensity").TypeConverterOption.Format("0.00");
                Map(x => x.MinElutionTime).Index(ColumnCount++).Name("MinElutionTime").TypeConverterOption.Format("0.000");
                Map(x => x.MaxElutionTime).Index(ColumnCount++).Name("MaxElutionTime").TypeConverterOption.Format("0.000");
                Map(x => x.ElutionLength).Index(ColumnCount++).Name("ElutionLength").TypeConverterOption.Format("0.000");
                Map(x => x.Envelope).Index(ColumnCount++).Name("Envelope");
                Map(x => x.LikelihoodRatio).Index(ColumnCount++).Name("LikelihoodRatio").TypeConverterOption.Format("0.0000");
            }
        }

        /// <summary>
        /// Class mapping <see cref="Ms1FtEntry"/> with <see cref="Ms1FtEntryExtendedData"/> properties to text file columns
        /// </summary>
        public class Ms1FtExtendedEntryMap : Ms1FtEntryMap
        {
            /// <summary>
            /// Constructor: Create the mapping, inheriting the base mapping
            /// </summary>
            [SuppressMessage("ReSharper", "VirtualMemberCallInConstructor")]
            public Ms1FtExtendedEntryMap()
            {
                //References<Ms1FtEntryExtendedData.Ms1FtEntryExtendedDataMap>(x => x.ExtendedData); // Not working, but should (in theory) provide the same result
                Map(x => x.ExtendedData.BestEvenCharge).Index(ColumnCount++).Name("BestEvenCharge");
                Map(x => x.ExtendedData.BestOddCharge).Index(ColumnCount++).Name("BestOddCharge");
                Map(x => x.ExtendedData.CorrEvenCharge).Index(ColumnCount++).Name("CorrEvenCharge").TypeConverterOption.Format("0.000");
                Map(x => x.ExtendedData.CorrOddCharge).Index(ColumnCount++).Name("CorrOddCharge").TypeConverterOption.Format("0.000");
                Map(x => x.ExtendedData.IntensityEvenCharge).Index(ColumnCount++).Name("IntensityEvenCharge").TypeConverterOption.Format("0.000");
                Map(x => x.ExtendedData.IntensityOddCharge).Index(ColumnCount++).Name("IntensityOddCharge").TypeConverterOption.Format("0.000");
                Map(x => x.ExtendedData.SummedCorrEvenCharge).Index(ColumnCount++).Name("SummedCorrEvenCharge").TypeConverterOption.Format("0.000");
                Map(x => x.ExtendedData.SummedCorrOddCharge).Index(ColumnCount++).Name("SummedCorrOddCharge").TypeConverterOption.Format("0.000");
                Map(x => x.ExtendedData.SummedIntensityEvenCharge).Index(ColumnCount++).Name("SummedIntensityEvenCharge").TypeConverterOption.Format("0.000");
                Map(x => x.ExtendedData.SummedIntensityOddCharge).Index(ColumnCount++).Name("SummedIntensityOddCharge").TypeConverterOption.Format("0.000");
                Map(x => x.ExtendedData.XicCorrBetCharges1).Index(ColumnCount++).Name("XicCorrBetCharges1").TypeConverterOption.Format("0.000");
                Map(x => x.ExtendedData.XicCorrBetCharges2).Index(ColumnCount++).Name("XicCorrBetCharges2").TypeConverterOption.Format("0.000");
                Map(x => x.ExtendedData.AbundanceRatioEvenCharge).Index(ColumnCount++).Name("AbundanceRatioEvenCharge").TypeConverterOption.Format("0.000");
                Map(x => x.ExtendedData.AbundanceRatioOddCharge).Index(ColumnCount++).Name("AbundanceRatioOddCharge").TypeConverterOption.Format("0.000");
            }
        }
    }
}
