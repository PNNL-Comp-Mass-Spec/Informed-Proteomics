using System;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.FeatureFindingResults;
using InformedProteomics.FeatureFinding.Clustering;
using InformedProteomics.FeatureFinding.Data;

namespace InformedProteomics.FeatureFinding.Util
{
    public static class Extensions
    {
        /// <summary>
        /// Converts a <see cref="LcMsPeakCluster"/> to a <see cref="Ms1FtEntry"/>
        /// </summary>
        /// <param name="feature"></param>
        /// <param name="featureId"></param>
        /// <returns>New Ms1FtEntry object</returns>
        public static Ms1FtEntry ToMs1FtEntry(this LcMsPeakCluster feature, int featureId = 0)
        {
            var intensity = feature.RepresentativeSummedEnvelop;
            var maxIntensity = intensity.Max();
            var sb = new StringBuilder();
            for (var i = 0; i < intensity.Length; i++)
            {
                if (i != 0)
                {
                    sb.Append(";");
                }

                sb.AppendFormat("{0},{1:0.000}", feature.TheoreticalEnvelope.Isotopes[i].Index, intensity[i] / maxIntensity);
            }

            var ms1FtEntry = new Ms1FtEntry
            {
                FeatureId = featureId,
                MinScan = feature.MinScanNum,
                MaxScan = feature.MaxScanNum,
                MinCharge = feature.MinCharge,
                MaxCharge = feature.MaxCharge,
                MonoMass = feature.RepresentativeMass,
                RepresentativeScan = feature.RepresentativeScanNum,
                RepresentativeCharge = feature.RepresentativeCharge,
                RepresentativeMz = feature.RepresentativeMz,
                Abundance = feature.Abundance,
                ApexScanNum = feature.ApexScanNum,
                ApexIntensity = feature.ApexIntensity,
                MinElutionTime = feature.MinElutionTime,
                MaxElutionTime = feature.MaxElutionTime,
                ElutionLength = feature.ElutionLength,
                Envelope = sb.ToString(),
                LikelihoodRatio = feature.Score,
                ExtendedData = feature.ToMs1FtEntryExtendedData()
            };

            return ms1FtEntry;
        }

        /// <summary>
        /// Converts a <see cref="LcMsPeakCluster"/> to a <see cref="Ms1FtEntryExtendedData"/>
        /// </summary>
        /// <param name="feature"></param>
        /// <returns>New Ms1FtEntryExtendedData object</returns>
        private static Ms1FtEntryExtendedData ToMs1FtEntryExtendedData(this LcMsPeakCluster feature)
        {
            var extended = new Ms1FtEntryExtendedData
            {
                BestEvenCharge = feature.BestCharge[LcMsPeakCluster.EvenCharge],
                BestOddCharge = feature.BestCharge[LcMsPeakCluster.OddCharge],
                CorrEvenCharge = feature.BestCorrelationScoreAcrossCharge[LcMsPeakCluster.EvenCharge],
                CorrOddCharge = feature.BestCorrelationScoreAcrossCharge[LcMsPeakCluster.OddCharge],
                IntensityEvenCharge = feature.BestIntensityScoreAcrossCharge[LcMsPeakCluster.EvenCharge],
                IntensityOddCharge = feature.BestIntensityScoreAcrossCharge[LcMsPeakCluster.OddCharge],
                SummedCorrEvenCharge = feature.EnvelopeCorrelationScoreAcrossCharge[LcMsPeakCluster.EvenCharge],
                SummedCorrOddCharge = feature.EnvelopeCorrelationScoreAcrossCharge[LcMsPeakCluster.OddCharge],
                SummedIntensityEvenCharge = feature.EnvelopeIntensityScoreAcrossCharge[LcMsPeakCluster.EvenCharge],
                SummedIntensityOddCharge = feature.EnvelopeIntensityScoreAcrossCharge[LcMsPeakCluster.OddCharge],
                XicCorrBetCharges1 = feature.XicCorrelationBetweenBestCharges[LcMsPeakCluster.EvenCharge],
                XicCorrBetCharges2 = feature.XicCorrelationBetweenBestCharges[LcMsPeakCluster.OddCharge],
                AbundanceRatioEvenCharge = feature.AbundanceDistributionAcrossCharge[LcMsPeakCluster.EvenCharge],
                AbundanceRatioOddCharge = feature.AbundanceDistributionAcrossCharge[LcMsPeakCluster.OddCharge]
            };

            return extended;
        }

        /// <summary>
        /// Converts a <see cref="LcMsFeature"/> to a <see cref="Ms1FtEntry"/>
        /// </summary>
        /// <param name="feature"></param>
        /// <param name="featureId"></param>
        /// <returns>New Ms1FtEntry object</returns>
        public static Ms1FtEntry ToMs1FtEntry(this LcMsFeature feature, int featureId = 0)
        {
            var ms1FtEntry = new Ms1FtEntry
            {
                FeatureId = featureId,
                MinScan = feature.MinScanNum,
                MaxScan = feature.MaxScanNum,
                MinCharge = feature.MinCharge,
                MaxCharge = feature.MaxCharge,
                MonoMass = feature.RepresentativeMass,
                RepresentativeScan = feature.RepresentativeScanNum,
                RepresentativeCharge = feature.RepresentativeCharge,
                RepresentativeMz = feature.RepresentativeMz,
                Abundance = feature.Abundance,
                MinElutionTime = feature.MinElutionTime,
                MaxElutionTime = feature.MaxElutionTime,
                ElutionLength = feature.ElutionLength,
                LikelihoodRatio = feature.Score
            };

            return ms1FtEntry;
        }

        /// <summary>
        /// Converts a <see cref="Ms1FtEntry"/> to a <see cref="LcMsFeature"/>
        /// </summary>
        /// <param name="ms1FtEntry"></param>
        /// <param name="totalDatasetRunTime">The total runtime of the dataset; generally the elution time of the last scan. Used to compute NET</param>
        /// <returns>New LcMsFeature object</returns>
        public static LcMsFeature ToLcMsFeature(this Ms1FtEntry ms1FtEntry, double totalDatasetRunTime = 0)
        {
            var repCharge = ms1FtEntry.RepresentativeCharge != 0 ? ms1FtEntry.RepresentativeCharge : (int)Math.Round(0.5 * (ms1FtEntry.MinCharge + ms1FtEntry.MaxCharge));
            var repMz = !ms1FtEntry.RepresentativeMz.Equals(0) ? ms1FtEntry.RepresentativeMz : (ms1FtEntry.MonoMass / repCharge) + Constants.Proton;
            var repScanNum = ms1FtEntry.RepresentativeScan > 0 ? ms1FtEntry.RepresentativeScan : ms1FtEntry.MinScan;

            var minNet = 0.0;
            var maxNet = 0.0;
            if (totalDatasetRunTime > 0)
            {
                minNet = ms1FtEntry.MinElutionTime / totalDatasetRunTime;
                maxNet = ms1FtEntry.MaxElutionTime / totalDatasetRunTime;
            }

            var feature = new LcMsFeature(ms1FtEntry.MonoMass, repCharge, repMz, repScanNum, ms1FtEntry.Abundance, ms1FtEntry.MinCharge, ms1FtEntry.MaxCharge, ms1FtEntry.MinScan, ms1FtEntry.MaxScan, ms1FtEntry.MinElutionTime, ms1FtEntry.MaxElutionTime, minNet, maxNet)
            {
                FeatureId = ms1FtEntry.FeatureId,
                Score = ms1FtEntry.LikelihoodRatio,
            };

            return feature;
        }
    }
}
