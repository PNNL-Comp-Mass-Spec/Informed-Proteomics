using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Quantification;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Ms2FeatureQuntification: Ms1FeatureMatrix
    {
        public Ms2FeatureQuntification(LcMsRun run, int minScanCharge = 2, int maxScanCharge = 60, int maxThreadCount = 0)
            : base(run, minScanCharge, maxScanCharge, maxThreadCount)
        {

        }

        

        public double GetMs1EvidenceScore(Ms2Feature ms2Feature)
        {
            SetQueryMass(ms2Feature.Mass);

            var totalElutionTime = Run.GetElutionTime(Run.MaxLcScan);
            var elutionTime = Run.GetElutionTime(ms2Feature.ScanNum);
            var charge = ms2Feature.Charge;
            var ms1ScanNums = Run.GetMs1ScanVector();

            var minElutionTime = elutionTime - totalElutionTime * 0.003;
            var maxElutionTime = elutionTime + totalElutionTime * 0.003;

            var ms1ScanIndex = Array.BinarySearch(ms1ScanNums, ms2Feature.ScanNum);
            if (ms1ScanIndex < 0) ms1ScanIndex = ~ms1ScanIndex; // next Ms1 scan num

            var bcDistances = new List<double>();
            var correlations = new List<double>();
            var enelopes    = new List<ObservedEnvelope>();
            var maxSearchScans = (int) Math.Max(ms1ScanNums.Length - ms1ScanIndex + 1, ms1ScanIndex);
            
            for (var i = 0; i <= maxSearchScans; i++)
            {
                for (var j = 0; j < 2; j++)
                {
                    if (i == 0 && j > 0) continue;

                    var col = (j < 1) ? ms1ScanIndex + i : ms1ScanIndex - i;
                    if (col >= ms1ScanNums.Length || col < 0) continue;
                    if (Run.GetElutionTime(ms1ScanNums[col]) > maxElutionTime || Run.GetElutionTime(ms1ScanNums[col]) < minElutionTime) continue;

                    var ms1Spec = Spectrums[col];
                    var observedPeaks = ms1Spec.GetAllIsotopePeaks(QueryMass, charge, TheoreticalEnvelope, MzTolerance);
                    var obsEnv = new ObservedEnvelope(charge - ChargeRange.Min, col, observedPeaks, TheoreticalEnvelope);
                    if (obsEnv.NumberOfPeaks < 3) continue;
                    var bcDist = obsEnv.GetBhattacharyyaDistance(TheoreticalEnvelope.EnvelopePdf);
                    var corrCoeff = obsEnv.GetPearsonCorrelation(TheoreticalEnvelope.Envelope);

                    if (bcDist < 0.1)
                    {
                        return 1 - bcDist;
                    }
                    if (corrCoeff > 0.85)
                    {
                        return corrCoeff;
                    }

                    bcDistances.Add(bcDist);
                    correlations.Add(corrCoeff);
                    enelopes.Add(obsEnv);                    
                }
            }

            if (bcDistances.Count < 1) return 0d;
            return correlations.Max();
        }
        
        protected override double GetBcDistTh(double normalizedElutionLen)
        {
            if (QueryMass < 15000)
            {
                /*if (normalizedElutionLen < 0.005) return 0.6;
                if (normalizedElutionLen < 0.01) return 0.4;
                if (normalizedElutionLen < 0.02) return 0.2;
                return 0.1;*/

                return 0.25;
            }
            else if (QueryMass < 25000)
            {
                //if (normalizedElutionLen < 0.005) return 1.0;
                //if (normalizedElutionLen < 0.01) return 0.5;
                //if (normalizedElutionLen < 0.02) return 0.25;
                //return 0.2;

                return 0.5;
            }
            else // > 25K 
            {
                //if (normalizedElutionLen < 0.005) return 1.2;
                //if (normalizedElutionLen < 0.01) return 0.8;
                //if (normalizedElutionLen < 0.02) return 0.3;
                //return 0.2;
                return 0.8;
            }
        }

        protected override double GetCorrTh(double normalizedElutionLen)
        {
            if (QueryMass < 15000)
            {
                //if (normalizedElutionLen < 0.005) return 0.3;
                //if (normalizedElutionLen < 0.01) return 0.4;
                //if (normalizedElutionLen < 0.02) return 0.6;
                //return 0.7;
                return 0.4;
            }
            else if (QueryMass < 25000)
            {
                //if (normalizedElutionLen < 0.005) return 0;
                //if (normalizedElutionLen < 0.01) return 0.2;
                //if (normalizedElutionLen < 0.02) return 0.4;
                //return 0.6;
                return 0.2;
            }
            else // 25K
            {
                //if (normalizedElutionLen < 0.005) return -1;
                //if (normalizedElutionLen < 0.01) return 0.1;
                //if (normalizedElutionLen < 0.02) return 0.4;
                //return 0.5;
                return -1;
            }
        }

        /*
        protected override void EvaluateCluster(ref Ms1FeatureCluster cluster)
        {
            var minCol = cluster.MinCol;
            var maxCol = cluster.MaxCol;
            var minRow = cluster.MinRow;
            var maxRow = cluster.MaxRow;
            var massTol = MzTolerance.GetToleranceAsTh(cluster.RepresentativeMass);

            cluster.ClearMember();
            for (var col = minCol; col <= maxCol; col++)
            {
                for (var row = minRow; row <= maxRow; row++)
                {
                    var mass = AccurateMass[row][col];

                    if (mass > 0 && Math.Abs(cluster.RepresentativeMass - mass) < massTol)
                    {
                        var obsEnv = new ObservedEnvelope(row, col, FeatureMatrix[row][col], _isotopeList);
                        cluster.AddMember(obsEnv);
                    }
                    else
                    {
                        var observedPeaks = Spectrums[col].GetAllIsotopePeaks(cluster.RepresentativeMass, row + ChargeRange.Min, _isotopeList, MzTolerance);
                        var obsEnv = new ObservedEnvelope(row, col, observedPeaks, _isotopeList);
                        if (obsEnv.NumberOfPeaks < 3) continue;
                        cluster.AddMember(obsEnv);
                    }
                }
            }

            if (cluster.Envelopes.Count < 1) return;
            cluster.NormalizedElutionLength = (Run.GetElutionTime(cluster.MaxScanNum) -
                                               Run.GetElutionTime(cluster.MinScanNum)) /
                                              Run.GetElutionTime(Run.MaxLcScan);
            cluster.UpdateScores(Spectrums);
            if (cluster.GoodEnvelopeCount > 0)
            {
                CalculateXicCorrelationOverTimeBetweenIsotopes(cluster);
            }
        }
        */

    }
}

