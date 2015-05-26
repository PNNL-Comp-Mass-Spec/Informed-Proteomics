using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using MathNet.Numerics.Statistics;

namespace InformedProteomics.Backend.Quantification
{
    public class TargetMs1FeatureMatrix: Ms1FeatureMatrix
    {
        public TargetMs1FeatureMatrix(LcMsRun run, int minScanCharge = 2, int maxScanCharge = 60, int maxThreadCount = 0)
            : base(run, minScanCharge, maxScanCharge, maxThreadCount, new MzComparerWithBinning(27), false)
        {
            _bcDistProfile = new double[NColumns];
            _bcDistProfileCharge = new double[NRows];
        }
        
        public Ms1Feature FindMs1Feature(TargetFeature feature)
        {
            _targetFeature = feature;
            SetQueryMass(_targetFeature.Mass);

            var ms1Feature = SearchTargetFeature();

            return ms1Feature;
        }
        
        public double GetMs1EvidenceScore(TargetFeature targetFeature)
        {
            SetQueryMass(targetFeature.Mass);

            var totalElutionTime = Run.GetElutionTime(Run.MaxLcScan);
            var elutionTime = Run.GetElutionTime(targetFeature.ScanNum);
            var charge = targetFeature.Charge;
            var ms1ScanNums = Run.GetMs1ScanVector();

            var minElutionTime = elutionTime - totalElutionTime * 0.003;
            var maxElutionTime = elutionTime + totalElutionTime * 0.003;

            var ms1ScanIndex = Array.BinarySearch(ms1ScanNums, targetFeature.ScanNum);
            if (ms1ScanIndex < 0) ms1ScanIndex = ~ms1ScanIndex; // next Ms1 scan num

            var bcDistances = new List<double>();
            var correlations = new List<double>();
            var enelopes = new List<ObservedEnvelope>();
            var maxSearchScans = (int)Math.Max(ms1ScanNums.Length - ms1ScanIndex + 1, ms1ScanIndex);

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
                    var obsEnv = new ObservedEnvelope(charge - MinScanCharge, col, observedPeaks, TheoreticalEnvelope);
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


        private Tuple<int, int> GetElutionWindow(int col, double halfWindowSize)
        {
            var minCol = col;
            var maxCol = col;
            var ms1ScanNums = Run.GetMs1ScanVector();
            
            for (var j = col; j >= 0; j--)
            {
                if (j < Cols.First() || j > Cols.Last()) continue;
                if (Run.GetElutionTime(_targetFeature.ScanNum) - Run.GetElutionTime(ms1ScanNums[j]) >
                    halfWindowSize) break;
                if (j < minCol) minCol = j;
            }

            for (var j = col + 1; j < NColumns; j++)
            {
                if (j < Cols.First() || j > Cols.Last()) continue;
                if (Run.GetElutionTime(ms1ScanNums[j]) - Run.GetElutionTime(_targetFeature.ScanNum) >
                    halfWindowSize) break;
                if (j > maxCol) maxCol = j;
            }            
            return new Tuple<int, int>(minCol, maxCol);
        }


        private static readonly SavitzkyGolaySmoother Smoother = new SavitzkyGolaySmoother(9, 2);
        private TargetFeature _targetFeature;

        private Ms1Feature SearchTargetFeature()
        {
            BuildFeatureMatrix(); // should be called first
            
            const double initialElutionPeriod = 0.5; // 1 minute
            var ms1ScanNums = Run.GetMs1ScanVector();

            var col = Array.BinarySearch(ms1ScanNums, _targetFeature.ScanNum);
            if (col < 0) col = ~col;

            var row = _targetFeature.Charge - MinScanCharge;
            var tempWindow = GetElutionWindow(col, initialElutionPeriod);
            var minCol = tempWindow.Item1;
            var maxCol = tempWindow.Item2;
            
            // 1) determine charge state range
            var minMaxRow = FindMinMaxRow(row, minCol, maxCol);
            var minRow = minMaxRow.Item1;
            var maxRow = minMaxRow.Item2;

            // 2) determine elution period
            var minMaxCol = FindMinMaxCol(col, minRow, maxRow);
            minCol = minMaxCol.Item1;
            maxCol = minMaxCol.Item2;
            
            // 3) determine abundance
            var xic = new double[NColumns];
            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (DistanceMap[i][j] < 0.3 || CorrelationMap[i][j] > 0.6)
                    {
                        for (var k = 0; k < TheoreticalEnvelope.Count; k++)
                        {
                            if (FeatureMatrix[i][j][k] != null) xic[j] += FeatureMatrix[i][j][k].Intensity;
                        }
                    }
                }
            }

            var smoothedXic = Smoother.Smooth(xic);
            var abunddance = 0d;
            for (var k = 0; k < smoothedXic.Length - 1; k++)
            {
                var centerIntensity = 0.5*(Math.Max(0, smoothedXic[k]) + Math.Max(0, smoothedXic[k + 1]));

                if (!(centerIntensity > 0)) continue;
                var timeInterval = Run.GetElutionTime(ms1ScanNums[k + 1]) - Run.GetElutionTime(ms1ScanNums[k]);

                abunddance += centerIntensity*timeInterval;
            }

            var feature = new Ms1Feature(Run, minRow + MinScanCharge, maxRow + MinScanCharge,
                                        ms1ScanNums[minCol], ms1ScanNums[maxCol], abunddance, 
                                        _targetFeature.Mass, _targetFeature.Charge, 0d,_targetFeature.ScanNum);

            var envelopes = GetEnvelopeList(minRow, maxRow, minCol, maxCol);
            feature.EnvelopeList = envelopes;

            /*
            var sb = new StringBuilder(1024);
            sb.Append(0);
            for (var j = minCol; j <= maxCol; j++)
            {
                sb.Append("\t");
                sb.Append(_bcDistProfile[j]);
            }
            sb.Append("\n");
            for (var i = minRow; i <= maxRow; i++)
            {
                sb.Append(_bcDistProfileCharge[i]);
                for (var j = minCol; j <= maxCol; j++)
                {
                    sb.Append("\t");
                    sb.Append(DistanceMap[i][j]);
                }
                sb.Append("\n");
            }
            feature.Desc = sb.ToString();
            */
            return feature;
        }

        private List<ObservedEnvelope> GetEnvelopeList(int minRow, int maxRow, int minCol, int maxCol)
        {
            var envelopes = new List<ObservedEnvelope>();
            var bcDistList = new List<double>();
            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!(AccurateMass[i][j] > 0)) continue;
                    var envelope = new ObservedEnvelope(i, j, FeatureMatrix[i][j], TheoreticalEnvelope);
                    envelopes.Add(envelope);
                    bcDistList.Add(DistanceMap[i][j]);
                }
            }

            return envelopes;
            /*
            //var bcDistSample = bcDistList.Where(d => d < bcDistList.Median()).ToList();
            //var meanStd = bcDistSample.MeanStandardDeviation();
            var meanStd = bcDistList.MeanStandardDeviation();
            var bcCutoff = Math.Max(meanStd.Item1 + 2 * meanStd.Item2, 0.05);
            bcCutoff = Math.Min(bcCutoff, 0.5);

            var goodEnvelopes = new List<ObservedEnvelope>();
            for (var k = 0; k < bcDistList.Count; k++)
            {
                if (bcDistList[k] < bcCutoff) goodEnvelopes.Add(envelopes[k]);
            }
            return goodEnvelopes;
             */
        }


        private readonly double[] _bcDistProfile;
        private readonly double[] _bcDistProfileCharge;
        private Tuple<int, int> FindMinMaxCol(int seedCol, int minRow, int maxRow)
        {
            const double searchStopBcThreshold = 0.2;
            var minColTemp = seedCol;
            var maxColTemp = seedCol;

            for (var j = seedCol; j <= Cols.Last(); j++)
            {
                var summedBcDist = GetBestSummedEnvelope(minRow, maxRow, j, j);
                if (summedBcDist > searchStopBcThreshold && j > seedCol && _bcDistProfile[j - 1] > searchStopBcThreshold) break;
                if (j > maxColTemp) maxColTemp = j;
                _bcDistProfile[j] = summedBcDist;
            }

            for (var j = seedCol - 1; j >= Cols.First(); j--)
            {
                var summedBcDist = GetBestSummedEnvelope(minRow, maxRow, j, j);
                if (summedBcDist > searchStopBcThreshold && j < seedCol - 1 && _bcDistProfile[j+1] > searchStopBcThreshold) break;
                if (j < minColTemp) minColTemp = j;
                _bcDistProfile[j] = summedBcDist;
            }

            // sampling bc distances around local minimum using 1 min. window size
            var centerCol = seedCol;
            var tempWindow = GetElutionWindow(centerCol, 0.5);
            var bcDistSample = new List<double>();
            var bcDistSample2 = new List<double>();
            for (var j = tempWindow.Item1; j <= tempWindow.Item2; j++)
            {
                bcDistSample.Add(_bcDistProfile[j]);
            }
            bcDistSample.Sort();
            for(var i = 0; i < 20 && i < bcDistSample.Count; i++) bcDistSample2.Add(bcDistSample[i]);

            var meanStd = bcDistSample2.MeanStandardDeviation();
            var bcCutoff = Math.Max(meanStd.Item1 + 2 * meanStd.Item2, 0.01);
            bcCutoff = Math.Min(bcCutoff, 0.1);

            var maxCol = centerCol;
            var minCol = centerCol;
            for (var j = centerCol + 1; j <= maxColTemp; j++)
            {
                if (_bcDistProfile[j] > bcCutoff && _bcDistProfile[j + 1] > bcCutoff) break;
                maxCol = j;
            }
            for (var j = centerCol - 1; j >= minColTemp; j--)
            {
                if (_bcDistProfile[j] > bcCutoff && _bcDistProfile[j - 1] > bcCutoff) break;
                minCol = j;
            }

            return new Tuple<int, int>(minCol, maxCol);
        }
        
        private Tuple<int, int> FindMinMaxRow(int seedRow, int minCol, int maxCol)
        {
            const double searchStopBcThreshold = 0.2;
            var minRow = seedRow;
            var maxRow = seedRow;

            var bcDistSample = new List<double>();
            foreach(var i in Rows)
            {
                var summedBcDist = GetBestSummedEnvelope(i, i, minCol, maxCol);
                _bcDistProfileCharge[i] = summedBcDist;
                if (summedBcDist < searchStopBcThreshold) bcDistSample.Add(summedBcDist);
            }
            
            var meanStd = bcDistSample.MeanStandardDeviation();
            var bcCutoff = Math.Max(meanStd.Item1 + 2 * meanStd.Item2, 0.05);

            foreach (var i in Rows)
            {
                if (_bcDistProfileCharge[i] > bcCutoff) continue;
                if (i < minRow) minRow = i;
                if (i > maxRow) maxRow = i;
            }
         
            return new Tuple<int, int>(minRow, maxRow);
        }

        private double GetBestSummedEnvelope(int minRow, int maxRow, int minCol, int maxCol)
        {
            var summedBcDist = 10.0d;
            var cellList = new List<Tuple<double, int, int>>();
            for (var i = minRow; i <= maxRow; i++)
            {
                for (var j = minCol; j <= maxCol; j++)
                {
                    if (!(AccurateMass[i][j] > 0)) continue;
                    var bcDist = DistanceMap[i][j];
                    var cell = new Tuple<double, int, int>(bcDist, i, j);
                    cellList.Add(cell);
                }
            }

            if (cellList.Count < 1) return summedBcDist;
            if (cellList.Count == 1) return cellList[0].Item1;

            var tempEnvelope = new double[TheoreticalEnvelope.Count];
            var summedEnvelope = new double[TheoreticalEnvelope.Count];
            
            foreach (var cell in cellList.OrderBy(c => c.Item1))
            {
                FeatureMatrix[cell.Item2][cell.Item3].SumEnvelopeTo(tempEnvelope);
                var tempBcDist = TheoreticalEnvelope.GetBhattacharyyaDistance(tempEnvelope);
                if (tempBcDist > summedBcDist)
                {
                    Array.Copy(summedEnvelope, tempEnvelope, tempEnvelope.Length);
                }
                else
                {
                    summedBcDist = tempBcDist;
                    Array.Copy(tempEnvelope, summedEnvelope, tempEnvelope.Length);
                }
            }
            return summedBcDist;
        }
        

    }
}

