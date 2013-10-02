using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.TopDownScoring
{
    public class TopDownScorer //TODO consider binning .... 
    {
        public static int MinCharge = 3;
        public static int MaxCharge = 60; 

        private readonly LcMsRun _run;
        private readonly Tolerance _tolerance;
        private readonly double[] _isotopeEnvelope;
        private readonly int _maxIsotopeIndex;
        private readonly Composition _composition;
        private readonly double[][][] _xicArray; // charge, num isotope, scan number index
        private readonly double[][][] _smoothedXicArray; // charge, num isotope, scan number index
        private readonly SubScoreFactory _factory;
         
        private const int IsotopeWindowsOneSideSizeInDa = 2; // windowsOneSideLength*2 + 1 isotopes will be considered

        public TopDownScorer(Composition proteinComposition, LcMsRun run, Tolerance tolernace, SubScoreFactory factory)
        {
            _run = run;
            _composition = proteinComposition;
            _tolerance = tolernace;
            _factory = factory;

            _maxIsotopeIndex = _composition.GetMostAbundantIsotopeZeroBasedIndex();
            _isotopeEnvelope = new double[IsotopeWindowsOneSideSizeInDa * 2 + 1];
            var thorethicalIsotopeEnvelope = _composition.GetIsotopomerEnvelop();
            for (var k = -IsotopeWindowsOneSideSizeInDa; k <= IsotopeWindowsOneSideSizeInDa; k++)
            {
                _isotopeEnvelope[k + IsotopeWindowsOneSideSizeInDa] = _maxIsotopeIndex + k>=0 ? thorethicalIsotopeEnvelope[_maxIsotopeIndex + k] : 0;
            }
            _xicArray = GetXicArray();
            _smoothedXicArray = GetSmoothedXicArray();
        }

        private Tuple<int, int> GetScanNumberRangeForCharge(int charge, int scanNumber) // if scanNumber >=0 scanNumber should be included in the range
        {
            var s = _smoothedXicArray[charge - MinCharge][IsotopeWindowsOneSideSizeInDa];
            var maxS = .0;
            var maxIndex = 0;

            if (scanNumber < 0)
            {
                for (var i = 0; i < s.Length; i++)
                {
                    if (maxS < s[i])
                    {
                        maxS = s[i];
                        maxIndex = i;
                    }
                }
            }
            else maxIndex = scanNumber;

            var l = Math.Max(0, maxIndex - 1);
            while (l > Math.Max(maxIndex - 7, 0))
            {
                if (l < maxIndex - 3 && s[l] < maxS * .3) break;
                l--;
            }

            var r = Math.Min(s.Length - 1, maxIndex + 1);
            while (r < Math.Min(maxIndex + 7, s.Length - 1))
            {
                if (r > maxIndex + 3 && s[r] > maxS * .3) break;
                r++;
            }
            return new Tuple<int, int>(l,r);
        }

        private bool filterXic(double[][] truncatedXics1, double[][] truncatedXics2)
        {
            const double threshold = 0.3; // TODO
            var abundantXic1 = truncatedXics1[IsotopeWindowsOneSideSizeInDa];
            var abundantXic2 = truncatedXics2[IsotopeWindowsOneSideSizeInDa];
            //Console.WriteLine(abundantXic1.Length + " " + abundantXic2.Length);

            //for (var i = 0; i < abundantXic1.Length; i++)
            //{
            //    Console.WriteLine("\t" + i + " " + abundantXic1[i] + " " + abundantXic2[i]);
            //}

            var corr = SimpleMath.GetCorrelation(abundantXic1, abundantXic2);
            //Console.WriteLine(" corr " + corr);
            return corr < threshold;
        }


        private int GetIsotopeCorrelationIntegerRawScore(double[][] truncatedXics)
        {
            const double threshold = 0.3; // TODO

            var abundantXic = truncatedXics[IsotopeWindowsOneSideSizeInDa];
            var apexIndex = 0;
            var max = .0;

            for (var i = 0; i < abundantXic.Length; i++)
            {
                if (!(max < abundantXic[i])) continue;
                max = abundantXic[i];
                apexIndex = i;
            }

            var isotopePeaks = new double[truncatedXics.Length];
            //Console.WriteLine(" truncatedXics.Length " + truncatedXics.Length);
            for (var k = 0; k < truncatedXics.Length; k++)
            {
                if (k != IsotopeWindowsOneSideSizeInDa)
                {
                    var corr = SimpleMath.GetCorrelation(abundantXic, truncatedXics[k]);

                    //for (var o = 0; o < abundantXic.Length;o++ )
                    //    Console.WriteLine(abundantXic[o] + "\t" + truncatedXics[k][o]);

                    //Console.WriteLine(corr);
                    if (corr < threshold) continue;
                }
                isotopePeaks[k] = truncatedXics[k][apexIndex];
                //Console.Write(isotopePeaks[k] + " ");
            }

            var isotopeCorr = SimpleMath.GetCorrelation(isotopePeaks, _isotopeEnvelope);


            //Console.WriteLine(" corr " + isotopeCorr + " score " + CorrToInt(isotopeCorr));
            return CorrToInt(isotopeCorr); //isotopeCorr to int score 
        }

        public static int CorrToInt(double corr) // the higher the better from 1 to 5
        {
            var m = 1.0;
            var score = -1;
            for (; score < 10; score++)
            {
                if (1 - m > corr)
                    break;
                m = m * 0.75; // 
            }
            return score;
        }

        public double[][][][] GetTruncatedXics(int scanNumber) // main charge, other charge, Lc, mz  i scanNumber >= 0 Xics should contain the scanNumber
        {
            var truncatedXics = new double[MaxCharge - MinCharge + 1][][][];
            for (var charge = MinCharge; charge <= MaxCharge; charge++)
            {
                var scanNumberRange = GetScanNumberRangeForCharge(charge, scanNumber);
                //if(charge == 13) Console.WriteLine(scanNumberRange.Item1 + " " + scanNumberRange.Item2);
                truncatedXics[charge - MinCharge] = new double[MaxCharge - MinCharge + 1][][];
                //truncatedXics[charge - MinCharge] = new double[_smoothedXicArray[charge - MinCharge].Length][];
                for (var charge2 = MinCharge; charge2 <= MaxCharge; charge2++)
                {
                    truncatedXics[charge - MinCharge][charge2 - MinCharge] = new double[_smoothedXicArray[charge2 - MinCharge].Length][];
                    for (var i = 0; i < truncatedXics[charge - MinCharge][charge2 - MinCharge].Length; i++)
                    {
                        truncatedXics[charge - MinCharge][charge2 - MinCharge][i] = new double[scanNumberRange.Item2 - scanNumberRange.Item1];
                        for (var j = scanNumberRange.Item1; j < scanNumberRange.Item2; j++)
                        {
                            //Console.WriteLine(j + " " + truncatedXics[charge - MinCharge][charge2 - MinCharge][i].Length + " " + _smoothedXicArray[charge2 - MinCharge][i].Length);
                            truncatedXics[charge - MinCharge][charge2 - MinCharge][i][j - scanNumberRange.Item1] =
                                _smoothedXicArray[charge2 - MinCharge][i][j];
                        }
                    }
                }
            }
            return truncatedXics;
        }

        public int GetScore()
        {
            var truncatedXics = GetTruncatedXics(-1);
            var maxScore = -10000.0;
            
            for (var charge = MinCharge; charge <= MaxCharge; charge++)
            {
                //Console.WriteLine("\t" + charge);
                
                //if(charge == 13)
                //    foreach(var t in truncatedXics[charge-MinCharge][charge-MinCharge][IsotopeWindowsOneSideSizeInDa])
                //        Console.WriteLine(t );

                var score = .0;
                for (var charge2 = MinCharge; charge2 <= MaxCharge; charge2++)
                {
                    double diff;
                    if (charge != charge2 &&
                        filterXic(truncatedXics[charge - MinCharge][charge - MinCharge],
                                  truncatedXics[charge - MinCharge][charge2 - MinCharge]))
                    {
                        diff = _factory == null ? -.2 : _factory.GetMissingXicScore(charge2);
                    }
                    else
                    {
                        var corrIntScore = GetIsotopeCorrelationIntegerRawScore(truncatedXics[charge - MinCharge][charge2 - MinCharge]);
                        diff = _factory == null ? corrIntScore : _factory.GetXicCorrScore(charge2, corrIntScore);
                    }
                    score += diff;
                    //if (diff!= -.2) Console.WriteLine(charge2 + " " + diff);
                }
                if (score > maxScore) maxScore = score;
            }
            return (int)maxScore;
        }


        public bool[] AreXicMissing(int charge, int scanNumber) // charge is fixed, and the specified scanNumber should be in scanRange : for training
        {
            var truncatedXics = GetTruncatedXics(scanNumber);
            var missing = new bool[MaxCharge - MinCharge + 1];
            for (var charge2 = MinCharge; charge2 <= MaxCharge; charge2++)
            {
                if (charge != charge2 &&
                    filterXic(truncatedXics[charge - MinCharge][charge - MinCharge],
                                truncatedXics[charge - MinCharge][charge2 - MinCharge]))
                {
                    missing[charge2 - MinCharge] = true;
                }
            }
           
            return missing;
        }

        public int[] GetIsotopeCorrelationIntensityRawScores(int charge, int scanNumber) // charge is fixed, and the specified scanNumber should be in scanRange : for training
        {
            var truncatedXics = GetTruncatedXics(scanNumber);
            var scores = new int[MaxCharge - MinCharge + 1];
            for (var charge2 = MinCharge; charge2 <= MaxCharge; charge2++)
            {
                if (charge != charge2 &&
                    filterXic(truncatedXics[charge - MinCharge][charge - MinCharge],
                              truncatedXics[charge - MinCharge][charge2 - MinCharge]))
                {
                    //missing[charge2 - MinCharge] = true;
                }
                else
                {
                    scores[charge2 - MinCharge] = GetIsotopeCorrelationIntegerRawScore(truncatedXics[charge - MinCharge][charge2 - MinCharge]);
                }
            }

            return scores;
        }

        private double[][][] GetSmoothedXicArray()
        {
            var smoothedXicArray = new double[_xicArray.Length][][];
            const int smoothingWindowOneSideSize = 2;

            for (var i = 0; i < _xicArray.Length; i++)
            {
                smoothedXicArray[i] = new double[_xicArray[i].Length][];
                for (var j = 0; j < _xicArray[i].Length; j++)
                {
                    var sum = .0;
                    var x = _xicArray[i][j];
                    var s = new double[x.Length];
                    //var maxS = .0;
                    for (var k = -smoothingWindowOneSideSize;
                         k < x.Length + smoothingWindowOneSideSize;
                         k++)
                    {
                        var r = k + smoothingWindowOneSideSize;
                        var l = k - smoothingWindowOneSideSize;
                        if (r < x.Length) sum += x[r];
                        if (k >= 0 && k < x.Length)
                        {
                            s[k] = sum;
                            //maxS = maxS < s[k] ? s[k] : maxS;
                        }
                        if (l >= 0) sum -= x[l];
                    }
                    /*if (maxS > 0)
                    {
                        for (var k = 0; k < s.Length; k++)
                            s[k] = s[k]/maxS;
                    }*/
                    smoothedXicArray[i][j] = s;
                }
            }
            return smoothedXicArray;
        }


        private double[][][] GetXicArray()
        {
            var scanNumberIndex = new Dictionary<int, int>();
            var i = 0;
            foreach (var scanNumber in _run.GetScanNumbers(1))
            {
                scanNumberIndex[scanNumber] = i++;
            }

            var xicArray = new double[MaxCharge - MinCharge + 1][][];

            for (var charge = MinCharge; charge <= MaxCharge; charge++)
            {
                var xicArrayForThisCharge = new double[IsotopeWindowsOneSideSizeInDa * 2 + 1][];
                for (var j = 0; j < xicArrayForThisCharge.Length; j++)
                    xicArrayForThisCharge[j] = new double[scanNumberIndex.Count];
                var precursorIon = new Ion(_composition + Composition.H2O, charge);

                for (var k = - IsotopeWindowsOneSideSizeInDa; k <= IsotopeWindowsOneSideSizeInDa; k++)
                {
                    var mz = precursorIon.GetIsotopeMz(_maxIsotopeIndex + k);
                    var xic = _run.GetExtractedIonChromatogram(mz, _tolerance);
                    foreach (var xicPeak in xic)
                    {
                        xicArrayForThisCharge[k + IsotopeWindowsOneSideSizeInDa][scanNumberIndex[xicPeak.ScanNum]] = xicPeak.Intensity;
                        //if(k == 0 && charge == 13) Console.WriteLine(xicPeak.ScanNum + " " + xicPeak.Intensity);
                        
                    }
                }
               // Console.WriteLine();

                xicArray[charge - MinCharge] = xicArrayForThisCharge;
                //if (charge == 13)
                //{
                 //   foreach(var t in xicArrayForThisCharge[IsotopeWindowsOneSideSizeInDa])   
                  //      Console.WriteLine(t);
                //}

            }
            return xicArray;
        }  
    }
}
