using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.MathAndStats;

namespace InformedProteomics.TopDown.Scoring
{
    public class TopDownScorer //TODO consider binning ....
    {
        // Ignore Spelling: readonly, Xic, Xics, foreach, iso, diff

        public static int MinCharge = 8;
        public static int MaxCharge = 25;

        private readonly LcMsRun _run;
        private readonly Tolerance _tolerance;
        private readonly double[] _isotopeEnvelope; // from _minIsotopeIndex

        private readonly Composition _proteinCompositionPlusWater;
        private readonly double[][][] _xicArray; // charge, num isotope (from _minIsotopeIndex) , scan number index
        private readonly double[][][] _smoothedXicArray; // charge, num isotope, scan number index

        // [Obsolete("This functionality of this class is not implemented")]
        // private readonly SubScoreFactory _factory;

        private readonly int _minIsotopeIndex;
        private readonly int _maxIntensityIsotopeIndex;

        private const int NumberAfterMaxIsotopeIndex = 3; //
        private const float MinIsotopeIntensity = 0.2f;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="proteinComposition"></param>
        /// <param name="run"></param>
        /// <param name="tolerance"></param>
        /// <param name="factory"></param>
        [Obsolete("Use the version of TopDownScore that does not take factory")]
        public TopDownScorer(Composition proteinComposition, LcMsRun run, Tolerance tolerance, SubScoreFactory factory) : this(proteinComposition, run, tolerance)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="proteinComposition"></param>
        /// <param name="run"></param>
        /// <param name="tolerance"></param>
        public TopDownScorer(Composition proteinComposition, LcMsRun run, Tolerance tolerance)
        {
            _run = run;
            _proteinCompositionPlusWater = proteinComposition + Composition.H2O;
            _tolerance = tolerance;

            _maxIntensityIsotopeIndex = _proteinCompositionPlusWater.GetMostAbundantIsotopeZeroBasedIndex();
            var theoreticalIsotopeEnvelope = _proteinCompositionPlusWater.GetIsotopomerEnvelopeRelativeIntensities();

            _minIsotopeIndex = 0;
            for (var i = 0; i < theoreticalIsotopeEnvelope.Length; i++)
            {
                if (!(theoreticalIsotopeEnvelope[i] > MinIsotopeIntensity))
                {
                    continue;
                }

                _minIsotopeIndex = i;
                break;
            }

            _isotopeEnvelope = new double[Math.Min(_maxIntensityIsotopeIndex + NumberAfterMaxIsotopeIndex, theoreticalIsotopeEnvelope.Length) - _minIsotopeIndex];
            for (var k = 0; k < _isotopeEnvelope.Length; k++)
            {
                _isotopeEnvelope[k] = theoreticalIsotopeEnvelope[k + _minIsotopeIndex];
            }

            /*foreach (var iso in theoreticalIsotopeEnvelope)
            {
                Console.WriteLine(iso);
            }
            Console.WriteLine();
            foreach (var iso in _isotopeEnvelope)
            {
                Console.WriteLine(iso);
            }

            System.Environment.Exit(1);
            */
            _xicArray = GetXicArray();
            _smoothedXicArray = GetSmoothedXicArray();
        }

        private Tuple<int, int> GetScanNumberRangeForCharge(int charge, int scanNumber) // if scanNumber >=0 scanNumber should be included in the range
        {
            var s = _smoothedXicArray[charge - MinCharge][_maxIntensityIsotopeIndex - _minIsotopeIndex];
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
            else
            {
                maxIndex = scanNumber;
            }

            var l = Math.Max(0, maxIndex - 1);
            while (l > Math.Max(maxIndex - 7, 0))
            {
                if (l < maxIndex - 3 && s[l] < maxS * .3)
                {
                    break;
                }

                l--;
            }

            var r = Math.Min(s.Length - 1, maxIndex + 1);
            while (r < Math.Min(maxIndex + 7, s.Length - 1))
            {
                if (r > maxIndex + 3 && s[r] > maxS * .3)
                {
                    break;
                }

                r++;
            }
            return new Tuple<int, int>(l, r);
        }

        private bool FilterXic(IReadOnlyList<double[]> truncatedXics1, IReadOnlyList<double[]> truncatedXics2)
        {
            const double threshold = 0.5; // TODO
            var abundantXic1 = truncatedXics1[_maxIntensityIsotopeIndex - _minIsotopeIndex];
            var abundantXic2 = truncatedXics2[_maxIntensityIsotopeIndex - _minIsotopeIndex];
            //Console.WriteLine(abundantXic1.Length + " " + abundantXic2.Length);

            var corr = SimpleMath.GetCorrelation(abundantXic1, abundantXic2);

            /*            if (corr > 0)
                        {
                                for (var i = 0; i < abundantXic1.Length; i++)
                                {
                                    Console.WriteLine("\t" + i + "\t" + abundantXic1[i] + "\t" + abundantXic2[i]);
                                }

                                Console.WriteLine(" corr " + corr + "\n__________________\n");
                            }
                        */
            return corr < threshold;
        }

        private int GetIsotopeCorrelationIntegerRawScore(IReadOnlyList<double[]> truncatedXics)
        {
            const double threshold = 0.5; // TODO

            var abundantXic = truncatedXics[_maxIntensityIsotopeIndex - _minIsotopeIndex];
            var apexIndex = 0;
            var max = .0;

            for (var i = 0; i < abundantXic.Length; i++)
            {
                if (!(max < abundantXic[i]))
                {
                    continue;
                }

                max = abundantXic[i];
                apexIndex = i;
            }

            var isotopePeaks = new double[truncatedXics.Count];
            //Console.WriteLine(" truncatedXics.Length " + truncatedXics.Length);
            for (var k = 0; k < truncatedXics.Count; k++)
            {
                if (k != _maxIntensityIsotopeIndex - _minIsotopeIndex)
                {
                    var corr = SimpleMath.GetCorrelation(abundantXic, truncatedXics[k]);

                    //for (var o = 0; o < abundantXic.Length;o++ )
                    //    Console.WriteLine(abundantXic[o] + "\t" + truncatedXics[k][o]);

                    //Console.WriteLine(corr);
                    if (corr < threshold)
                    {
                        continue;
                    }
                }

                isotopePeaks[k] = truncatedXics[k][apexIndex];
                //
            }
            //for (var k = 0; k < truncatedXics.Length; k++)
            //{
            //   Console.WriteLine(k + "\t" + isotopePeaks[k] + "\t" + _isotopeEnvelope[k] + "\t" + apexIndex);
            //}

            var isotopeCorr = SimpleMath.GetCorrelation(isotopePeaks, _isotopeEnvelope);

            //Console.WriteLine(" corr " + isotopeCorr + " score " + CorrToInt(isotopeCorr));
            return CorrToInt(isotopeCorr); //isotopeCorr to int score
        }

        public static int CorrToInt(double corr) // the higher the better from 1 to 5
        {
            var m = 1.0;
            var score = -2;
            for (; score < 10; score++)
            {
                if (m <= -corr + 1)
                {
                    break;
                }

                m *= 0.6; // the higher the higher the score
            }
            return score;
        }

        public double[][][][] GetTruncatedXics(int scanNumber) // main charge, other charge, LC, mz  i scanNumber >= 0 Xics should contain the scanNumber
        {
            var truncatedXics = new double[MaxCharge - MinCharge + 1][][][];
            for (var charge = MinCharge; charge <= MaxCharge; charge++)
            {
                var scanNumberRange = GetScanNumberRangeForCharge(charge, scanNumber);
                //Console.WriteLine(charge + "\t" + scanNumberRange.Item1 + "\t" + scanNumberRange.Item2);
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
                            //.WriteLine(j + " " + truncatedXics[charge - MinCharge][charge2 - MinCharge][i].Length + " " + _smoothedXicArray[charge2 - MinCharge][i].Length);
                            truncatedXics[charge - MinCharge][charge2 - MinCharge][i][j - scanNumberRange.Item1] =
                                _smoothedXicArray[charge2 - MinCharge][i][j];
                            //  if(charge == 15 && i == _maxIntensityIsotopeIndex - _minIsotopeIndex)
                            //      Console.WriteLine(charge2 + "\t" + i + "\t" + j + "\t" + _smoothedXicArray[charge2 - MinCharge][i][j]);
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
                //    foreach(var t in truncatedXics[charge-MinCharge][charge-MinCharge][NumberAfterMaxIsotopeIndex])
                //        Console.WriteLine(t );

                var score = .0;
                for (var charge2 = MinCharge; charge2 <= MaxCharge; charge2++)
                {
                    double diff;
                    if (charge != charge2 &&
                        FilterXic(truncatedXics[charge - MinCharge][charge - MinCharge],
                                  truncatedXics[charge - MinCharge][charge2 - MinCharge]))
                    {
                        // diff = _factory?.GetMissingXicScore(charge2) ?? -.3;
                        diff = -.3;

                        //  Console.WriteLine(charge + "\t" + charge2 + "\t" + "filtered" + "\t**");
                    }
                    else
                    {
                        var corrIntScore = GetIsotopeCorrelationIntegerRawScore(truncatedXics[charge - MinCharge][charge2 - MinCharge]);
                        // diff = _factory?.GetXicCorrScore(charge2, corrIntScore) ?? corrIntScore;
                        diff = corrIntScore;

                        //diff = 1;
                        // Console.WriteLine(charge + "\t" + charge2 + "\t" + diff + "\t**");
                    }

                    score += diff;
                }
                if (score > maxScore)
                {
                    maxScore = score;
                }
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
                    FilterXic(truncatedXics[charge - MinCharge][charge - MinCharge],
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
                    FilterXic(truncatedXics[charge - MinCharge][charge - MinCharge],
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
            const int smoothingWindowOneSideSize = 1;

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
                        if (r < x.Length)
                        {
                            sum += x[r];
                        }

                        if (k >= 0 && k < x.Length)
                        {
                            s[k] = sum;
                            //if(x[k]>0)Console.WriteLine(s[k] + "\t" + _xicArray[i][j][k]);

                            //maxS = maxS < s[k] ? s[k] : maxS;
                        }
                        if (l >= 0)
                        {
                            sum -= x[l];
                        }
                    }
                    /*if (maxS > 0)
                    {
                        for (var k = 0; k < s.Length; k++)
                            s[k] = s[k]/maxS;
                    }*/
                    smoothedXicArray[i][j] = s;
                }
            }
            // System.Environment.Exit(1);
            //  for (var i = 0; i < smoothedXicArray[15 - MinCharge][_maxIntensityIsotopeIndex-_minIsotopeIndex].Length; i++ )
            // {
            //    Console.WriteLine(i + "\t" + smoothedXicArray[15 - MinCharge][_maxIntensityIsotopeIndex - _minIsotopeIndex][i] + "\t" + _xicArray[15 - MinCharge][_maxIntensityIsotopeIndex - _minIsotopeIndex][i]);
            // }
            return smoothedXicArray;
        }

        private double[][][] GetXicArray()
        {
            var scanNumberIndex = new Dictionary<int, int>();
            var i = 0;
            foreach (var scanNumber in _run.GetScanNumbers(1))
            {
                // Console.WriteLine(scanNumber + " " + i);
                scanNumberIndex[scanNumber] = i++;
            }

            var xicArray = new double[MaxCharge - MinCharge + 1][][];
            // Console.WriteLine("***");
            for (var charge = MinCharge; charge <= MaxCharge; charge++)
            {
                var xicArrayForThisCharge = new double[_isotopeEnvelope.Length][];
                for (var j = 0; j < xicArrayForThisCharge.Length; j++)
                {
                    xicArrayForThisCharge[j] = new double[scanNumberIndex.Count];
                }

                var precursorIon = new Ion(_proteinCompositionPlusWater, charge);

                // Console.Write("precursorCharge " + charge);

                for (var k = 0; k < xicArrayForThisCharge.Length; k++)
                {
                    var mz = precursorIon.GetIsotopeMz(k + _minIsotopeIndex);
                    var xic = _run.GetPrecursorExtractedIonChromatogram(mz, _tolerance);

                    //Console.WriteLine("Mz" + mz);

                    foreach (var xicPeak in xic)
                    {
                        // Console.WriteLine(xicPeak.MostAbundantIsotopeMz);
                        xicArrayForThisCharge[k][scanNumberIndex[xicPeak.ScanNum]] = xicPeak.Intensity;
                        // if (k == _maxIntensityIsotopeIndex-_minIsotopeIndex)
                        //       Console.WriteLine(charge + "\t" + xicPeak.MostAbundantIsotopeMz + "\t" + xicPeak.Intensity + "\t" + scanNumberIndex[xicPeak.MostAbundantIsotopeMz]);
                        //if(k == 0 && charge == 13) Console.WriteLine(xicPeak.MostAbundantIsotopeMz + " " + xicPeak.Intensity);
                    }
                    //  Console.WriteLine();
                }
                // System.Environment.Exit(1);
                //
                xicArray[charge - MinCharge] = xicArrayForThisCharge;
            }
            return xicArray;
        }
    }
}
