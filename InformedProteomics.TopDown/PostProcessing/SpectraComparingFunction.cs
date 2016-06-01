using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Scoring.TopDown;

namespace InformedProteomics.TopDown.PostProcessing
{
    public class SpectraComparingFunction
    {

        public double DotProduct(Peak[] spectrum1, Peak[] spectrum2, FilteredProteinMassBinning comparer)
        {

            spectrum1 = GuassianFilter(spectrum1,.5);
            spectrum2 = GuassianFilter(spectrum2,.5);
            var vectorLength = comparer.GetBinNumber(10000.0);
            var featureVector1 = ConvertToFullIntensityVector(spectrum1,vectorLength,comparer);
            var featureVector2 = ConvertToFullIntensityVector(spectrum2, vectorLength,comparer);


            var sum = 0d;
            for (int i = 0; i < vectorLength; i++)
            {
                sum += featureVector1[i]*featureVector2[i];
            }

            var norm1 = featureVector1.Sum(x => x*x);
            var norm2 = featureVector2.Sum(x => x * x);
            return sum/Math.Sqrt(norm1*norm2);
        }

        public double PearsonCorrelation(Peak[] spectrum1, Peak[] spectrum2,FilteredProteinMassBinning comparer)
        {

            var spec1Bar = 0d;
            var spec2Bar = 0d;

            spectrum1 = GuassianFilter(spectrum1, .5);
            spectrum2 = GuassianFilter(spectrum2, .5);
            var vectorLength = comparer.GetBinNumber(10000.0);

            spec1Bar = spectrum1.Sum(x => x.Intensity)/vectorLength;
            spec2Bar = spectrum1.Sum(y => y.Intensity)/vectorLength;


            var intensityVector1 = ConvertToFullIntensityVector(spectrum1, vectorLength, comparer);
            var intensityVector2 = ConvertToFullIntensityVector(spectrum2, vectorLength, comparer);

            var cov = 0.0;
            var s1 = 0.0;
            var s2 = 0.0;

            for (var i = 0; i < vectorLength; i++)
            {
                var d1 = intensityVector1[i] - spec1Bar;
                var d2 = intensityVector2[i] - spec2Bar;
                cov += d1*d2;
                s1 += d1*d1;
                s2 += d2*d2;
            }

            if (s1 <= 0 || s2 <= 0) return 0;
            return cov < 0 ? 0f : cov/Math.Sqrt(s1*s2);
        }

        public double RootMeanSquareDeviation(Peak[] spectrum1, Peak[] spectrum2, FilteredProteinMassBinning comparer)
        {

            spectrum1 = GuassianFilter(spectrum1, .5);
            spectrum2 = GuassianFilter(spectrum2, .5);
            var vectorLength = comparer.GetBinNumber(10000.0);
            var intensityVector1 = ConvertToFullIntensityVector(spectrum1, vectorLength, comparer);
            var intensityVector2 = ConvertToFullIntensityVector(spectrum2, vectorLength, comparer);
            var mean1 = spectrum1.Sum(p => p.Intensity)/spectrum1.Length;
            var mean2 = spectrum1.Sum(p => p.Intensity)/spectrum2.Length;

            var sum = 0d;
            for (int i = 0; i < vectorLength; i++)
            {
                var diff = intensityVector1[i] - intensityVector2[i];
                sum += diff*diff;
            }

            return Math.Sqrt(sum/vectorLength);
        }


        public double[] ConvertToFullIntensityVector(Peak[] spectrum, int length, FilteredProteinMassBinning comparer)
        {
            var intensityVector = new double[length];
            Array.Clear(intensityVector,0,length);
            spectrum = RemovePeaks(10000.0, spectrum);

            for (var i = 0; i < spectrum.Length; i++)
            {
                var binNumber = comparer.GetBinNumber(spectrum[i].Mz);
                if (binNumber >= 0)
                {
                    intensityVector[binNumber - 1] = spectrum[i].Intensity;
                }
            }
            return intensityVector;
        }

        public double GetRegressionSlope(Peak[] spectrum)
        {
            var n = spectrum.Length;
            var sumXY = spectrum.Sum(p => p.Mz*p.Intensity);
            var sumX = spectrum.Sum(p => p.Mz);
            var sumY = spectrum.Sum(p => p.Intensity);
            var sumXSquared = spectrum.Sum(p => p.Mz*p.Mz);

            return (n*sumXY - sumX*sumY)/(n*sumXSquared - sumX*sumX);
        }

        public double GetRegressionIntercept(Peak[] spectrum, double b)
        {
            var n = spectrum.Length;
            var sumX = spectrum.Sum(p => p.Mz);
            var sumY = spectrum.Sum(p => p.Intensity);
            return (sumY - b*sumX)/n;
        }

        public double StandardDeviation(double[] spectrum, double mean)
        {
            var sum = 0d;
            for (int i = 0; i < spectrum.Length; i++)
            {
                var diff = spectrum[i] - mean;
                sum += diff*diff;
            }
            return Math.Sqrt(sum/spectrum.Length);
        }

        public Peak[] GuassianFilter(Peak[] spectrum, double std = 1)
        {
            var n = spectrum.Length;
            var newPeak = new Peak[n];
            for (var i = 0; i < n; i++)
            {
                var sum = 0d;
                for (var j = 0; j < n; j++)
                {
                    var intensity = spectrum[j].Intensity;
                    var massDiff = spectrum[j].Mz - spectrum[i].Mz;
                    sum += intensity*Math.Exp(-(massDiff*massDiff)/(2*std*std));
                }
                newPeak[i] = new Peak(spectrum[i].Mz,sum/(n*2*std*std));
            }
            return newPeak;
        }


        public Peak[] RemovePeaks(double mass, Peak[] peakList)
        {
            return peakList.Where(p => p.Mz <= mass).ToArray();
        }

    }
}
