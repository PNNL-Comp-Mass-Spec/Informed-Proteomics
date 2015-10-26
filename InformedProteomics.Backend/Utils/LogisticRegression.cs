using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Integration;


namespace InformedProteomics.Backend.Utils
{
    public class LogisticRegression
    {

        public double[] ComputeWeights(double[][] obsMatrix, int[] classVector, double[] weights = null)
        {
            var dataMatrix = new double[obsMatrix.Length][];
            var learningRate = .1 / dataMatrix.Length;

            for (var i = 0; i < obsMatrix.Length; i++)
            {
                var obsList = obsMatrix[i].ToList();
                obsList.Add(classVector[i]);
                dataMatrix[i] = obsList.ToArray();
            }

            if (weights == null) weights = Enumerable.Repeat(0.0, obsMatrix[0].Length+1).ToArray();

            var newWeights = new double[weights.Length];
            var error = -10.0; //error is between [-1,1]
            var likelihood = double.NegativeInfinity;
            var count = 0;
            while (error < .99999999 || Double.IsNaN(error)) //!IsTerminated(dataMatrix, newWeights)
            {
                newWeights = GradientDescent(dataMatrix, weights, learningRate);
                var newError = ComputeError(weights, newWeights);
                //if(count % 100 == 0) Console.WriteLine("Count {0} Error {1}",count, newError);
                error = newError;
                weights = newWeights;
                count++;
                dataMatrix.Shuffle();
            }

            return newWeights;
        }

        private double[] GradientDescent(double[][] observed, double[] wVector, double learningRate)
        {

            var updatedWVector = new double[wVector.Length];
            wVector.CopyTo(updatedWVector,0);
            for (var i = 0; i < observed.Length; i++)
            {
                var computed = ComputeOut(observed[i], updatedWVector);
                var targetIndex = observed[i].Length - 1;
                var target = observed[i][targetIndex];
                updatedWVector[0] += learningRate*(target - computed)*1;
                for (var j = 1; j <  updatedWVector.Length; j++)
                {
                    updatedWVector[j] += learningRate*(target - computed)*observed[i][j - 1];
                }
            }
            return updatedWVector;
        }

        private double ComputeOut(double[] data, double[] weights)
        {
            var z = 0.0;
            z += weights[0];
            for (int i = 0; i < weights.Length - 1; i++)
            {
                z += (weights[i + 1]*data[i]);
            }
            return 1.0/(1 + Math.Exp(-z));
        }

        private double ComputeError(double[] w1, double[] w2)
        {
            var dotProd = w1.Zip(w2, (a, b) => a*b).Sum();
            var aEuclidLength = Math.Sqrt(w1.Zip(w1, (a, b) => a*b).Sum());
            var bEuclidLength = Math.Sqrt(w2.Zip(w2, (a, b) => a * b).Sum());
            return dotProd / (aEuclidLength * bEuclidLength);
        }

        private double LogLikelihood(double[][] obs, int[] classification, double[] weight)
        {
            var sum = 0.0;
            for (int i = 0; i < obs.Length; i++)
            {
                var y = classification[i];
                var h = ComputeOut(obs[i], weight);
                sum += y*Math.Log(h) + (1 - y)*Math.Log(1 - h);
            }
            return sum;
        }

        private bool IsTerminated(double[][] dataMatrix, double[] weight)
        {
            var correctCount = 0;
            var neededCorrect = (int)Math.Ceiling(dataMatrix.Length*.95);
            for (var i = 0; i < dataMatrix.Length; i++)
            {
                var sum = 0.0;
                sum += weight[0];
                for (var j = 0; j < dataMatrix[j].Length - 1; j++)
                {
                    sum += dataMatrix[i][j] * weight[j + 1];
                }
                var realClass = (int)dataMatrix[i][dataMatrix[0].Length - 1];
                var calcClass = GetClassification(sum);
                if (realClass == calcClass) correctCount++;
            }

            return (correctCount >= neededCorrect ? true : false);

        }

        private int GetClassification(double x)
        {
            return (x > 0) ? 1 : 0;
        }
    }
}
