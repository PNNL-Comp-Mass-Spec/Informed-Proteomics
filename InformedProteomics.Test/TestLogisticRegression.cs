using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestLogisticRegression
    {
        [Test]
        public void LogisticRegressionTest()
        {

            var obsverationFile = @"\\protoapps\UserData\Jungkap\Joshua\score_training\observations.txt";
            var categoriesFile = @"\\protoapps\UserData\Jungkap\Joshua\score_training\category_label.txt";

            var obs = new List<double[]>();
            var cat = new List<int>();
            var line = "";
            var logistic = new LogisticRegression();

            System.IO.StreamReader file = new System.IO.StreamReader(obsverationFile);
            while ((line = file.ReadLine()) != null)
            {
                var lineElements = line.Split('\t');
                var doubleElements = new double[lineElements.Length];
                for (var i = 0; i < lineElements.Length; i++)
                {
                    doubleElements[i] = Convert.ToDouble(lineElements[i]);
                }
                obs.Add(doubleElements);
            }

            file = new System.IO.StreamReader(categoriesFile);
            while ((line = file.ReadLine()) != null)
            {
                cat.Add(Convert.ToInt32(line));
            }

            var obsMatrix = obs.ToArray();
            var classVector = cat.ToArray();

            var weights = new[]
            {
                -3.94700432551079,
                -0.304056338862474,
                -0.908148491499137,
                0.232858075077728,
                0.340450694800618,
                0.949161406416399,
                0.268541376564235,
                -0.0198709768288523,
                -0.842410987704276,
                0.196949733855582,
                0.398403018300038,
                0.860192152276247,
                -0.0287180268260492,
                -0.0178256644915106
            };

            var stopwatch = Stopwatch.StartNew();
            weights = logistic.ComputeWeights(obsMatrix, classVector,weights);
            stopwatch.Stop();

            Console.WriteLine();
            Console.WriteLine("Gradient Descent Time: {0}", stopwatch.ElapsedMilliseconds / 1000.0);
            foreach (var x in weights)
            {
                Console.WriteLine(x);
            }
            Console.WriteLine();

            var oneCount = 0;
            var zeroCount = 0;
            var totalRight = 0;
            for (var i = 0; i < 10000; i++)
            {
                var score = CalculateClassificationScore(obsMatrix[i], weights);
                var binary = GetClassification(score);
                if (binary == cat[i]) totalRight++;
                if (binary == 1) oneCount++;
                if (binary == 0) zeroCount++;
            }

            Console.WriteLine("Scanning for 0");
            Console.WriteLine("The number of 1: {0}",oneCount);
            Console.WriteLine("The number of 0: {0}",zeroCount);

            oneCount = 0;
            zeroCount = 0;

            for (var i = 10000; i < 20000; i++)
            {
                var score = CalculateClassificationScore(obsMatrix[i], weights);
                var binary = GetClassification(score);
                if (binary == cat[i]) totalRight++;
                if (binary == 1) oneCount++;
                if (binary == 0) zeroCount++;
            }

            Console.WriteLine("Scanning for 1");
            Console.WriteLine("The number of 1: {0}", oneCount);
            Console.WriteLine("The number of 0: {0}", zeroCount);

            Console.WriteLine("Total Right: {0}",totalRight);
        }

        private int GetClassification(double x)
        {
            return (x > 0) ? 1 : 0;
        }

        private double CalculateClassificationScore(double[] data, double[] weights)
        {
            var sum = 0.0;
            sum += weights[0];
            for (var i = 0; i < data.Length; i++)
            {
                sum += data[i] * weights[i+1];
            }
            return sum;
        }
    }
}
