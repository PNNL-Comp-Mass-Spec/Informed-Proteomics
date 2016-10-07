using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.Backend.Utils
{
    public class ProbabilityDistributionCalculator
    {
        public ProbabilityDistributionCalculator()
        {
            _massBins = Enumerable.Range(0, 31).Select(x => 500 + x * 1000).ToArray();
            _massBins[_massBins.Length - 1] = 30000;
            _densityEstimators = new KernelDensityEstimator[_massBins.Length];

            for (var i = 0; i < _densityEstimators.Length; i++)
            {
                _densityEstimators[i] = new KernelDensityEstimator();
            }
        }

        public void GetDistributions(string folderLocation)
        {
            var samplingPoints = Enumerable.Range(0, 1001).Select(x => x*.001).ToArray();
            var variableNames = new [] {"D1", "C1", "I1", "D2", "C2", "I2", "R", "X1", "X2"};
            const double bandwith = .02;

            for (var i = 0; i < variableNames.Length; i++)
            {
                var tsv = new StringBuilder();
                for (var j = 0; j < _densityEstimators.Length; j++)
                {
                    var results = _densityEstimators[j].GetDensityEstimation(samplingPoints, bandwith, i + 1);
                    var sumResults = results.Sum();
                    for (var k = 0; k < results.Length - 1; k++)
                    {
                        tsv.Append(string.Format("{0} \t", results[k] / sumResults));
                    }
                    tsv.Append(string.Format("{0}", results[results.Length-1] / sumResults));
                    tsv.Append("\n");
                }
                File.WriteAllText(folderLocation + variableNames[i] + ".tsv", tsv.ToString());
            }
        }

        public void GetLikelihoods(string targetFileLocation, string decoyFileLocation, string outputLocation)
        {
            string targetLine;
            string decoyLine;

            var targetFile = new StreamReader(targetFileLocation);
            var decoyFile = new StreamReader(decoyFileLocation);

            var tsv = new StringBuilder();

            while (((targetLine = targetFile.ReadLine()) != null) && ((decoyLine = decoyFile.ReadLine()) != null))
            {
                var targetVals = targetLine.Split('\t');
                var decoyVals = decoyLine.Split('\t');

                double targetVal;
                double decoyVal;
                for (var i = 0; i < targetVals.Length - 1 ; i++)
                {
                    targetVal = double.Parse(targetVals[i]);
                    decoyVal = double.Parse(decoyVals[i]);
                    tsv.Append(string.Format("{0} \t", Math.Log(targetVal/decoyVal)));
                }

                targetVal = double.Parse(targetVals[targetVals.Length-1]);
                decoyVal = double.Parse(decoyVals[decoyVals.Length-1]);
                tsv.Append(string.Format("{0} \t", Math.Log(targetVal / decoyVal)));
                tsv.Append("\n");
            }
            File.WriteAllText(outputLocation, tsv.ToString());
        }

        public void AddObservations(IEnumerable<double[]> observations)
        {
            foreach (var obs in observations)
            {
                for (var i = 0; i < _densityEstimators.Length; i++)
                {
                    var minMass = _massBins[i]*.7;
                    var maxMass = _massBins[i]*1.3;
                    var obsMass = obs[0];
                    if (i == 0)
                    {
                        minMass = 0;
                        maxMass = 1000;
                    }
                    if (i == _massBins.Length - 1)
                    {
                        minMass = 23000;
                        maxMass = int.MaxValue;
                    }

                    if(obsMass >= minMass && obsMass <= maxMass) _densityEstimators[i].AddObservation(obs);
                }
            }
        }

        public void ClearObservations()
        {
            foreach (var kernel in _densityEstimators)
            {
                kernel.ClearObservations();
            }
        }

        private readonly int[] _massBins;
        private readonly KernelDensityEstimator[] _densityEstimators;
    }
}
