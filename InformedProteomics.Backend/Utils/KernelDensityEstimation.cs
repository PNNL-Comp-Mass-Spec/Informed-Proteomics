using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InformedProteomics.Backend.Utils
{
    public class KernelDensityEstimator
    {

        public KernelDensityEstimator()
        {
            _observations = new List<double[]>();
        }

        public double[] GetDensityEstimation(double[] vector, double bandwidth, int pointIndex)
        {
            var n = _observations.Count;
            var h = bandwidth;
            var x = new double[vector.Length];
            var c = 1/(n*h);

            if (n == 0)
            {
                Console.WriteLine("No Observation given");
                return null;
            }
            if (pointIndex > _observations[0].Length-1)
            {
                Console.WriteLine("Cannot have more points then there are elements in an observation");
                return null;
            }
           
            for (var j = 0; j < x.Length; j++)
            {
                var sum = 0.0;
                for (var i = 0; i < n; i++)
                {
                    var diff = vector[j] - _observations[i][pointIndex];
                    sum += GuassianKernel(diff/bandwidth);
                }
                x[j] = c*sum;
            }
             
            return x;
        }

        public void AddObservation(double[] obs)
        {
            if (_observations.Count != 0)
            {
                if (obs.Length != _observations[0].Length)
                {
                    Console.WriteLine("Observation length must be same for all observations will not be added");
                    return;
                }
            }
            _observations.Add(obs);
        }

        public bool HasObservation()
        {
            return (_observations.Count != 0) && true;
        }

        private double GuassianKernel(double x)
        {
            return Math.Exp(-.5*Math.Pow(x, 2))/Math.Sqrt(2*Math.PI);
        }

        public void ClearObservations()
        {
            _observations.Clear();
        }

        private readonly List<double[]> _observations;

    }
}
