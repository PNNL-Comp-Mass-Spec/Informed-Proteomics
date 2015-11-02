using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestKernelDensityEstimator
    {
        [Test]
        public void TestKernel()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var x = new [] {.2, .3};
            var y = new [] { .1, .4};
            var z = new [] {.3, .4};

            var kernelEstimator = new KernelDensityEstimator();
            //kernelEstimator.AddObservation(x);
            //kernelEstimator.AddObservation(y);
            Console.WriteLine(kernelEstimator.HasObservation());
            var estimation = kernelEstimator.GetDensityEstimation(z, .1, 0);
            if (estimation == null) return;
            foreach (var e in estimation)
            {
                Console.WriteLine(e);
            }
        }
    }
}
