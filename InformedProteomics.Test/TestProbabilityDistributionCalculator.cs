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
    public class TestProbabilityDistributionCalculator
    {

        [Test]
        public void TestProbabilityDistribution()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);
            var probabilityCalc = new ProbabilityDistributionCalculator();
            //var probabilityCalcTarget = new ProbabilityDistributionCalculator();
            var observations = new List<double[]>();
            var testDecoyObsFile = @"\\protoapps\UserData\Jungkap\Joshua\LcMsFeatureScoreTrain\decoy_trainset.tsv";
            var testTargetObsFile = @"\\protoapps\UserData\Jungkap\Joshua\LcMsFeatureScoreTrain\target_trainset.tsv";


            string line; 
            var file = new System.IO.StreamReader(testDecoyObsFile);
            while ((line = file.ReadLine()) != null)
            {
                var values = line.Split('\t');
                var valArray = new double[values.Length];
                for (var i = 0; i < valArray.Length; i++)
                {
                    valArray[i] = double.Parse(values[i]);
                }
                observations.Add(valArray);
            }

            file.Close();

            probabilityCalc.AddObservations(observations);
            probabilityCalc.GetDistributions(@"C:\Users\mend645\Desktop\Decoy_Prob\");
            probabilityCalc.ClearObservations();

            file = new System.IO.StreamReader(testTargetObsFile);
            observations.Clear();

            while ((line = file.ReadLine()) != null)
            {
                var values = line.Split('\t');
                var valArray = new double[values.Length];
                for (var i = 0; i < valArray.Length; i++)
                {
                    valArray[i] = double.Parse(values[i]);
                }
                observations.Add(valArray);
            }

            file.Close();

            probabilityCalc.AddObservations(observations);
            probabilityCalc.GetDistributions(@"C:\Users\mend645\Desktop\Target_Prob\");
            //probabilityCalcTarget.AddObservations(observations);
            //probabilityCalcTarget.GetDistributions(@"C:\Users\mend645\Desktop\Target_Prob\");

            var variableNames = new[] { "D1", "C1", "I1", "D2", "C2", "I2", "R", "X1", "X2" };

            var decoyFileLocation = @"C:\Users\mend645\Desktop\Decoy_Prob\";
            var targetFileLocation = @"C:\Users\mend645\Desktop\Target_Prob\";

            for (var i = 0; i < variableNames.Length; i++)
            {
                var decoyFile = decoyFileLocation + variableNames[i]+ ".tsv";
                var targetFile = targetFileLocation + variableNames[i] + ".tsv";
                var outPutFile = @"C:\Users\mend645\Desktop\Likelihood\" + variableNames[i] + "_likelihood.tsv";
                probabilityCalc.GetLikelihoods(targetFile, decoyFile, outPutFile);
            }

            


        }
    }
}
