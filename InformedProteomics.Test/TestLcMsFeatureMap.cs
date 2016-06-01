using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Graphics;
using NUnit.Framework;

namespace InformedProteomics.Test
{

    [TestFixture]
    class TestLcMsFeatureMap
    {

        [Test]
        [STAThread]
        public void TestFeatureMapGeneration()
        {
            Console.WriteLine("Testing Working");
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            const string testFile = @"\\protoapps\UserData\Jungkap\Joshua\FeatureMap\QC_Shew_Intact_26Sep14_Bane_C2Column3.ms1ft";
            const string outputFile = @"\\protoapps\UserData\Jungkap\Joshua\";

            var map = new LcMsFeatureMap(testFile,185);
            map.SaveImage(outputFile + "test.png",100);
        }

    }
}
