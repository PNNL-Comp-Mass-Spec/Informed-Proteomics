using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.TopDown.Quantification;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class UnidentifiedFeatureAnalysisTest
    {
        [Test]
        public void TestUnidentifiedFeatureAnalysis()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            var crossTabFile = @"C:\Users\mend645\Desktop\TopDown_Analyis\CPTAC_Intact_CompRef_Analysis.txt";
            var outputFolder = @"C:\Users\mend645\Desktop\TopDown_Analyis\";
            var rawFiles = new[]
            {
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR32A_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR32B_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR32C_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR32D_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR32E_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR32F_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR32G_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR33A_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR33B_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR33C_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR33D_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR33E_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR33F_24Aug15_Bane_15-02-06-RZ.pbf",
                @"\\protoapps\UserData\Jungkap\CompRef\raw\CPTAC_Intact_CR33G_24Aug15_Bane_15-02-06-RZ.pbf"
            };
            var databaseFile = @"\\protoapps\UserData\Jungkap\CompRef\db\H_sapiens_M_musculus_Trypsin_NCBI_Build37_2011-12-02.fasta";


            var unFeatureAnalyzer = new UnidentifiedFeatureAnalysis(rawFiles,crossTabFile,databaseFile);
            Console.WriteLine("Filtering Features.............");
            unFeatureAnalyzer.FilterFeatures(outputFolder);
            unFeatureAnalyzer.DoAnalysis();
            unFeatureAnalyzer.PrintAnalysisToFile(outputFolder);
        }
    }
}
