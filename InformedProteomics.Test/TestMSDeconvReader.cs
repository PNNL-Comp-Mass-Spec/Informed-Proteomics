using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.TopDown.Quantification;

namespace InformedProteomics.Test
{
    [TestFixture]
    class TestMSDeconvReader
    {
        [Test]
        public void MsDeconvReaderTest()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName);

            /*var txtFiles = new[]
            {
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep10_15Jan15_Bane_C2-14-08-02RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep8_15Jan15_Bane_C2-14-08-02RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep7_15Jan15_Bane_C2-14-08-02RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep5_15Jan15_Bane_C2-14-08-02RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep4_15Jan15_Bane_C2-14-08-02RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep3_15Jan15_Bane_C2-14-08-02RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep1_15Jan15_Bane_C2-14-08-02RZ_msdeconv.txt",
            };*/

            /**var csvFiles = new[]
            {
                @"\\Proto-5\VOrbiETD02\2015_1\CPTAC_Intact_rep10_15Jan15_Bane_C2-14-08-02RZ\DLS201501211509_Auto1153791\CPTAC_Intact_rep10_15Jan15_Bane_C2-14-08-02RZ_isos.csv",
                @"\\Proto-5\VOrbiETD02\2015_1\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ\DLS201501211509_Auto1153781\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ_isos.csv",
                @"\\Proto-5\VOrbiETD02\2015_1\CPTAC_Intact_rep8_15Jan15_Bane_C2-14-08-02RZ\DLS201501211509_Auto1153789\CPTAC_Intact_rep8_15Jan15_Bane_C2-14-08-02RZ_isos.csv",
                @"\\Proto-5\VOrbiETD02\2015_1\CPTAC_Intact_rep7_15Jan15_Bane_C2-14-08-02RZ\DLS201501211509_Auto1153783\CPTAC_Intact_rep7_15Jan15_Bane_C2-14-08-02RZ_isos.csv",
                @"\\Proto-5\VOrbiETD02\2015_1\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ\DLS201501211509_Auto1153799\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ_isos.csv",
                @"\\Proto-5\VOrbiETD02\2015_1\CPTAC_Intact_rep5_15Jan15_Bane_C2-14-08-02RZ\DLS201501211509_Auto1153801\CPTAC_Intact_rep5_15Jan15_Bane_C2-14-08-02RZ_isos.csv",
                @"\\Proto-5\VOrbiETD02\2015_1\CPTAC_Intact_rep4_15Jan15_Bane_C2-14-08-02RZ\DLS201501211509_Auto1153795\CPTAC_Intact_rep4_15Jan15_Bane_C2-14-08-02RZ_isos.csv",
                @"\\Proto-5\VOrbiETD02\2015_1\CPTAC_Intact_rep3_15Jan15_Bane_C2-14-08-02RZ\DLS201501211509_Auto1153793\CPTAC_Intact_rep3_15Jan15_Bane_C2-14-08-02RZ_isos.csv",
                @"\\Proto-5\VOrbiETD02\2015_1\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ\DLS201501211509_Auto1153797\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ_isos.csv",
                @"\\Proto-5\VOrbiETD02\2015_1\CPTAC_Intact_rep1_15Jan15_Bane_C2-14-08-02RZ\DLS201501211509_Auto1153775\CPTAC_Intact_rep1_15Jan15_Bane_C2-14-08-02RZ_isos.csv"
            };**/

            /*var rawFiles = new[]
            {
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep10_15Jan15_Bane_C2-14-08-02RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep9_15Jan15_Bane_C2-14-08-02RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep8_15Jan15_Bane_C2-14-08-02RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep7_15Jan15_Bane_C2-14-08-02RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep6_15Jan15_Bane_C2-14-08-02RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep5_15Jan15_Bane_C2-14-08-02RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep4_15Jan15_Bane_C2-14-08-02RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep3_15Jan15_Bane_C2-14-08-02RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep2_15Jan15_Bane_C2-14-08-02RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files\CPTAC_Intact_rep1_15Jan15_Bane_C2-14-08-02RZ.raw"
            };*/

            var txtFiles = new[]
            {
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_1x_1_27Apr15_Bane_14-09-03RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_1x_2_27Apr15_Bane_14-09-03RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_1x_3_27Apr15_Bane_14-09-03RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_1x_4_27Apr15_Bane_14-09-03RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_1x_5_27Apr15_Bane_14-09-03RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_5x_1_27Apr15_Bane_14-09-03RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_5x_2_27Apr15_Bane_14-09-03RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_5x_3b_30Apr15_Bane_14-09-01RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_5x_4_27Apr15_Bane_14-09-01RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_5x_5_27Apr15_Bane_14-09-01RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_10x_1_27Apr15_Bane_14-09-01RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_10x_2_27Apr15_Bane_14-09-01RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_10x_3_2May15_Bane_14-09-01RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_10x_4_27Apr15_Bane_14-09-01RZ_msdeconv.txt",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_10x_5_27Apr15_Bane_14-09-01RZ_msdeconv.txt"
            };

        var csvFiles = new[]
            {
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_1x_1_27Apr15_Bane_14-09-03RZ\DLS201504290957_Auto1191334\CPTAC_Intact_Spike_1x_1_27Apr15_Bane_14-09-03RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_1x_2_27Apr15_Bane_14-09-03RZ\DLS201504290958_Auto1191337\CPTAC_Intact_Spike_1x_2_27Apr15_Bane_14-09-03RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_1x_3_27Apr15_Bane_14-09-03RZ\DLS201504290956_Auto1191331\CPTAC_Intact_Spike_1x_3_27Apr15_Bane_14-09-03RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_1x_4_27Apr15_Bane_14-09-03RZ\DLS201504290958_Auto1191340\CPTAC_Intact_Spike_1x_4_27Apr15_Bane_14-09-03RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_1x_5_27Apr15_Bane_14-09-03RZ\DLS201505011038_Auto1192336\CPTAC_Intact_Spike_1x_5_27Apr15_Bane_14-09-03RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_5x_1_27Apr15_Bane_14-09-03RZ\DLS201505011041_Auto1192348\CPTAC_Intact_Spike_5x_1_27Apr15_Bane_14-09-03RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_5x_2_27Apr15_Bane_14-09-03RZ\DLS201505011040_Auto1192345\CPTAC_Intact_Spike_5x_2_27Apr15_Bane_14-09-03RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_5x_3b_30Apr15_Bane_14-09-01RZ\DLS201505011039_Auto1192339\CPTAC_Intact_Spike_5x_3b_30Apr15_Bane_14-09-01RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_5x_4_27Apr15_Bane_14-09-01RZ\DLS201505061337_Auto1193464\CPTAC_Intact_Spike_5x_4_27Apr15_Bane_14-09-01RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_5x_5_27Apr15_Bane_14-09-01RZ\DLS201505061337_Auto1193461\CPTAC_Intact_Spike_5x_5_27Apr15_Bane_14-09-01RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_10x_1_27Apr15_Bane_14-09-01RZ\DLS201505061337_Auto1193455\CPTAC_Intact_Spike_10x_1_27Apr15_Bane_14-09-01RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_10x_2_27Apr15_Bane_14-09-01RZ\DLS201505061337_Auto1193458\CPTAC_Intact_Spike_10x_2_27Apr15_Bane_14-09-01RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_10x_3_2May15_Bane_14-09-01RZ\DLS201505061337_Auto1193449\CPTAC_Intact_Spike_10x_3_2May15_Bane_14-09-01RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_10x_4_27Apr15_Bane_14-09-01RZ\DLS201505061337_Auto1193467\CPTAC_Intact_Spike_10x_4_27Apr15_Bane_14-09-01RZ_isos.csv",
                @"\\proto-4\VOrbiETD02\2015_2\CPTAC_Intact_Spike_10x_5_27Apr15_Bane_14-09-01RZ\DLS201505061337_Auto1193446\CPTAC_Intact_Spike_10x_5_27Apr15_Bane_14-09-01RZ_isos.csv"
            };

            var rawFiles = new[]
            {
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_1x_1_27Apr15_Bane_14-09-03RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_1x_2_27Apr15_Bane_14-09-03RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_1x_3_27Apr15_Bane_14-09-03RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_1x_4_27Apr15_Bane_14-09-03RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_1x_5_27Apr15_Bane_14-09-03RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_5x_1_27Apr15_Bane_14-09-03RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_5x_2_27Apr15_Bane_14-09-03RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_5x_3b_30Apr15_Bane_14-09-01RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_5x_4_27Apr15_Bane_14-09-01RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_5x_5_27Apr15_Bane_14-09-01RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_10x_1_27Apr15_Bane_14-09-01RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_10x_2_27Apr15_Bane_14-09-01RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_10x_3_2May15_Bane_14-09-01RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_10x_4_27Apr15_Bane_14-09-01RZ.raw",
                @"C:\Users\mend645\Documents\Raw_Files2\CPTAC_Intact_Spike_10x_5_27Apr15_Bane_14-09-01RZ.raw"
            };

            var filePath = @"C:\Users\mend645\Documents\Cluster_Files2\";

            var msDeconvReader = new MSDeconvReader();
            msDeconvReader.MsLevel = 1;
            msDeconvReader.MinCharge = 2;
            msDeconvReader.MinMass = 2000.0;

            for (int i = 0; i < rawFiles.Length; i++)
            {
                var nodeList = msDeconvReader.GetDeconvNodesForMsDeconv(txtFiles[i]);
                //var nodeList = msDeconvReader.GetDeconvNodesForDecon2Ls(csvFiles[i]);
                var run = PbfLcMsRun.GetLcMsRun(rawFiles[i]);
                var tol = new Tolerance(10);
                var elutionInterval = .5;
                var stringPieces = txtFiles[i].Split(new char[] { '\\', '.' });
                //var stringPieces = csvFiles[i].Split(new char[] {'\\', '.'});
                var fileName = filePath + stringPieces[stringPieces.Length - 2];
                var clusterer = new MSDeconvClusterer(run);
                Console.WriteLine("start file {0}", i);
                var connectedComponents = clusterer.GetClustersList(tol, elutionInterval, nodeList);
                Console.WriteLine("end file {0}", i);
                clusterer.SaveClusterInfoToFile(fileName,connectedComponents);
            }
        }
    }
}
