using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.UnitTests
{
    [TestFixture]
    internal class TestMzMLReader
    {
        [Test]
        [TestCase(@"TEST_FOLDER\QC_Shew_Intact_26Sep14_Bane_C2Column3_Excerpt.mzML", 1001)]
        public void TestReadMzMLStd(string filePath, int expectedSpectra)
        {
            var methodName = MethodBase.GetCurrentMethod()?.Name;
            var mzMLFile = Utils.GetTestFile(methodName, filePath.Replace("TEST_FOLDER", Utils.DEFAULT_SPEC_FILES_FOLDER));

            TestReadMzML(mzMLFile, expectedSpectra);
        }

        [Test]
        [TestCase(@"TEST_FOLDER\MZML\Online_Dig_v17_QC_Shew_Stab-02_c0-5_01_04Aug14_Alder_14-06-11.mzML", 17288)] // Centroid, Thermo/XCalibur
        [TestCase(@"TEST_FOLDER\MZML\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.mzML", 31163)] // Profile, Thermo/XCalibur
        [TestCase(@"TEST_FOLDER\MZML\QC_Shew_13_06_500ng_CID_2_9Aug14_Lynx_14-04-08.mzML", 5649)] // Centroid, Thermo/XCalibur
        [TestCase(@"TEST_FOLDER\MZML\CTRL_Dam_17022011_1.mzML", 18696)] // Centroid, Waters/MassLynx, mzML 1.0.0, requires the referenceable Param Group
        [TestCase(@"TEST_FOLDER\MZML\napedro_L120224_005_SW-400AQUA no background 2ul dilution 6.mzML", 78012)] // Centroid, ABI SCIex WIFF files
        [TestCase(@"TEST_FOLDER\MZML\VA139IMSMS.mzML", 3145)] // Centroid, Agilent QTOF
        [TestCase(@"TEST_FOLDER\MZML\VA139IMSMS.mzML.gz", 3145)] // Centroid, Agilent QTOF, gzipped
        [TestCase(@"TEST_FOLDER\MZML\VA139IMSMS_compressed.mzML", 3145)] // Centroid, Agilent QTOF, compressed binary data
        [Category("PNL_Domain")]
        public void TestReadMzMLAddnl(string filePath, int expectedSpectra)
        {
            var methodName = MethodBase.GetCurrentMethod()?.Name;
            var mzMLFile = Utils.GetTestFile(methodName, filePath.Replace("TEST_FOLDER", Utils.DEFAULT_TEST_FILE_FOLDER));

            TestReadMzML(mzMLFile, expectedSpectra);
        }

        private void TestReadMzML(FileSystemInfo mzMLFile, int expectedSpectra)
        {
            var methodName = MethodBase.GetCurrentMethod()?.Name;
            Utils.ShowStarting(methodName, mzMLFile.Name);

            var timer = new Stopwatch();
            timer.Start();

            var reader = new MzMLReader(mzMLFile.FullName);
            var constTime = timer.Elapsed;
            Console.WriteLine("Constructor time: " + constTime);

            var numSpectra = reader.NumSpectra;
            var metaTime = timer.Elapsed - constTime;
            Console.WriteLine("Metadata read time: " + metaTime);

            var spectra = reader.ReadAllSpectra();
            var spectraCount = spectra.Count();
            timer.Stop();

            reader.Close();

            Console.WriteLine("Spectra Read time: " + (timer.Elapsed - metaTime));
            Console.WriteLine("Time: " + timer.Elapsed);

            Assert.AreEqual(expectedSpectra, numSpectra, "NumSpectra");
            Assert.AreEqual(expectedSpectra, spectraCount, "SpectraCount");
        }
    }
}
