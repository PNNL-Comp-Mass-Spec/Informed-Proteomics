using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    internal class TestMzMLReader
    {
        [Test]
        [TestCase(@"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MZML\Online_Dig_v17_QC_Shew_Stab-02_c0-5_01_04Aug14_Alder_14-06-11.mzML", 17288)] // Centroid, Thermo/XCaliber
        [TestCase(@"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MZML\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.mzML", 31163)] // Profile, Thermo/XCaliber
        [TestCase(@"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MZML\QC_Shew_13_06_500ng_CID_2_9Aug14_Lynx_14-04-08.mzML", 5649)] // Centroid, Thermo/XCaliber
        [TestCase(@"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MZML\CTRL_Dam_17022011_1.mzML", 18696)] // Centroid, Waters/MassLynx, mzML 1.0.0, requires the referenceable Param Group
        [TestCase(@"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MZML\napedro_L120224_005_SW-400AQUA no background 2ul dilution 6.mzML", 78012)] // Centroid, ABI SCIex WIFF files
        [TestCase(@"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MZML\VA139IMSMS.mzML", 3145)] // Centroid, Agilent QTOF
        [TestCase(@"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MZML\VA139IMSMS.mzML.gz", 3145)] // Centroid, Agilent QTOF, gzipped
        [TestCase(@"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\MZML\VA139IMSMS_compressed.mzML", 3145)] // Centroid, Agilent QTOF, compressed binary data
        public void TestReadMzML(string filePath, int expectedSpectra)
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            TestUtils.ShowStarting(methodName + " (" + Path.GetFileName(filePath) + ")");

            Stopwatch timer = new Stopwatch();
            timer.Start();
            var reader = new MzMLReader(filePath);
            var constTime = timer.Elapsed;
            Console.WriteLine(@"Constructor time: " + constTime);
            var numSpectra = reader.NumSpectra;
            var metaTime = timer.Elapsed - constTime;
            Console.WriteLine(@"Metadata read time: " + metaTime);
            var spectra = reader.ReadAllSpectra();
            var spectraCount = spectra.Count();
            timer.Stop();
            Console.WriteLine(@"Spectra Read time: " + (timer.Elapsed - metaTime));
            Console.WriteLine(@"Time: " + timer.Elapsed);
            Assert.AreEqual(expectedSpectra, numSpectra, "NumSpectra");
            Assert.AreEqual(expectedSpectra, spectraCount, "SpectraCount");
            Assert.AreEqual(numSpectra, spectraCount, "NumSpectra vs. SpectraCount");
            reader.Close();
        }
    }
}
