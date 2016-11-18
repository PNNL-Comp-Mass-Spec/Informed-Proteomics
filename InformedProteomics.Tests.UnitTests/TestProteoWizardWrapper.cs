using System.Reflection;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Tests.Base;
using NUnit.Framework;

namespace InformedProteomics.Tests.UnitTests
{
    [TestFixture]
    internal class TestProteoWizardWrapper
    {
        public const string TestRawFilePath = @"\\proto-2\UnitTest_Files\InformedProteomics_TestFiles\SpecFiles\QC_Shew_12_02_2_1Aug12_Cougar_12-06-11.raw";

        [Test]
        public void TestLoadingProteoWizardWrapper()
        {
            var methodName = MethodBase.GetCurrentMethod().Name;
            Utils.ShowStarting(methodName);

            //try
            //{
            var reader = new ProteoWizardReader(TestRawFilePath);
            reader.ReadMassSpectrum(1);

            System.Console.WriteLine("NumSpectra: {0}", reader.NumSpectra);

            //}
            //catch (FileNotFoundException e)
            //{
            //    Console.WriteLine("FileNotFound: {0}", e.FileName);
            //}
            //var mzs = new[] {1.0};
            //var intensities = new[] {10000.0};
            //var centroider = new Centroider(mzs, intensities);
            //centroider.GetCentroidedData(out mzs, out intensities);
        }
    }
}
