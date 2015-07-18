using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Test.FunctionalTests
{
    [TestFixture]
    public class TestMassSpecDataReaderFactory
    {
        [Test]
        public void TestThermoRawAvailable()
        {
            Assert.AreEqual(true, MassSpecDataReaderFactory.IsThermoRawAvailable());
        }

        [Test]
        public void TestPwizAvailable()
        {
            Assert.AreEqual(true, MassSpecDataReaderFactory.IsPwizAvailable());
        }
    }
}
