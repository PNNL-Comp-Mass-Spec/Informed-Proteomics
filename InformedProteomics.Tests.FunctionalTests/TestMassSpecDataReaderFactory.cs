using InformedProteomics.Backend.MassSpecData;
using NUnit.Framework;

namespace InformedProteomics.Tests.FunctionalTests
{
    [TestFixture]
    public class TestMassSpecDataReaderFactory
    {
        [Test]
        [Category("PNL_Domain")]
        public void TestThermoRawAvailable()
        {
            Assert.AreEqual(true, MassSpecDataReaderFactory.IsThermoRawAvailable());
        }

        [Test]
        [Category("PNL_Domain")]
        public void TestPwizAvailable()
        {
            Assert.AreEqual(true, MassSpecDataReaderFactory.IsPwizAvailable());
        }
    }
}
