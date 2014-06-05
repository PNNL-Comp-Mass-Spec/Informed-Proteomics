using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.TopDownViewer;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestIcFileReader
    {
        [Test]
        public void IcFileReader()
        {
            const string rawFile = @"\\protoapps\UserData\Wilkins\TopDown\Vlad\raw\Alz_RA_C1_HCD_11012013_SW_03Nov2013.raw";
            const string tsvFile = @"\\protoapps\UserData\Wilkins\TopDown\Vlad\tsv\Alz_RA_C1_HCD_11012013_SW_03Nov2013_IcTda.tsv";
            var reader = new IcFileReader(tsvFile, rawFile);
            var prsms = reader.Read();
//            var charge = prsms[0].GetCharge(2);
        }
    }
}
