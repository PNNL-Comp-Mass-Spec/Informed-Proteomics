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

            var end = Math.Min(10, prsms.Count);
            for (int i = 0; i < 10; i++)
            {
                var prsm = prsms[i];
                Console.WriteLine(prsm.Sequence);
                Console.WriteLine("\tCharges");
                var charges = prsm.Charges;
                foreach (var charge in charges)
                {
                    Console.WriteLine("\t"+charge);
                    var scans = prsm.GetScanNums(charge);
                    Console.WriteLine("\t\tScans");
                    foreach (var scan in scans) Console.WriteLine("\t\t" + scan);
                }
            }
        }

        [Test]
        public void PrSm()
        {
            const string rawFile = @"\\protoapps\UserData\Wilkins\TopDown\Vlad\raw\Alz_RA_C1_HCD_11012013_SW_03Nov2013.raw";
            const string tsvFile = @"\\protoapps\UserData\Wilkins\TopDown\Vlad\tsv\Alz_RA_C1_HCD_11012013_SW_03Nov2013_IcTda.tsv";
            var reader = new IcFileReader(tsvFile, rawFile);
            var prsms = reader.Read();

            var prsm = prsms[5];
            var charges = prsm.Charges;
            var scans = prsm.GetScanNums(charges[0]);
            var data = prsm.GetData(charges[0], scans[0]);
            Console.WriteLine(data.Sequence);
        }
    }
}
