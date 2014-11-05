using System;
using System.CodeDom;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.SequenceTag;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;
using MathNet.Numerics.Integration.Algorithms;

namespace InformedProteomics.Test
{
    [TestFixture]
    class TestLcMs1FeatureExtraction
    {
        //private const string TestRawFile = @"D:\\MassSpecFiles\\QC_Shew_Intact_26Sep14_Bane_C2Column3.raw";
        //private const string TestResultFile = @"D:\\MassSpecFiles\\QC_Shew_Intact_26Sep14_Bane_C2Column3_IcTda.tsv";
        public const string TestRawFile = @"\\protoapps\UserData\Sangtae\Yufeng\raw\yufeng_column_test2.raw";

        private const int MinCharge = 2;
        private const int MaxCharge = 50;
        private const int NumBits = 27;

        private static readonly IEnumerable<int> Ms2ScanNums = new[] { 46454, 46475, 46484, 46506, 46562, 46661 };

        [Test]
        public void TestExtractMs2FromMass()
        {
            var proteinMass = 43875.202027d;
            var run = PbfLcMsRun.GetLcMsRun(TestRawFile, MassSpecDataType.XCaliburRun, 0, 0.0);
            var ms1ScanNumbers = run.GetMs1ScanVector();

            var csm = new ChargeLcScanMatrix(run, NumBits, MinCharge, MaxCharge, false);

            var ms2Nums = csm.GetMatchingMs2ScanNums(proteinMass);

            foreach (var ms2 in ms2Nums)
            {
                Console.WriteLine(ms2 + "\t" + Ms2ScanNums.Contains(ms2));
            }
                

            
        }




    }
}
