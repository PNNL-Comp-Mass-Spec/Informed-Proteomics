using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Execution;
using InformedProteomics.TopDown.Scoring;
using NUnit.Framework;
using InformedProteomics.Backend.SequenceTag;

namespace InformedProteomics.Test
{
    [TestFixture]
    public class TestSequenceTagFinder
    {
        [Test]
        public void TestSequenceTag()
        {
            const string TestRawFile = @"D:\\Vlad_TopDown\\raw\\yufeng_column_test2.raw";
            const string TestResultFile = @"D:\\Vlad_TopDown\\results\\yufeng_column_test2_IcTda.tsv";

            var tsvParser = new TsvFileParser(TestResultFile);
            var headerList = tsvParser.GetHeaders();
            var tsvData = tsvParser.GetAllData();

            var scanNum = 46562;

            var seqStr = "";
            var modStr = "";
            double qValue = 0.0;

            var rawReader = new XCaliburReader(TestRawFile);

            Spectrum spectrum = rawReader.ReadMassSpectrum(scanNum);
            Sequence sequence = IdentifiedSequenceTag.GenerateSequence(seqStr, modStr);
            var tagFinder = new SequenceTagFinder(spectrum);

            //var existingTags = tagFinder.ExtractExistingSequneceTags(sequence);
            //Console.Write(scanNum + "\t" + existingTags.Count);

            foreach (var tag in tagFinder.FindSequenceTags())
            {
                Console.WriteLine(tag.GetTagStrings());
            }
                

        }    
    }
}
