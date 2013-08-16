using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Enum;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.IMS;
using InformedProteomics.Backend.IMSScoring;
using MultiDimensionalPeakFinding;
using NUnit.Framework;
using UIMFLibrary;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Test
{
    [TestFixture]
    internal class TestIcIms
    {
        [Test]
        public void TestDbSearchMultiThreading()
        {
            // Search parameters
            const int minPeptideLength = 6;
            const int maxPeptideLength = 30;
            const int numTolerableTermini = 2;
            const int numMissedCleavages = 1;
            const int numMaxModsPepPeptide = 0;
            const int minPrecursorCharge = 1;
            const int maxPrecursorCharge = 4;
            var enzyme = Enzyme.Trypsin;

            var sw = new System.Diagnostics.Stopwatch();

            // Configure amino acids
            //var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var fixCarbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C', SequenceLocation.Everywhere, true);
            var searchModifications = new List<SearchModification> { fixCarbamidomethylC };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPepPeptide);

            // Read database
//            const string dbFilePath = @"..\..\..\TestFiles\BSA.fasta";
            const string dbFilePath = @"..\..\..\TestFiles\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            Console.Write("Reading a database...");
            Console.Out.Flush();
            sw.Reset();
            sw.Start();
            var targetDb = new FastaDatabase(dbFilePath);
            var decoyDb = targetDb.Decoy(enzyme);
            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);

            // Initialize IMS data
            //const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            const string uimfFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo_Precursors_Removed_Collision_Energy_Collapsed.UIMF";
            var imsData = new ImsDataCached(uimfFilePath);

            // Find features
            Console.Write("Creating precursor features...");
            Console.Out.Flush();
            sw.Reset();
            sw.Start();
            //imsData.CreatePrecursorFeaturesMultiThreads();
            imsData.CreatePrecursorFeatures();
            sw.Stop();
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);

            Console.Write("Creating fragment features...");
            Console.Out.Flush();
            sw.Reset();
            sw.Start();
            //imsData.CreateFragmentFeaturesMultiThreads();
            imsData.CreateFragmentFeatures();
            sw.Stop();
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);
            Console.Out.Flush();

            // Initialize scoring
            const string paramFile = @"..\..\..\TestFiles\HCD_train.mgf_para.txt";
            var imsScorerFactory = new ImsScorerFactory(paramFile);

            var writer = TextWriter.Synchronized(new StreamWriter(Path.ChangeExtension(uimfFilePath, "tsv")));
            writer.WriteLine(
            "Peptide\tPrecursorMz\tPrecursorCharge\tLCBegin\tLCEnd\tIMSBegin\tIMSEnd\tApexLC\tApexIMS\tSumIntensities\tNumPoints\tPrecursorScore\tScore\tIsDecoy");

            for (var isDecoy = 0; isDecoy < 2; isDecoy++)
            {
                var db = isDecoy == 0 ? targetDb : decoyDb;

                var indexedDb = new IndexedDatabase(db);
                var numPeptides = indexedDb.NumSequences(minPeptideLength, maxPeptideLength, numTolerableTermini, numMissedCleavages, enzyme);
                Console.WriteLine("Scoring {0} {1} peptides", numPeptides, isDecoy == 0 ? "target" : "decoy");
                var decoy = isDecoy;
                long numProcessedPeptides = 0;
                Parallel.ForEach(
                    indexedDb.SequencesAsStrings(
                        minPeptideLength, maxPeptideLength, numTolerableTermini, numMissedCleavages, enzyme),
                    annotation =>
                    {
                        // annotation: pre + "." + peptide + "." + post (e.g. R.PEPTIDER.G)
                        var peptide = annotation.Substring(2, annotation.Length - 4);

                        var seqGraph = new SequenceGraph(aaSet, peptide);
                        foreach (var scoringGraph in seqGraph.GetScoringGraphs())
                        {
                            scoringGraph.RegisterImsData(imsData, imsScorerFactory);
                            for (var precursorCharge = minPrecursorCharge; precursorCharge <= maxPrecursorCharge; precursorCharge++)
                            {
                                var precursorMz = scoringGraph.GetPrecursorIon(precursorCharge).GetMz();
                                if (precursorMz > imsData.MaxPrecursorMz || precursorMz < imsData.MinPrecursorMz)
                                    continue;
                                var best = scoringGraph.GetBestFeatureAndScore(precursorCharge);
                                var feature = best.Item1;
                                if (best.Item1 != null && best.Item3 > -50)
                                {
                                    writer.WriteLine(
                                        "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}"
                                        , annotation
                                        , precursorMz
                                        , precursorCharge
                                        , feature.ScanLcStart
                                        , feature.ScanLcStart + feature.ScanLcLength - 1
                                        , feature.ScanImsStart
                                        , feature.ScanImsStart + feature.ScanImsLength - 1
                                        , feature.ScanLcStart + feature.ScanLcRepOffset
                                        , feature.ScanImsStart + feature.ScanImsRepOffset
                                        , feature.SumIntensities
                                        , feature.NumPoints
                                        , best.Item2 // precursorScore
                                        , best.Item3 // score
                                        , decoy
                                        );
                                }
                            }
                        }
                        if(++numProcessedPeptides % 10000 == 0)
                        {
                            Console.WriteLine("Database search progress: {0}%", (double)numProcessedPeptides/numPeptides*100);
                        }
                    }
                    );
            }
            writer.Close();
        }

        [Test]
        public void TestParallelFeatureDetection()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            var uimfUtil = new UimfUtil(uimfFilePath);
            var featureDetectionUtil = new FeatureDetectionUtil(uimfFilePath, 11, 4);
            var tolerance = new Tolerance(25);
            const double minPrecursorMz = 400.0;
            const double maxPrecursorMz = 1000.0;
            var minTargetBin = uimfUtil.GetBinFromMz(400.0);
            var maxTargetBin = uimfUtil.GetBinFromMz(1000.0);

            var sw = new System.Diagnostics.Stopwatch();
            double sec = 0;

            // Parallel version
            sw.Reset();
            sw.Start();
            var targetMzList = Enumerable.Range(minTargetBin, maxTargetBin - minTargetBin + 1).Select(uimfUtil.GetMzFromBin).ToList();
            featureDetectionUtil.GetFeatureStatistics(targetMzList, tolerance.GetValue(), DataReader.FrameType.MS1, tolerance.GetUnit());
            sw.Stop();
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Parallel: {0:f4} sec", sec);

            // Serial version
            sw.Reset();
            sw.Start();
            var imsData = new ImsDataCached(uimfFilePath, minPrecursorMz, maxPrecursorMz, 10, 1500, tolerance, tolerance);
            imsData.CreatePrecursorFeatures();
            sw.Stop();
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Serial: {0:f4} sec", sec);
        }



        [Test]
        public void TestParallelization()
        {
            var sw = new System.Diagnostics.Stopwatch();
            sw.Start();

            // Search parameters
            const int minPeptideLength = 7;
            const int maxPeptideLength = 30;
            const int numTolerableTermini = 2;
            const int numMissedCleavages = 1;
            const int numMaxModsPepPeptide = 0;
            var enzyme = Enzyme.Trypsin;
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\SuffixArray\BSA.fasta";
            //const string dbFilePath = @"C:\cygwin\home\kims336\Data\SuffixArray\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";

            var db = new FastaDatabase(dbFilePath);
            var indexedDb = new IndexedDatabase(db);

            const string outputPath = @"C:\cygwin\home\kims336\Data\SuffixArray\test.txt";
            var writer = new StreamWriter(outputPath);

            //var annotations = indexedDb.SequencesAsStrings(
            //    minPeptideLength, maxPeptideLength, numTolerableTermini, numMissedCleavages, enzyme).ToArray();

            var histThreads = new Dictionary<int, int>();


            Parallel.ForEach(
                indexedDb.SequencesAsStrings(
                    minPeptideLength, maxPeptideLength, numTolerableTermini, numMissedCleavages, enzyme),
                //annotations,
                annotation =>
                {
                    writer.WriteLine(annotation);
                    //Console.WriteLine("{0}\t{1}", annotation, Thread.CurrentThread.ManagedThreadId);
                    var threadId = Thread.CurrentThread.ManagedThreadId;
                    if (histThreads.ContainsKey(threadId)) histThreads[threadId]++;
                    else histThreads[threadId] = 1;
                }
                );

            //foreach (var annotation in
            //    indexedDb.SequencesAsStrings(
            //    minPeptideLength, maxPeptideLength, numTolerableTermini, numMissedCleavages, enzyme)
            //    )
            //{
            //    writer.WriteLine(annotation);
            //    var threadId = Thread.CurrentThread.ManagedThreadId;
            //    if (histThreads.ContainsKey(threadId)) histThreads[threadId]++;
            //    else histThreads[threadId] = 1;
            //}

            foreach (var pair in histThreads)
            {
                Console.WriteLine("Thread {0} occurred {1} times", pair.Key, pair.Value);
            }

            writer.Close();

            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"{0:f4} sec", sec);
        }

        [Test]
        public void TestDbSearch()
        {
            // Search parameters
            const int minPeptideLength = 7;
            const int maxPeptideLength = 30;
            const int numTolerableTermini = 2;
            const int numMissedCleavages = 1;
            const int numMaxModsPepPeptide = 0;
            const int minPrecursorCharge = 2;
            const int maxPrecursorCharge = 2;

            const double minPrecursorMz = 600;
            const double maxPrecursorMz = 610;
            const double minFragmentMz = 500;
            const double maxFragmentMz = 510;
            var precursorTolerance = new Tolerance(25);
            var fragmentTolerance = new Tolerance(25);

            var enzyme = Enzyme.Trypsin;

            // Configure amino acids
            //var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var fixCarbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C', SequenceLocation.Everywhere, true);
            var searchModifications = new List<SearchModification> { fixCarbamidomethylC };
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPepPeptide);

            // Initialize IMS data
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            //const string uimfFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo_Precursors_Removed_Collision_Energy_Collapsed.UIMF";
            //var imsData = new ImsDataCached(uimfFilePath, minPrecursorMz, maxPrecursorMz, minFragmentMz, maxFragmentMz, precursorTolerance, fragmentTolerance);
            var imsData = new ImsDataCached(uimfFilePath);

            var sw = new System.Diagnostics.Stopwatch();
            Console.Out.Flush();
            Console.Write("Creating precursor features...");
            sw.Start();
            imsData.CreatePrecursorFeatures();
            sw.Stop();
            var sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);

            Console.Write("Creating fragment features...");
            sw.Start();
            imsData.CreateFragmentFeatures();
            sw.Stop();
            sec = (double)sw.ElapsedTicks / (double)System.Diagnostics.Stopwatch.Frequency;
            Console.WriteLine(@"Done. {0:f4} sec", sec);
            Console.Out.Flush();

            // Initialize scoring
            const string paramFile = @"..\..\..\TestFiles\HCD_train.mgf_para.txt";
            var imsScorerFactory = new ImsScorerFactory(paramFile);

            // Read database
            const string dbFilePath = @"C:\cygwin\home\kims336\Data\SuffixArray\BSA.fasta";
            //const string dbFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\H_sapiens_Uniprot_SPROT_2013-05-01_withContam.fasta";
            var targetDb = new FastaDatabase(dbFilePath);
            var decoyDb = targetDb.Decoy(enzyme);

            var output = new StringBuilder();
            using (var writer = new StreamWriter(Path.ChangeExtension(uimfFilePath, "tsv")))
            {
                writer.WriteLine(
                "Peptide\tPrecursorMz\tPrecursorCharge\tLCBegin\tLCEnd\tIMSBegin\tIMSEnd\tApexLC\tApexIMS\tSumIntensities\tNumPoints\tPrecursorScore\tScore\tIsDecoy");

                for (var isDecoy = 0; isDecoy < 2; isDecoy++)
                {
                    var db = isDecoy == 0 ? targetDb : decoyDb;

                    var indexedDb = new IndexedDatabase(db);
                    Console.WriteLine("Scoring {0} peptides", isDecoy == 0 ? "target" : "decoy");
                    foreach(var annotation in indexedDb.SequencesAsStrings(
                            minPeptideLength, maxPeptideLength, numTolerableTermini, numMissedCleavages, enzyme))
                    {
                        // annotation: pre + "." + peptide + "." + post (e.g. R.PEPTIDER.G)
                        var peptide = annotation.Substring(2, annotation.Length - 4);

                        //var scoringGraph = seqGraph.GetScoringGraph(0);

                        var seqGraph = new SequenceGraph(aaSet, peptide);
                        foreach (var scoringGraph in seqGraph.GetScoringGraphs())
                        {
                            scoringGraph.RegisterImsData(imsData, imsScorerFactory);
                            for (var precursorCharge = minPrecursorCharge; precursorCharge <= maxPrecursorCharge; precursorCharge++)
                            {
                                var precursorMz = scoringGraph.GetPrecursorIon(precursorCharge).GetMz();
                                if (precursorMz > imsData.MaxPrecursorMz || precursorMz < imsData.MinPrecursorMz)
                                    continue;
                                var best = scoringGraph.GetBestFeatureAndScore(precursorCharge);
                                var feature = best.Item1;
                                if (best.Item1 != null && best.Item3 > -20)
                                {
                                    writer.WriteLine(
                                        "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}"
                                        , annotation
                                        , precursorMz
                                        , precursorCharge
                                        , feature.ScanLcStart
                                        , feature.ScanLcStart + feature.ScanLcLength - 1
                                        , feature.ScanImsStart
                                        , feature.ScanImsStart + feature.ScanImsLength - 1
                                        , feature.ScanLcStart + feature.ScanLcRepOffset
                                        , feature.ScanImsStart + feature.ScanImsRepOffset
                                        , feature.SumIntensities
                                        , feature.NumPoints
                                        , best.Item2 // precursorScore
                                        , best.Item3 // score
                                        , isDecoy
                                        );
                                }
                            }
                        }
                    }
                }
            }
        }

        [Test]
        public void TestImsScoring()
        {
            //const string uimfFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo_Collision_Energy_Collapsed.UIMF";
            const string uimfFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo_Precursors_Removed_Collision_Energy_Collapsed.UIMF";
            //const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            //var imsData = new ImsDataCached(uimfFilePath);
            var imsData = new ImsDataCached(uimfFilePath, 300.0, 2000.0, 10.0, 2500.0,
                new Tolerance(25, DataReader.ToleranceType.PPM), new Tolerance(25, DataReader.ToleranceType.PPM));
            const string paramFile = @"..\..\..\TestFiles\HCD_train.mgf_para.txt";
            var imsScorerFactory = new ImsScorerFactory(paramFile);

            //const string targetPeptide = "CCAADDKEACFAVEGPK";
            //const string targetPeptide = "ECCHGDLLECADDRADLAK";
            //const string targetPeptide = "VTLTCVAPLSGVDFQLR";       // charge 2 must show up at (196,173)
            const string targetPeptide = "CSPHLVLSALTSDNHGATYAFSGTHYWR";
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);

            var seqGraph = new SequenceGraph(aaSet, targetPeptide);
            var scoringGraph = seqGraph.GetScoringGraph(0);
            //foreach (var composition in scoringGraph.GetCompositions())
            //{
            //    Console.WriteLine(composition);
            //}
            scoringGraph.RegisterImsData(imsData, imsScorerFactory);
            for (var precursorCharge = 2; precursorCharge <= 4; precursorCharge++)
            {
                var best = scoringGraph.GetBestFeatureAndScore(precursorCharge);
                Console.WriteLine("PrecursorMz: " + scoringGraph.GetPrecursorIon(precursorCharge).GetMz());
                Console.WriteLine("Charge: " + precursorCharge);
                Console.WriteLine("Feature: " + best.Item1);
                Console.WriteLine("PrecursorScore: " + best.Item2);
                Console.WriteLine("Score: " + best.Item3);
                Console.WriteLine();
            }
        }
        [Test]
        public void TestBsaSearch()
        {
            var oxM = new SearchModification(Modification.Oxidation, 'M', SequenceLocation.Everywhere, false);
            var fixCarbamidomethylC = new SearchModification(Modification.Carbamidomethylation, 'C', SequenceLocation.Everywhere, true);

            var searchModifications = new List<SearchModification> { fixCarbamidomethylC };
            const int numMaxModsPepPeptide = 2;
            var aaSet = new AminoAcidSet(searchModifications, numMaxModsPepPeptide);

            const string uimfFilePath =
                @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            //const string uimfFilePath =
            //    @"C:\cygwin\home\kims336\Data\IMS_Ian_BSA_New\BSA_digest_IMS03_AgTOF01_1ugperml_2_tenths_higher_gain.UIMF";
            //const string uimfFilePath =
            //    @"C:\cygwin\home\kims336\Data\IMS_Ian_BSA\BSA digest IMS03 AgQTOF01 1ugperml.UIMF";
            var imsData = new ImsDataCached(uimfFilePath, 300.0, 2000.0, 10.0, 2500.0,
                new Tolerance(25, DataReader.ToleranceType.PPM), new Tolerance(25, DataReader.ToleranceType.PPM));

            const string paramFile = @"..\..\..\TestFiles\HCD_train.mgf_para.txt";
            var imsScorerFactory = new ImsScorerFactory(paramFile);


            var dbFilePaths = new string[] 
                {
                    @"..\..\..\TestFiles\BSAPeptides_ST.txt",
                    @"..\..\..\TestFiles\BSA_ReversePeptides_ST.txt"
                };

            var isDecoy = 0;

            //using (var writer = new StreamWriter(@"C:\cygwin\home\kims336\Data\IMS_Sarc\ICResults.txt"))
            //using (var writer = new StreamWriter(@"..\..\..\TestFiles\IC_BSA_Results.txt"))
            //using (var writer = new StreamWriter(@"C:\cygwin\home\kims336\Data\IMS_BSA\IC_BSA_NewResults.txt"))
            using (var writer = new StreamWriter(Path.GetDirectoryName(uimfFilePath)+@"\IC_BSA_OldResults.txt"))
            {
                writer.WriteLine(
                    "Peptide\tPrecursorMz\tPrecursorCharge\tLCBegin\tLCEnd\tIMSBegin\tIMSEnd\tApexLC\tApexIMS\tSumIntensities\tNumPoints\tPrecursorScore\tScore\tIsDecoy");
                foreach (var dbFilePath in dbFilePaths)
                {
                    foreach (var annotation in File.ReadLines(dbFilePath))
                    {
                        var peptide = annotation.Substring(2, annotation.Length - 4);

                        //var scoringGraph = seqGraph.GetScoringGraph(0);
                        var seqGraph = new SequenceGraph(aaSet, peptide);
                        foreach (var scoringGraph in seqGraph.GetScoringGraphs())
                        {
                            scoringGraph.RegisterImsData(imsData, imsScorerFactory);
                            for (var precursorCharge = 1; precursorCharge <= 5; precursorCharge++)
                            {
                                var precursorMz = scoringGraph.GetPrecursorIon(precursorCharge).GetMz();
                                if (precursorMz > imsData.MaxPrecursorMz || precursorMz < imsData.MinPrecursorMz)
                                    continue;
                                var best = scoringGraph.GetBestFeatureAndScore(precursorCharge);
                                var feature = best.Item1;
                                if (best.Item1 != null && best.Item3 > -20)
                                {
                                    writer.WriteLine(
                                        "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}"
                                        , annotation
                                        , precursorMz
                                        , precursorCharge
                                        , feature.ScanLcStart
                                        , feature.ScanLcStart + feature.ScanLcLength - 1
                                        , feature.ScanImsStart
                                        , feature.ScanImsStart + feature.ScanImsLength - 1
                                        , feature.ScanLcStart + feature.ScanLcRepOffset
                                        , feature.ScanImsStart + feature.ScanImsRepOffset
                                        , feature.SumIntensities
                                        , feature.NumPoints
                                        , best.Item2    // precursorScore
                                        , best.Item3    // score
                                        , isDecoy
                                        );
                                }
                            }
                        }
                    }
                    ++isDecoy;
                }
            }
        }
        private class Match
        {
            public string Peptide { get; set; }
            public double PrecursorMz { get; set; }
            public int PrecursorCharge { get; set; }
            public int LcBegin { get; set; }
            public int LcEnd { get; set; }
            public int ImsBegin { get; set; }
            public int ImsEnd { get; set; }
            public int ApexLc { get; set; }
            public int ApexIms { get; set; }
            public double Score { get; set; }
        }

        [Test]
        public void ExtractBestMatchPerFeature()
        {
            var isotopeMass = Atom.Get("13C").Mass - Atom.Get("C").Mass;

            const string resultFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\ICResults_Human.txt";
            const string outputFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\ICResults_Filtered_Human.tsv";

            var bestFeatureDictionary = new SortedDictionary<double, List<Match>>();

            using (var writer = new StreamWriter(outputFilePath))
            {
                foreach (var result in File.ReadLines(resultFilePath))
                {
                    if (result.StartsWith("Peptide"))
                    {
                        writer.WriteLine(result);
                        continue;
                    }
                    var token = result.Split('\t');
                    var match = new Match
                        {
                            Peptide = token[0],
                            PrecursorMz = double.Parse(token[1]),
                            PrecursorCharge = int.Parse(token[2]),
                            LcBegin = int.Parse(token[3]),
                            LcEnd = int.Parse(token[4]),
                            ImsBegin = int.Parse(token[5]),
                            ImsEnd = int.Parse(token[6]),
                            ApexLc = int.Parse(token[7]),
                            ApexIms = int.Parse(token[8]),
                            Score = double.Parse(token[11])
                        };
                    for (var isotope = 0; isotope < 4; isotope++)   // consider up to 4th isotope
                    {
                        var precursorMz = match.PrecursorMz + isotopeMass * isotope / match.PrecursorCharge;
                        var tolerance = precursorMz * 5 / 1e6;  // 5 ppm tolerance
                        foreach (var m in bestFeatureDictionary.Where(
                            entry => entry.Key > precursorMz - tolerance && entry.Key < precursorMz + tolerance))
                        {
                        }
                    }
                }
            }
        }

        [Test]
        public void ResumeSarcSearch()
        {
            const string uimfFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo_Collision_Energy_Collapsed.UIMF";
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);

            var imsData = new ImsDataCached(uimfFilePath, 300.0, 2000.0, 10.0, 2500.0,
                new Tolerance(25, DataReader.ToleranceType.PPM), new Tolerance(25, DataReader.ToleranceType.PPM));

            const string paramFile = @"..\..\..\TestFiles\HCD_train.mgf_para.txt";
            var imsScorerFactory = new ImsScorerFactory(paramFile);

            var dbFilePaths = new string[] 
                {
                    @"C:\cygwin\home\kims336\Data\IMS_Sarc\HumanPeptides.txt",
                    @"C:\cygwin\home\kims336\Data\IMS_Sarc\HumanPeptides_Reverse.txt"
                };

            const string outputFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\ICResults.txt";
            var lastLine = File.ReadLines(outputFilePath).Last();

            var lastAnnotation = lastLine.Split('\t')[0];

            const string newOutputFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\ICResults_Resumed.txt";

            var newPeptide = false;
            var lineNum = 0;
            const int isDecoy = 1;
            using (var writer = new StreamWriter(newOutputFilePath))
            {
                foreach (var annotation in File.ReadLines(dbFilePaths[1]))
                {
                    ++lineNum;
                    if (lineNum%10000 == 1)
                    {
                        Console.WriteLine("Processing {0} peptides...", lineNum);
                    }

                    if (!newPeptide && annotation.Equals(lastAnnotation))
                    {
                        newPeptide = true;
                        Console.WriteLine("Found " + annotation);
                    }
                    if (!newPeptide)
                        continue;
                    var peptide = annotation.Substring(2, annotation.Length - 4);

                    var seqGraph = new SequenceGraph(aaSet, peptide);
                    foreach (var scoringGraph in seqGraph.GetScoringGraphs())
                    {
                        scoringGraph.RegisterImsData(imsData, imsScorerFactory);
                        for (var precursorCharge = 1; precursorCharge <= 5; precursorCharge++)
                        {
                            double precursorMz = scoringGraph.GetPrecursorIon(precursorCharge).GetMz();
                            if (precursorMz > imsData.MaxPrecursorMz || precursorMz < imsData.MinPrecursorMz)
                                continue;
                            var best = scoringGraph.GetBestFeatureAndScore(precursorCharge);
                            var feature = best.Item1;
                            if (best.Item1 != null && best.Item2 > 0)
                            {
                                writer.WriteLine(
                                    "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}"
                                    , annotation
                                    , precursorMz
                                    , precursorCharge
                                    , feature.ScanLcStart
                                    , feature.ScanLcStart + feature.ScanLcLength - 1
                                    , feature.ScanImsStart
                                    , feature.ScanImsStart + feature.ScanImsLength - 1
                                    , feature.ScanLcStart + feature.ScanLcRepOffset
                                    , feature.ScanImsStart + feature.ScanImsRepOffset
                                    , feature.SumIntensities
                                    , feature.NumPoints
                                    , best.Item2
                                    , isDecoy
                                    );
                            }
                        }
                    }
                }
            }
        }

        [Test]
        public void TestSarcSearch()
        {
            const string uimfFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo_Precursors_Removed_Collision_Energy_Collapsed.UIMF";

            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);

            var imsData = new ImsDataCached(uimfFilePath, 300.0, 2000.0, 10.0, 2500.0,
                new Tolerance(25, DataReader.ToleranceType.PPM), new Tolerance(25, DataReader.ToleranceType.PPM));

            const string paramFile = @"..\..\..\TestFiles\HCD_train.mgf_para.txt";
            var imsScorerFactory = new ImsScorerFactory(paramFile);

            Console.WriteLine("Finished reading the UIMF file.");

            var dbFilePaths = new string[] 
                {
                    @"C:\cygwin\home\kims336\Data\IMS_Sarc\HumanPeptides.txt",
                    @"C:\cygwin\home\kims336\Data\IMS_Sarc\HumanPeptides_Reverse.txt"
                };

            var isDecoy = 0;

            const string outputFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\ICResults_PrecursorRemoved.txt";
            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine(
                    "Peptide\tPrecursorMz\tPrecursorCharge\tLCBegin\tLCEnd\tIMSBegin\tIMSEnd\tApexLC\tApexIMS\tSumIntensities\tNumPoints\tScore\tIsDecoy");
                foreach (var dbFilePath in dbFilePaths)
                {
                    var lineNum = 0;
                    foreach (var annotation in File.ReadLines(dbFilePath))
                    {
                        ++lineNum;
                        if (lineNum%10000 == 1)
                        {
                            Console.WriteLine("Processing {0} peptides...", lineNum);
                        }
                        var peptide = annotation.Substring(2, annotation.Length - 4);

                        var seqGraph = new SequenceGraph(aaSet, peptide);
                        foreach (var scoringGraph in seqGraph.GetScoringGraphs())
                        {
                            scoringGraph.RegisterImsData(imsData, imsScorerFactory);
                            for (var precursorCharge = 1; precursorCharge <= 5; precursorCharge++)
                            {
                                double precursorMz = scoringGraph.GetPrecursorIon(precursorCharge).GetMz();
                                if (precursorMz > imsData.MaxPrecursorMz || precursorMz < imsData.MinPrecursorMz)
                                    continue;
                                var best = scoringGraph.GetBestFeatureAndScore(precursorCharge);
                                var feature = best.Item1;
                                if (best.Item1 != null && best.Item2 > 0)
                                {
                                    writer.WriteLine(
                                        "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}"
                                        , annotation
                                        , precursorMz
                                        , precursorCharge
                                        , feature.ScanLcStart
                                        , feature.ScanLcStart + feature.ScanLcLength - 1
                                        , feature.ScanImsStart
                                        , feature.ScanImsStart + feature.ScanImsLength - 1
                                        , feature.ScanLcStart + feature.ScanLcRepOffset
                                        , feature.ScanImsStart + feature.ScanImsRepOffset
                                        , feature.SumIntensities
                                        , feature.NumPoints
                                        , best.Item2
                                        , isDecoy
                                        );
                                }
                            }
                        }
                    }
                    ++isDecoy;
                }
            }
        }

        [Test]
        public void TestImsFeatureFinding()
        {
            const string uimfFilePath = @"C:\cygwin\home\kims336\Data\IMS_Sarc\SarcCtrl_P21_1mgml_IMS6_AgTOF07_210min_CID_01_05Oct12_Frodo_Collision_Energy_Collapsed.UIMF";
            //const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            var imsData = new ImsDataCached(uimfFilePath);

            //Console.WriteLine("Generating precursor features (MinMz: " + imsData.MinPrecursorMz + " MaxMz: " + imsData.MaxPrecursorMz + ")");
            //int numPrecursorFeatures = imsData.CreatePrecursorFeatures();
            //Console.WriteLine("TotalNumPrecursorFeatures: " + numPrecursorFeatures);            

            const string targetPeptide = "VTLTCVAPLSGVDFQLR";//EYANQFMWEYSTNYGQAPLSLLVSYTK  CCHGDLLECADDRADLAK
            var aaSet = new AminoAcidSet(Modification.Carbamidomethylation);
            var precursorComposition = aaSet.GetComposition(targetPeptide);
            for (int charge = 2; charge <= 2; charge++)
            {
                var precursorIon = new Ion(precursorComposition + Composition.H2O, charge);
                for (var isotopeIndex = 0; isotopeIndex <= 5; isotopeIndex++)
                {
                    double precursorMz = precursorIon.GetIsotopeMz(isotopeIndex);
                    FeatureSet precursorFeatures = imsData.GetPrecursorFeatures(precursorMz);
                    Console.WriteLine("PrecursorMz: {0}, Charge: {1}, Isotope: {2}\n", precursorMz, charge, isotopeIndex);
                    foreach (Feature precursorFeature in precursorFeatures)
                    {
                        Console.WriteLine(precursorFeature);
                        Console.WriteLine("LC Apex profile");
                        Console.WriteLine(string.Join("\t", precursorFeature.LcApexPeakProfile));
                        Console.WriteLine("IMS Apex profile");
                        Console.WriteLine(string.Join("\t", precursorFeature.ImsApexPeakProfile));
                        Console.WriteLine();
                    }
                }
            }

            //Console.WriteLine("Generating fragment features (MinMz: " + imsData.MinFragmentMz + " MaxMz: " + imsData.MaxFragmentMz + ")");
            //int numFragmentFeatures = imsData.CreateFragmentFeatures();
            //Console.WriteLine("TotalNumFragmentFeatures: " + numFragmentFeatures);

            //const string targetFragment = "AADDKEACFAVEGPK"; //"YGQAPLSLLVSYTK";
            //var fragmentComposition = aaSet.GetComposition(targetFragment);
            //    for (int charge = 2; charge <= 2; charge++)
            //{
            //    var y15C2 = new Ion(fragmentComposition + Composition.H2O, 2);
            //    double fragmentMz = y15C2.GetMz();
            //    FeatureSet fragmentFeatures = imsData.GetFragmentFeatures(fragmentMz);
            //    Console.WriteLine("Fragment: {0}, Charge: {1}", fragmentMz, charge);
            //    foreach (Feature fragmentFeature in fragmentFeatures)
            //    {
            //        Console.WriteLine(fragmentFeature);
            //        Console.WriteLine("LC Apex profile");
            //        Console.WriteLine(string.Join("\t", fragmentFeature.LcApexPeakProfile));
            //        Console.WriteLine("IMS Apex profile");
            //        Console.WriteLine(string.Join("\t", fragmentFeature.ImsApexPeakProfile));
            //        Console.WriteLine();
            //    }
            //}
        }

        [Test]
        public void TestGeneratingAllPrecursorXicsForEachMzBin()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";

            var imsData = new ImsData(uimfFilePath);
            var precursorXicMap = new Dictionary<int, FeatureSet>();

            int totalNumFeatures = 0;

            const double minMz = 400.0;
            const double maxMz = 2500.0;

            int minTargetBin = imsData.GetBinFromMz(minMz);
            int maxTargetBin = imsData.GetBinFromMz(maxMz);
            if (maxTargetBin >= imsData.GetNumberOfBins())
                maxTargetBin = imsData.GetNumberOfBins()-1;

            //Console.WriteLine("TargetMz: " + theoMz);
            Console.WriteLine("MinMz: " + minMz + " MaxMz: " + maxMz);
            Console.WriteLine("MinTargetBin: " + minTargetBin + " MaxTargetBin: " + maxTargetBin);
            for (int targetBin = minTargetBin; targetBin <= maxTargetBin; targetBin++)
            {
                double mz = imsData.GetMzFromBin(targetBin);
                FeatureSet featureSet = imsData.GetFeatures(targetBin, DataReader.FrameType.MS1);
                int numFeatures = featureSet.GetFeatures().Count();
                Console.WriteLine("Bin: " + targetBin + "\tMZ: " + mz + "\tNumFeatures: " + numFeatures);
                if (featureSet.GetFeatures().Any())
                    precursorXicMap.Add(targetBin, featureSet);
                totalNumFeatures += numFeatures;
            }
            Console.WriteLine("TotalNumPrecursorFeatures: " + totalNumFeatures);


        }

        [Test]
        public void TestGeneratingAllProductXicsForEachMzBin()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";

            var imsData = new ImsData(uimfFilePath);
            var fragmentXicMap = new Dictionary<int, FeatureSet>();

            int totalNumFragmentFeatures = 0;

            const double minMz = 0.0;
            const double maxMz = 2500.0;

            int minFragmentTargetMz = imsData.GetBinFromMz(minMz);
            int maxFragmentTargetBin = imsData.GetBinFromMz(maxMz);
            if (maxFragmentTargetBin >= imsData.GetNumberOfBins())
                maxFragmentTargetBin = imsData.GetNumberOfBins() - 1;

            Console.WriteLine("MinMz: " + minMz + " MaxMz: " + maxMz);
            Console.WriteLine("MinTargetBin: " + minFragmentTargetMz + " MaxTargetBin: " + maxFragmentTargetBin);
            for (int targetBin = minFragmentTargetMz; targetBin <= maxFragmentTargetBin; targetBin++)
            {
                double mz = imsData.GetMzFromBin(targetBin);
                FeatureSet featureSet = imsData.GetFeatures(targetBin, DataReader.FrameType.MS2);
                int numFeatures = featureSet.GetFeatures().Count();
                Console.WriteLine("Bin: " + targetBin + "\tMZ: " + mz + "\tNumFeatures: " + numFeatures);
                if (featureSet.GetFeatures().Any())
                    fragmentXicMap.Add(targetBin, featureSet);
                totalNumFragmentFeatures += numFeatures;
            }
            Console.WriteLine("TotalNumFragmentFeatures: " + totalNumFragmentFeatures);
        }

        [Test]
        public void TestSimple()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            var uimfReader = new DataReader(uimfFilePath);         
            Console.WriteLine("NumFrames: " + uimfReader.GetGlobalParameters().NumFrames);
            Console.WriteLine("NumScans: " + uimfReader.GetFrameParameters(1).Scans);
        }

        [Test]
        public void TestGeneratingSpectrum()
        {
            const string uimfFilePath = @"..\..\..\TestFiles\BSA_10ugml_IMS6_TOF03_CID_27Aug12_Frodo_Collision_Energy_Collapsed.UIMF";
            var dataReader = new DataReader(uimfFilePath);
            const int startFrameNumber = 579;
            const int endFrameNumber = 581;
            const int startScanNumber = 124;
            const int endScanNumber = 135;

            double[] mzArray;
            int[] intensityArray;
            dataReader.GetSpectrum(startFrameNumber, endFrameNumber, DataReader.FrameType.MS2, startScanNumber,
                                   endScanNumber, out mzArray, out intensityArray);

            var size = mzArray.Length;
            using (var writer = new StreamWriter(@"C:\cygwin\home\kims336\Data\IMS_BSA\testSpec.mgf"))
            {
                writer.WriteLine("BEGIN IONS");
                writer.WriteLine("TITLE=test");
                writer.WriteLine("PEPMASS=582.3447");
                writer.WriteLine("CHARGE=2+");
                const int windowSize = 10;
                for (var i = 0; i < size; i++)
                {
                    bool isPeak = true;
                    int rank = 0;
                    for (var j = i - windowSize; j < i + windowSize; j++)
                    {
                        if (j >= 0 && j < intensityArray.Length)
                        {
                            if (intensityArray[i] < intensityArray[j])
                            {
                                isPeak = false;
                                break;
                            }
                            if (mzArray[j] > mzArray[i] - 0.2 && mzArray[j] < mzArray[i] + 0.2)
                            {
                                if (intensityArray[j] <= intensityArray[i]) rank++;
                            }
                        }
                    }
                    //if (rank <= 1)
                    {
                        writer.WriteLine(mzArray[i] + "\t" + intensityArray[i]);
                    }
                }
                writer.WriteLine("END IONS");
            }
        }
    }
}
