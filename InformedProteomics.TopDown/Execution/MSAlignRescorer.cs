using System;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Scoring;

namespace InformedProteomics.TopDown.Execution
{
    public class MsAlignRescorer
    {
        public MsAlignRescorer(string specFilePath, string msAlignFilePath, string outputFilePath, Tolerance tolerance, double ms2CorrThreshold = 0.7
            , int minProductIonCharge = 1, int maxProductIonCharge = 10)
        {
            var run = InMemoryLcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
            _topDownScorer = new InformedTopDownScorer(run, new AminoAcidSet(), minProductIonCharge, maxProductIonCharge, tolerance, ms2CorrThreshold);
            Rescore(msAlignFilePath, outputFilePath);
        }

        private readonly InformedTopDownScorer _topDownScorer;

        private void Rescore(string msAlignFilePath, string outputFilePath)
        {
            var parser = new TsvFileParser(msAlignFilePath);
            var sequences = parser.GetData("Peptide");
            var scanNums = parser.GetData("Scan(s)").Select(s => Convert.ToInt32(s)).ToArray();
            var charges = parser.GetData("Charge").Select(c => Convert.ToInt32(c)).ToArray();

            var rows = parser.GetRows();
            var headers = parser.GetHeaders();

            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine("{0}\t{1}", string.Join("\t", headers), IcScores.GetScoreNames());
                for (var i = 0; i < parser.NumData; i++)
                {
                    var row = rows[i];
                    var seqStr = SimpleStringProcessing.GetStringBetweenDots(sequences[i]);
                    if (seqStr == null || seqStr.Contains("(")) continue; //TODO: currently ignore ids with modifications

                    var composition = AASet.GetComposition(seqStr);
                    //var sequence = new Sequence(seqStr, AASet);
                    //if (sequence == null)
                    //{
                    //    Console.WriteLine("Ignore illegal sequence: {0}", seqStr);
                    //    continue;
                    //}
                    var charge = charges[i];
                    var scanNum = scanNums[i];

                    var scores = _topDownScorer.GetScores(AminoAcid.ProteinNTerm, seqStr, AminoAcid.ProteinCTerm, composition, charge, scanNum);
                    if (scores == null) continue;

                    writer.WriteLine("{0}\t{1}", row, scores);
                }
            }
        }

        private static readonly AminoAcidSet AASet = new AminoAcidSet();

    }
}
