using System;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;
using InformedProteomics.TopDown.Scoring;

namespace InformedProteomics.TopDown.Execution
{
    public class IcRescorer
    {
        public IcRescorer(string specFilePath, string icResultFilePath, string outputFilePath, AminoAcidSet aaSet, Tolerance tolerance, double ms2CorrThreshold = 0.7
            , int minProductIonCharge = 1, int maxProductIonCharge = 10)
        {
            var run = LcMsRun.GetLcMsRun(specFilePath, MassSpecDataType.XCaliburRun, 1.4826, 1.4826);
            _scorer = new InformedScorer(run, aaSet, minProductIonCharge, maxProductIonCharge, tolerance, ms2CorrThreshold);
            Rescore(icResultFilePath, outputFilePath);
        }

        private readonly InformedScorer _scorer;

        private void Rescore(string icResultFilePath, string outputFilePath)
        {
            var parser = new TsvFileParser(icResultFilePath);
            var sequences = parser.GetData("Sequence");
            var scanNums = parser.GetData("ScanNum").Select(s => Convert.ToInt32(s)).ToArray();
            var charges = parser.GetData("Charge").Select(c => Convert.ToInt32(c)).ToArray();
            var compositions = parser.GetData("Composition").Select(Composition.Parse).ToArray();
            var modIndex = parser.GetHeaders().IndexOf("Modifications");

            var rows = parser.GetRows();
            var headers = parser.GetHeaders();

            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine("{0}\t{1}", string.Join("\t", headers), IcScores.GetScoreNames());
                for (var i = 0; i < parser.NumData; i++)
                {
                    var row = rows[i];
                    var seqStr = sequences[i];
                    var charge = charges[i];
                    var scanNum = scanNums[i];
                    var composition = compositions[i];

                    var scores = _scorer.GetScores(seqStr, composition, charge, scanNum);

                    var token = row.Split('\t');
                    for (var j = 0; j < token.Length; j++)
                    {
                        if(j != modIndex) writer.Write(token[j]+"\t");
                        else writer.Write("["+scores.Modifications+"]"+"\t");
                    }
                    writer.WriteLine(scores);

                    //var scores = _scorer.GetScores(sequence, charge, scanNum);
                    //if (scores == null) continue;
                    //writer.WriteLine("{0}\t{1}", row, scores);
                }
            }
        }
    }
}
