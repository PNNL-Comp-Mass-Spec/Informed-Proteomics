using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Database;

namespace InformedProteomics.Backend.Utils
{
    public class FdrCalculator
    {
        private IList<string> _headers;
        private string[] _results;

        public FdrCalculator(string targetResultFilePath, string decoyResultFilePath)
        {
            if (!CalculateQValues(targetResultFilePath, decoyResultFilePath))
            {
                throw new Exception("Illegal file format at FdrCalculator");
            }
            if (!CalculatePepQValues(targetResultFilePath, decoyResultFilePath))
            {
                throw new Exception("Illegal file format at FdrCalculator");
            }
        }

        public void WriteTo(string outputFilePath, bool includeDecoy = false)
        {
            var sequenceIndex = _headers.IndexOf("ProteinName");

            using (var writer = new StreamWriter(outputFilePath))
            {
                writer.WriteLine(string.Join("\t", _headers));
                foreach (var r in _results)
                {
                    if (!includeDecoy && r.Split('\t')[sequenceIndex].StartsWith(FastaDatabase.DecoyProteinPrefix))
                    {
                        continue;
                    }
                    writer.WriteLine(r);
                }
            }
        }

        private bool CalculateQValues(string targetResultFilePath, string decoyResultFilePath)
        {
            var targetData = File.ReadAllLines(targetResultFilePath);
            var decoyData = File.ReadAllLines(decoyResultFilePath);

            if (targetData.Length <= 1 || decoyData.Length <= 1) return false;

            var targetHeaders = targetData[0].Split('\t');
            var decoyHeaders = decoyData[0].Split('\t');

            if (targetHeaders.Length != decoyHeaders.Length) return false;
            if (targetHeaders.Where((t, i) => !t.Equals(decoyHeaders[i])).Any()) return false;

            _headers = targetHeaders;

            var concatenated = decoyData.Skip(1).Concat(targetData.Skip(1)).ToArray();
            var scoreIndex = _headers.IndexOf("IcScore");
            if (scoreIndex < 0) scoreIndex = _headers.IndexOf("Score");
            if (scoreIndex < 0) scoreIndex = _headers.IndexOf("#MatchedFragments");
            var scanNumIndex = _headers.IndexOf("Scan");
            var proteinIndex = _headers.IndexOf("ProteinName");
            if (scoreIndex < 0 || scanNumIndex < 0 || proteinIndex < 0) return false;

            var distinctSorted = concatenated.OrderByDescending(r => Convert.ToDouble(r.Split('\t')[scoreIndex]))
                .GroupBy(r => Convert.ToDouble(r.Split('\t')[scanNumIndex]))
                .Select(grp => grp.First())
                .ToArray();

            // Calculate q values
            _headers = _headers.Concat(new[] {"QValue"}).ToArray();
            var numDecoy = 0;
            var numTarget = 0;
            var fdr = new double[distinctSorted.Length];
            for (var i = 0; i < distinctSorted.Length; i++)
            {
                var row = distinctSorted[i];
                var columns = row.Split('\t');
                var protein = columns[proteinIndex];
                if (protein.StartsWith(FastaDatabase.DecoyProteinPrefix)) numDecoy++;
                else numTarget++;
                fdr[i] = Math.Min(numDecoy/(double)numTarget, 1.0);
            }
            var qValue = new double[fdr.Length];
            qValue[fdr.Length - 1] = fdr[fdr.Length - 1];
            for (var i = fdr.Length - 2; i >= 0; i--)
            {
                qValue[i] = Math.Min(qValue[i + 1], fdr[i]);
            }

            _results = distinctSorted.Select((r, i) => r + "\t" + qValue[i]).ToArray();

            return true;
        }

        private bool CalculatePepQValues(string targetResultFilePath, string decoyResultFilePath)
        {
            var targetData = File.ReadAllLines(targetResultFilePath);
            var decoyData = File.ReadAllLines(decoyResultFilePath);

            if (targetData.Length <= 1 || decoyData.Length <= 1) return false;

            var targetHeaders = targetData[0].Split('\t');
            var decoyHeaders = decoyData[0].Split('\t');

            if (targetHeaders.Length != decoyHeaders.Length) return false;
            if (targetHeaders.Where((t, i) => !t.Equals(decoyHeaders[i])).Any()) return false;

            var concatenated = decoyData.Skip(1).Concat(targetData.Skip(1)).ToArray();
            var scoreIndex = _headers.IndexOf("IcScore");
            if (scoreIndex < 0) scoreIndex = _headers.IndexOf("Score");
            if (scoreIndex < 0) scoreIndex = _headers.IndexOf("#MatchedFragments");
            var sequenceIndex = _headers.IndexOf("Sequence");
            var scanNumIndex = _headers.IndexOf("Scan");
            var preIndex = _headers.IndexOf("Pre");
            var postIndex = _headers.IndexOf("Post");
            var proteinIndex = _headers.IndexOf("ProteinName");
            if (scoreIndex < 0 || sequenceIndex < 0 || preIndex < 0 || postIndex < 0 || proteinIndex < 0) return false;

            var distinctSorted = concatenated.OrderByDescending(r => Convert.ToDouble(r.Split('\t')[scoreIndex]))
                .GroupBy(r => Convert.ToDouble(r.Split('\t')[scanNumIndex]))
                .Select(grp => grp.First())
                .GroupBy(r => r.Split('\t')[preIndex] + r.Split('\t')[sequenceIndex] + r.Split('\t')[postIndex])
                .Select(grp => grp.First())
                .ToArray();

            // Calculate q values
            _headers = _headers.Concat(new[] { "PepQValue" }).ToArray();
            var numDecoy = 0;
            var numTarget = 0;
            var fdr = new double[distinctSorted.Length];
            var peptide = new string[distinctSorted.Length];
            for (var i = 0; i < distinctSorted.Length; i++)
            {
                var row = distinctSorted[i];
                var columns = row.Split('\t');
                var protein = columns[proteinIndex];
                var annotation = columns[preIndex] + "." + columns[sequenceIndex] + "." + columns[postIndex];
                if (protein.StartsWith(FastaDatabase.DecoyProteinPrefix)) numDecoy++;
                else numTarget++;
                fdr[i] = Math.Min(numDecoy / (double)numTarget, 1.0);
                peptide[i] = annotation;
            }

            var pepQValue = new double[fdr.Length];
            pepQValue[fdr.Length - 1] = fdr[fdr.Length - 1];
            for (var i = fdr.Length - 2; i >= 0; i--)
            {
                pepQValue[i] = Math.Min(pepQValue[i + 1], fdr[i]);
            }

            var annotationToPepQValue = new Dictionary<string, double>();
            for (var i = 0; i < distinctSorted.Length; i++)
            {
                annotationToPepQValue[peptide[i]] = pepQValue[i];
            }

            for (var i = 0; i < _results.Length; i++)
            {
                var columns = _results[i].Split('\t');
                var annotation = columns[preIndex] + "." + columns[sequenceIndex] + "." + columns[postIndex];
                _results[i] = _results[i] + "\t" + annotationToPepQValue[annotation];
            }

            return true;
        }


    }
}
