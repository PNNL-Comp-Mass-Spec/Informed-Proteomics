using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.TopDownViewer.Models;

namespace InformedProteomics.TopDownViewer
{
    public class IcFileReader
    {
        public IcParameters Config { get; private set; }

        public IcFileReader(string tsvFile, string paramFile)
        {
            Config = new IcParameters(paramFile);
            _tsvFile = tsvFile;
        }

        public List<ProteinId> Read()
        {
            var proteins = new Dictionary<string, ProteinId>();
            var file = File.ReadLines(_tsvFile);

            var lineCount = 0;
            foreach (var line in file)
            {
                lineCount++;
                if (lineCount == 1) continue; // first line
                var idData = new PrSm(line, Config);
                if (idData.QValue > QValueThreshold) continue;
                var sequence = idData.Protein;
                if (!proteins.ContainsKey(sequence))
                    proteins.Add(sequence, new ProteinId(idData));
                else proteins[sequence].Add(idData);
            }

            return proteins.Values.ToList();
        }

        private readonly string _tsvFile;
        private const double QValueThreshold = 0.01;
    }
}
