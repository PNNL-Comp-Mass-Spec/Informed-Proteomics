using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.TopDownViewer
{
    public class IcFileReader
    {
        public IcFileReader(string tsvFile, string rawFile)
        {
            _lcms = LcMsRun.GetLcMsRun(rawFile, MassSpecDataType.XCaliburRun, 0, 0);
            _tsvFile = tsvFile;
        }

        public List<PrSm> Read()
        {
            var prsms = new Dictionary<string, PrSm>();
            var file = File.ReadLines(_tsvFile);

            var lineCount = 0;
            foreach (var line in file)
            {
                lineCount++;
                if (lineCount == 1) continue; // first line
                var chargeStateData = new ChargeStateData(line, _lcms);
                var sequence = chargeStateData.Sequence;
                if (!prsms.ContainsKey(sequence))
                    prsms.Add(sequence, new PrSm(sequence));
                prsms[sequence].AddCharge(chargeStateData);
            }
            return prsms.Values.ToList();
        }

        private readonly LcMsRun _lcms;
        private readonly string _tsvFile;
    }
}
