using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Utils
{
    public class TsvFileParser
    {
        public TsvFileParser(string fileName)
        {
            FileName = fileName;
            Parse();
        }

        public string FileName { get; private set; }

        public IList<string> GetHeaders()
        {
            return _header;
        }

        public Dictionary<string, List<string>> GetAllData()
        {
            return _data;
        }

        public IList<string> GetData(string columnName)
        {
            List<string> columnData;
            return _data.TryGetValue(columnName, out columnData) ? columnData : null;
        }

        public ISet<string> GetPeptides(double pepQValueThreshold)
        {
            var peptideColumnIndex = -1;
            var pepQValueColumnIndex = -1;
            var proteinColumnIndex = -1;

            for (var i = 0; i < _header.Length; i++)
            {
                if (_header[i].Equals("Peptide") || _header[i].Equals("#Peptide")) peptideColumnIndex = i;
                if (_header[i].Equals("Protein")) proteinColumnIndex = i;
                if (_header[i].Equals("PepQValue")) pepQValueColumnIndex = i;
            }

            var peptideSet = new HashSet<string>();
            if (pepQValueColumnIndex < 0 || peptideColumnIndex < 0 || proteinColumnIndex < 0) return null;

            var peptides = _data[_header[peptideColumnIndex]];
            var proteins = _data[_header[proteinColumnIndex]];
            var pepQValues = _data[_header[pepQValueColumnIndex]];
            for (var i = 0; i < pepQValues.Count; i++)
            {
                var pepQValue = Convert.ToDouble(pepQValues[i]);
                if(pepQValue <= pepQValueThreshold && !proteins[i].StartsWith("XXX")) peptideSet.Add(peptides[i]);
            }
            return peptideSet;
        }

        private string[] _header;
        private Dictionary<string, List<string>> _data;
        private void Parse()
        {
            _data = new Dictionary<string, List<string>>();
            // parse header
            var firstRow = true;
            foreach (var line in File.ReadLines(FileName))
            {
                var token = line.Split('\t');
                if (firstRow)
                {
                    _header = new string[token.Length];
                    for (var i = 0; i < token.Length; i++)
                    {
                        _header[i] = token[i];
                        _data[_header[i]] = new List<string>();
                    }
                    firstRow = false;
                    continue;
                }

                if (token.Length != _header.Length) continue;
                for (var i = 0; i < token.Length; i++)
                {
                    _data[_header[i]].Add(token[i]);
                }
            }
        }
    }
}
