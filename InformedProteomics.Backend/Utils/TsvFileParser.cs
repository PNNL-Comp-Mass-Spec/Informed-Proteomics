using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace InformedProteomics.Backend.Utils
{
    /// <summary>
    /// Parses the data in a tab-delimited file, caching the data in memory
    /// </summary>
    /// <remarks>
    /// The data is stored in memory both as full rows (one string per row), and parsed by column
    /// Column names are case-sensitive
    /// </remarks>
    public class TsvFileParser
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="filePath">File to read</param>
        /// <param name="delimeter">Delimeter, tab by default</param>
        public TsvFileParser(string filePath, char delimeter = '\t')
        {
            FileName = filePath;
            _delimeter = delimeter;

            _header = new Dictionary<int, string>();
            _rows = new List<string>();
            _data = new Dictionary<string, List<string>>();

            Parse();
        }

        /// <summary>
        /// Filename
        /// </summary>
        public string FileName { get; }

        /// <summary>
        /// Number of rows
        /// </summary>
        public int NumData => _rows.Count;

        /// <summary>
        /// Get the TSV file headers
        /// </summary>
        /// <returns></returns>
        public IList<string> GetHeaders()
        {
            return (from item in _header orderby item.Key select item.Value).ToList();
        }

        /// <summary>
        /// Get all data in the TSV
        /// </summary>
        /// <returns></returns>
        public Dictionary<string, List<string>> GetAllData()
        {
            return _data;
        }

        /// <summary>
        /// Get the data in the specified column
        /// </summary>
        /// <param name="columnName"></param>
        /// <returns></returns>
        public IList<string> GetData(string columnName)
        {
            return _data.TryGetValue(columnName, out var columnData) ? columnData : null;
        }

        /// <summary>
        /// Get all rows of data
        /// </summary>
        /// <returns></returns>
        public IList<string> GetRows()
        {
            return _rows;
        }

        /// <summary>
        /// Get peptides that pass the specified threshold
        /// </summary>
        /// <param name="pepQValueThreshold"></param>
        /// <returns></returns>
        public ISet<string> GetPeptides(double pepQValueThreshold)
        {
            FindPeptideHeaders(out var peptideColumnIndex, out var proteinColumnIndex, out var pepQValueColumnIndex, "PepQValue");

            if (peptideColumnIndex < 0 || proteinColumnIndex < 0 || pepQValueColumnIndex < 0) return null;

            var peptideSet = GetPeptidesAboveThreshold(peptideColumnIndex, proteinColumnIndex, pepQValueColumnIndex, pepQValueThreshold);
            return peptideSet;
        }

        /// <summary>
        /// Get peptides that pass the specified threshold
        /// </summary>
        /// <param name="qValueThreshold"></param>
        /// <returns></returns>
        public ISet<string> GetPeptidesAboveQValueThreshold(double qValueThreshold)
        {
            FindPeptideHeaders(out var peptideColumnIndex, out var proteinColumnIndex, out var qValueColumnIndex, "QValue");

            if (peptideColumnIndex < 0 || proteinColumnIndex < 0 || qValueColumnIndex < 0) return null;

            var peptideSet = GetPeptidesAboveThreshold(peptideColumnIndex, proteinColumnIndex, qValueColumnIndex, qValueThreshold);
            return peptideSet;
        }

        /// <summary>
        /// Map from header column index to header column name
        /// </summary>
        private readonly Dictionary<int, string> _header;

        private readonly List<string> _rows;

        private readonly Dictionary<string, List<string>> _data;

        private readonly char _delimeter;

        private int FindColumnIgnoreCase(Dictionary<int, string> headers, string headerName)
        {
            foreach (var i in headers.Keys)
            {
                if (headers[i].IndexOf(headerName, StringComparison.OrdinalIgnoreCase) >= 0)
                {
                    return i;
                }
            }

            return -1;

        }

        private void FindPeptideHeaders(out int peptideColumnIndex, out int proteinColumnIndex, out int qValueColumnIndex, string qValueColumnName)
        {
            peptideColumnIndex = -1;
            proteinColumnIndex = -1;
            qValueColumnIndex = -1;

            foreach (var i in _header.Keys)
            {
                if (_header[i].Equals("Peptide") || _header[i].Equals("#Peptide") || _header[i].Equals("Sequence"))
                    peptideColumnIndex = i;
                else if (_header[i].Equals("Protein") || _header[i].Equals("ProteinName"))
                    proteinColumnIndex = i;
                else if (_header[i].Equals(qValueColumnName))
                    qValueColumnIndex = i;
            }

            if (peptideColumnIndex < 0)
            {
                peptideColumnIndex = FindColumnIgnoreCase(_header, "peptide");
            }

            if (proteinColumnIndex < 0)
            {
                proteinColumnIndex = FindColumnIgnoreCase(_header, "protein");
            }

            if (qValueColumnIndex < 0)
            {
                qValueColumnIndex = FindColumnIgnoreCase(_header, qValueColumnName);
            }

            if (qValueColumnIndex < 0)
            {
                qValueColumnIndex = FindColumnIgnoreCase(_header, "qvalue");
            }
        }

        private ISet<string> GetPeptidesAboveThreshold(int peptideColumnIndex, int proteinColumnIndex, int qValueColumnIndex, double qValueThreshold)
        {
            var peptides = _data[_header[peptideColumnIndex]];
            var proteins = _data[_header[proteinColumnIndex]];
            var qValues = _data[_header[qValueColumnIndex]];

            var peptideSet = new HashSet<string>();

            for (var i = 0; i < qValues.Count; i++)
            {
                if (double.TryParse(qValues[i], out var qValue))
                {
                    if (qValue <= qValueThreshold && !proteins[i].StartsWith("XXX"))
                        peptideSet.Add(peptides[i]);
                }
            }

            return peptideSet;

        }

        private void Parse()
        {
            _rows.Clear();
            _data.Clear();
            _header.Clear();

            var headerParsed = false;

            using (var reader = new StreamReader(new FileStream(FileName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
            {
                while (!reader.EndOfStream)
                {
                    var line = reader.ReadLine();
                    if (string.IsNullOrWhiteSpace(line))
                        continue;

                    var token = line.Split(_delimeter);
                    if (!headerParsed)
                    {
                        for (var i = 0; i < token.Length; i++)
                        {
                            if (_data.ContainsKey(token[i]))
                            {
                                Console.WriteLine("Warning: header line has duplicate column names; ignoring duplicate " + token[i]);
                                continue;
                            }
                            _header.Add(i, token[i]);
                            _data[token[i]] = new List<string>();
                        }
                        headerParsed = true;
                        continue;
                    }

                    if (token.Length > _header.Count)
                        continue;

                    for (var i = 0; i < token.Length; i++)
                    {
                        if (!_header.ContainsKey(i))
                            continue;

                        _data[_header[i]].Add(token[i]);
                    }
                    _rows.Add(line);
                }
            }
        }

    }
}
