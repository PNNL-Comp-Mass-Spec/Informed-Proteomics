using System;
using System.Collections.Generic;
using System.IO;

namespace ProMex
{
    public class TsvFileReader
    {
        public TsvFileReader(string fileName, char delimiter = '\t')
        {
            FileName = fileName;
            _delimeter = delimiter;

            _fileReader = new StreamReader(fileName);
            var line = _fileReader.ReadLine();
            _header = line.Split(_delimeter);
        }

        public void Close()
        {
            _fileReader.Close();
        }
        public string FileName { get; private set; }

        public IList<string> GetHeaders()
        {
            return _header;
        }

        public int GetColumnIndex(string header)
        {
            return Array.FindIndex(_header, x => x.Equals(header));
        }

        private readonly string[] _header;
        private readonly char _delimeter;
        private readonly StreamReader _fileReader;

        public int LineCounter { get; private set; }

        public IEnumerable<string[]> ReadLine()
        {
            // parse header
            //var firstRow = true;
            string line;
            while((line = _fileReader.ReadLine()) != null)
            {
                var token = line.Split(_delimeter);
                LineCounter++;
                if (token.Length == _header.Length) yield return token;
            }
        }
    }
}
