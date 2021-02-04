using System;
using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Scoring.LikelihoodScoring.Config
{
    public class Node
    {
        public string Header { get; }      // The header tag

        // All of the key/value pairs:
        public Dictionary<string, string> Contents { get; }

        // constructor
        public Node(string header, Dictionary<string, string> contents)
        {
            Header = header;
            Contents = contents;
        }
    }
    public class ConfigFileReader
    {
        public List<Node> Nodes { get; protected set; }

        // Exception
        public class InvalidHeaderException : Exception
        {
            public InvalidHeaderException() : base()
            {
            }

            public InvalidHeaderException(string message) : base(message)
            {
            }

            public InvalidHeaderException(string message, Exception innerException) : base(message, innerException)
            {
            }
        }

        private bool ValidHeader(string header)
        {
            if (string.IsNullOrWhiteSpace(header) || header.Length < 2)
                return false;

            return header[0] == '[' && header[header.Length - 1] == ']';
        }

        /* read() Read the file and store the result
         * in Nodes.
         */
        private void Read(string fileName)
        {
            Node currentNode = null;
            var keyValue = new Dictionary<string, string>();
            var lines = System.IO.File.ReadAllLines(fileName);
            char[] headerBrackets = { '[', ']' };
            var headerName = string.Empty;

            foreach (var line in lines)
            {
                var commentsStripped = line.Split('#')[0];      // remove comments
                var parts = commentsStripped.Split('=');       // split key/value

                if (parts.Length < 1)
                    continue;

                if (parts.Length < 2)
                {
                    // The line is either a header, an empty line, or invalid
                    parts[0] = parts[0].Trim().ToLower();
                    if (string.IsNullOrWhiteSpace(parts[0]))
                    {
                        // empty line
                        continue;
                    }

                    if (currentNode == null)
                    {
                        // First node in the file
                        currentNode = new Node(null, null);
                        headerName = parts[0].Trim(headerBrackets);
                    }
                    else
                    {
                        // This isn't the first node in the file
                        // Store the values for the current header
                        currentNode = new Node(headerName, keyValue);
                        Nodes.Add(currentNode);

                        // Initialize the new header
                        keyValue = new Dictionary<string, string>();
                        headerName = parts[0].Trim(headerBrackets);
                    }

                    if (!ValidHeader(parts[0]))
                    {
                        // Invalid header; should be in the form
                        // [HeaderName]
                        throw new InvalidHeaderException();
                    }
                    continue;
                }


                // key value pair
                var key = parts[0].Trim().ToLower();
                var value = parts[1].Trim();
                keyValue.Add(key, value);
            }

            if (string.IsNullOrEmpty(headerName) || keyValue.Count == 0)
                return;

            // Store the values for the current header
            currentNode = new Node(headerName, keyValue);
            Nodes.Add(currentNode);
        }

        /*
         *  Constructor
         */
        public ConfigFileReader(string fileName)
        {
            Nodes = new List<Node>();
            Read(fileName);
        }

        /*
         * getNodes() return a list of all the nodes with a particular
         * header tag.
         */
        public List<Node> GetNodes(string headerTag)
        {
            return (from i in Nodes
                    where i.Header == headerTag
                    select i).ToList();
        }
    }
}
