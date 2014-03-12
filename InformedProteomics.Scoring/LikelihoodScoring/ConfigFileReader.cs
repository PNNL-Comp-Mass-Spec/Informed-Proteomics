using System;
using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.Scoring.LikelihoodScoring
{
   class Node
    {
        public string Header { get; private set; }      // The header tag

        // All of the key/value pairs:
        public Dictionary<String, String> Contents { get; private set; }

        // constructor
        public Node(string header, Dictionary<String, String> contents)
        {
            Header = header;
            Contents = contents;
        }
    }
    class ConfigFileReader
    {
        public List<Node> Nodes { get; protected set; }

        // Exception
        public class InvalidHeader : Exception { }

        private bool ValidHeader(string header)
        {
            return (header[0] == '[' && header[header.Length - 1] == ']');
        }

        /* read() Read the file and store the result
         * in Nodes.
         */
        private void Read(string fileName)
        {
            Node currentNode = null;
            var keyvalue = new Dictionary<String, String>();
            string[] lines = System.IO.File.ReadAllLines(fileName);
            char[] headerbrackets = { '[', ']' };
            string header = "";
            foreach (var line in lines)
            {
                string commentsStripped = line.Split('#')[0];      // remove comments
                string[] parts = commentsStripped.Split('=');       // split key/value
                if (parts.Length < 2)
                {
                    // The line is either a header, empty line, or invalid
                    parts[0] = parts[0].Trim().ToLower();
                    if (parts[0] == "")
                        // empty line
                        continue;
                    if (currentNode == null)
                    {
                        // first node in the file
                        currentNode = new Node(null, null);
                        header = parts[0].Trim(headerbrackets);
                    }
                    else
                    {
                        // this isn't the first node in the file
                        currentNode = new Node(header, keyvalue);
                        keyvalue = new Dictionary<String, String>();
                        Nodes.Add(currentNode);
                        header = parts[0].Trim(headerbrackets);
                    }

                    if (!ValidHeader(parts[0]))
                        // invalid header
                        throw new InvalidHeader();
                }
                else
                {
                    // key value pair
                    string key = parts[0].Trim().ToLower();
                    string value = parts[1].Trim();
                    keyvalue.Add(key, value);
                }
            }
            currentNode = new Node(header, keyvalue);
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
