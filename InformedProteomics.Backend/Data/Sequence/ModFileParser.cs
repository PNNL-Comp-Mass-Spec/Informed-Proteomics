using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class ModFileParser
    {
        public ModFileParser(string modFilePath)
        {
            _modFilePath = modFilePath;
            _searchModifications = Parse(modFilePath, out _maxNumDynModsPerSequence);
        }

        public string ModFilePath
        {
            get { return _modFilePath;  }
        }

        public IEnumerable<SearchModification> SearchModifications
        {
            get { return _searchModifications; }
        }

        public int MaxNumDynModsPerSequence
        {
            get { return _maxNumDynModsPerSequence;  }
        }

        private readonly string _modFilePath;
        private readonly IEnumerable<SearchModification> _searchModifications;
        private readonly int _maxNumDynModsPerSequence;

        private IEnumerable<SearchModification> Parse(string modFilePath, out int maxNumDynModsPerPeptide)
        {
            var searchModList = new List<SearchModification>();

            TextReader textReader = new StreamReader(modFilePath);

            var buf = "";
            while ((buf = textReader.ReadLine()) != null)
            {
                // Comment 
                if (buf.Length == 0 || buf.StartsWith("#")) continue;
                if(buf.StartsWith("NumMods"))
                {

                }
                //
                var token = buf.Split(null);
                // TODO: Parse Search Modifications
                // TODO: set up maxNumDynModsPerPeptide
                //var searchModification = new SearchModification(mod, targetResidue, location, isFixedModification);
            }

            textReader.Close();

            //return searchModList.ToArray();
            maxNumDynModsPerPeptide = 0;

            // TODO: * cannot come with any
            return null;
        }
    }
}
