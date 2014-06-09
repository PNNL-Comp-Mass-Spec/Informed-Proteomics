using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Enum;

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

        private static IEnumerable<SearchModification> Parse(string modFilePath, out int maxNumDynModsPerPeptide)
        {
            var searchModList = new List<SearchModification>();

            var numMods = 0;
            var lineNum = 0;
            foreach(var line in File.ReadLines(modFilePath))
		    {
			    lineNum++;
                var tokenArr = line.Split('#');
			    if(tokenArr.Length == 0) continue;
			    var s = tokenArr[0].Trim();
			    if(s.Length == 0) continue;

			    if(s.StartsWith("NumMods="))
			    {
				    try {
                        numMods = Convert.ToInt32(s.Split('=')[1].Trim());
				    } catch (FormatException)
				    {
					    Console.WriteLine("{0}: Illegal NumMods parameter at line {1} - {2}", modFilePath, lineNum, s);
				        maxNumDynModsPerPeptide = -1;
				        return null;
				    }
			    }
			    else
			    {
				    var token = s.Split(',');
				    if(token.Length != 5) continue;

				    // Composition
				    
				    var compStr = token[0].Trim();
			        var composition = Composition.Composition.ParseFromPlainString(compStr) ??
			                          Composition.Composition.Parse(compStr);
			        if (composition == null)
			        {
                        Console.WriteLine("{0}: Illegal Composition at line {1} - {2}", modFilePath, lineNum, s);
                        maxNumDynModsPerPeptide = -1;
                        return null;
			        }
				
				    // Residues
				    var residueStr = token[1].Trim();
				    var isResidueStrLegitimate = residueStr.Equals("*") 
                        || residueStr.Any() && residueStr.All(AminoAcid.IsStandardAminoAcidResidue);
                    if(!isResidueStrLegitimate)
				    {
                        Console.WriteLine("{0}: Illegal residues at line {1} - {2}", modFilePath, lineNum, s);
                        maxNumDynModsPerPeptide = -1;
                        return null;
                    }
				
				    // isFixedModification
				    bool isFixedModification;
				    if(token[2].Trim().Equals("fix", StringComparison.InvariantCultureIgnoreCase)) isFixedModification = true;
				    else if(token[2].Trim().Equals("opt", StringComparison.InvariantCultureIgnoreCase)) isFixedModification = false;
				    else
				    {
                        Console.WriteLine("{0}: Illegal modification type (fix or opt) at line {1} - {2}", modFilePath, lineNum, s);
                        maxNumDynModsPerPeptide = -1;
                        return null;
                    }
					
				    // Location
				    SequenceLocation location;
				    var locStr = token[3].Split()[0].Trim();
			        if (locStr.Equals("any", StringComparison.InvariantCultureIgnoreCase))
			            location = SequenceLocation.Everywhere;
                    else if (locStr.Equals("N-Term", StringComparison.InvariantCultureIgnoreCase) ||
                             locStr.Equals("NTerm", StringComparison.InvariantCultureIgnoreCase))
                        location = SequenceLocation.PeptideNTerm;
                    else if (locStr.Equals("C-Term", StringComparison.InvariantCultureIgnoreCase) ||
                             locStr.Equals("CTerm", StringComparison.InvariantCultureIgnoreCase))
                        location = SequenceLocation.PeptideCTerm;
                    else if (locStr.Equals("Prot-N-Term", StringComparison.InvariantCultureIgnoreCase) ||
                             locStr.Equals("ProtNTerm", StringComparison.InvariantCultureIgnoreCase))
                        location = SequenceLocation.ProteinNTerm;
                    else if (locStr.Equals("Prot-C-Term", StringComparison.InvariantCultureIgnoreCase) ||
                             locStr.Equals("ProtCTerm", StringComparison.InvariantCultureIgnoreCase))
                        location = SequenceLocation.ProteinCTerm;
				    else
				    {
                        Console.WriteLine("{0}: Illegal modification location (fix or opt) at line {1} - {2}", modFilePath, lineNum, s);
                        maxNumDynModsPerPeptide = -1;
                        return null;
                    }

                    // Check if it's valid
			        if (residueStr.Equals("*") && location == SequenceLocation.Everywhere)
			        {
                        Console.WriteLine("{0}: Invalid modification: * should not be applied to \"any\"");
                        maxNumDynModsPerPeptide = -1;
                        return null;
			        }

				    var name = token[4].Split()[0].Trim();

			        var mod = Modification.Get(name) ?? Modification.RegisterAndGetModification(name, composition);

			        searchModList.AddRange(residueStr.Select(
                        residue => new SearchModification(mod, residue, location, isFixedModification)
                        ));
			    }
		    }
            maxNumDynModsPerPeptide = numMods;
            return searchModList;
        }
    }
}
