using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Enum;

namespace InformedProteomics.Backend.Data.Sequence
{
    /// <summary>
    /// Set of amino acids used in a search
    /// </summary>
    public class AminoAcidSet
    {
        #region Constructors

        /// <summary>
        /// Generate an amino acid map with 20 standard amino acids
        /// </summary>
        public AminoAcidSet()
        {
            // Initialization
            _modificationParams = new ModificationParams();
            SearchModifications = new List<SearchModification>();

            _locationSpecificResidueMap = new Dictionary<SequenceLocation, Dictionary<char, AminoAcid>>();
            _locationSpecificResidueModMap = new Dictionary<SequenceLocation, Dictionary<char, int[]>>();
            foreach (var loc in AllSequenceLocations)
            {
                _locationSpecificResidueModMap[loc] = new Dictionary<char, int[]>();
            }

            // Initialize standard amino acid set
            foreach (var loc in AllSequenceLocations)
            {
                var residueMap = new Dictionary<char, AminoAcid>();
                switch (loc)
                {
                    case SequenceLocation.PeptideNTerm:
                        residueMap.Add(AminoAcid.PeptideNTerm.Residue, AminoAcid.PeptideNTerm);
                        break;
                    case SequenceLocation.PeptideCTerm:
                        residueMap.Add(AminoAcid.PeptideCTerm.Residue, AminoAcid.PeptideCTerm);
                        break;
                    case SequenceLocation.ProteinNTerm:
                        residueMap.Add(AminoAcid.ProteinNTerm.Residue, AminoAcid.ProteinNTerm);
                        break;
                    case SequenceLocation.ProteinCTerm:
                        residueMap.Add(AminoAcid.ProteinCTerm.Residue, AminoAcid.ProteinCTerm);
                        break;
                }

                // loc is SequenceLocation.Everywhere
                foreach (var aa in AminoAcid.StandardAminoAcidArr)
                {
                    residueMap.Add(aa.Residue, aa);
                }
                _locationSpecificResidueMap[loc] = residueMap;
            }
        }

        /// <summary>
        /// Generate an amino acid map with a Cysteine static modification
        /// </summary>
        /// <param name="cysMod"></param>
        public AminoAcidSet(Modification cysMod) : this()
        {
            foreach (var loc in AllSequenceLocations)
            {
                var residueMap = _locationSpecificResidueMap[loc];
                residueMap['C'] = new ModifiedAminoAcid(residueMap['C'], cysMod);
            }
        }

        /// <summary>
        /// Construct an amino acid map using the 20 standard amino acids and any modifications in mod file at path <paramref name="modFilePath"/>
        /// </summary>
        /// <param name="modFilePath"></param>
        public AminoAcidSet(string modFilePath) : this(new ModFileParser(modFilePath))
        {
        }

        /// <summary>
        /// Construct an amino acid map using the 20 standard amino acids and any modifications in mod file parsed by <paramref name="modFileParser"/>
        /// </summary>
        /// <param name="modFileParser"></param>
        private AminoAcidSet(ModFileParser modFileParser)
            : this(modFileParser.SearchModifications, modFileParser.MaxNumDynModsPerSequence)
        {
        }

        /// <summary>
        /// Construct an amino acid map using the 20 standard amino acids and any modifications in <paramref name="searchModifications"/>
        /// </summary>
        /// <param name="searchModifications"></param>
        /// <param name="maxNumModsPerSequence"></param>
        public AminoAcidSet(IEnumerable<SearchModification> searchModifications, int maxNumModsPerSequence) : this()
        {
            if (searchModifications == null) return;

            var modifications = searchModifications as SearchModification[] ?? searchModifications.ToArray();
            SearchModifications = new List<SearchModification>(modifications);

            // apply fixed modifications
            foreach (var searchModification in modifications)
            {
                if (!searchModification.IsFixedModification) continue;

                var location = searchModification.Location;
                var targetResidue = searchModification.TargetResidue;

                foreach (var loc in AffectedLocations[location])
                {
                    var appliedResidue = targetResidue != '*' ? targetResidue : SequenceLocationToLocationResidue[loc];
                    var residueMap = _locationSpecificResidueMap[loc];
                    residueMap[appliedResidue] = new ModifiedAminoAcid(residueMap[appliedResidue], searchModification.Modification);
                }
            }

            if (maxNumModsPerSequence <= 0) return;

            // apply dynamic modifications
            var dynamicModifications = new HashSet<Modification>();
            var locationSpecificResidueVariableModMap = new Dictionary<SequenceLocation, Dictionary<char, List<Modification>>>();
            foreach (var loc in AllSequenceLocations)
            {
                locationSpecificResidueVariableModMap[loc] = new Dictionary<char, List<Modification>>();
            }

            foreach (var searchModification in modifications)
            {
                if (searchModification.IsFixedModification) continue;

                dynamicModifications.Add(searchModification.Modification);
                var location = searchModification.Location;
                var targetResidue = searchModification.TargetResidue;

                foreach (var loc in AffectedLocations[location])
                {
                    var residueModMap = locationSpecificResidueVariableModMap[loc];
                    var appliedResidue = targetResidue != '*' ? targetResidue : SequenceLocationToLocationResidue[loc];

                    if (residueModMap.TryGetValue(appliedResidue, out var modList))
                    {
                        modList.Add(searchModification.Modification);
                    }
                    else
                    {
                        residueModMap.Add(appliedResidue, new List<Modification> { searchModification.Modification });
                    }
                }
            }

            var dynModArray = dynamicModifications.ToArray();
            _modificationParams = new ModificationParams(dynModArray, maxNumModsPerSequence);

            foreach (var loc in AllSequenceLocations)
            {
                var residueModMap = locationSpecificResidueVariableModMap[loc];
                foreach (var entry in residueModMap)
                {
                    var residue = entry.Key;
                    var modList = entry.Value;

                    var modIndexArr = new int[modList.Count];
                    var index = -1;
                    foreach (var mod in modList)
                    {
                        var modIndex = Array.IndexOf(dynModArray, mod);
                        modIndexArr[++index] = modIndex;
                    }
                    _locationSpecificResidueModMap[loc].Add(residue, modIndexArr);
                }
            }
        }

        #endregion

        #region Properties

        /// <summary>
        /// The list of search modifications provided via the mods file.
        /// </summary>
        public IEnumerable<SearchModification> SearchModifications { get; }

        #endregion

        /// <summary>
        /// Get the amino acid specified by character <paramref name="residue"/>, with any static modifications
        /// </summary>
        /// <param name="residue"></param>
        /// <returns></returns>
        public AminoAcid GetAminoAcid(char residue)
        {
            return GetAminoAcid(residue, SequenceLocation.Everywhere);
        }
        /// <summary>
        /// Get the amino acid specified by character <paramref name="residue"/> that occurs at location <paramref name="location"/>, with any static modifications
        /// </summary>
        /// <param name="residue"></param>
        /// <param name="location"></param>
        /// <returns></returns>
        public AminoAcid GetAminoAcid(char residue, SequenceLocation location)
        {
            var residueMap = _locationSpecificResidueMap[location];

            if (residueMap.TryGetValue(residue, out var aminoAcid))
            {
                return aminoAcid;
            }
            return null;
        }

        /// <summary>
        /// Get the indices of modifications that apply to residue <paramref name="residue"/>
        /// </summary>
        /// <param name="residue"></param>
        /// <returns></returns>
        public int[] GetModificationIndices(char residue)
        {
            return GetModificationIndices(residue, SequenceLocation.Everywhere);
        }

        /// <summary>
        /// Get the indices of modifications that apply to residue <paramref name="residue"/> at location <paramref name="location"/>
        /// </summary>
        /// <param name="residue"></param>
        /// <param name="location"></param>
        /// <returns></returns>
        public int[] GetModificationIndices(char residue, SequenceLocation location)
        {
            var residueModMap = _locationSpecificResidueModMap[location];

            if (residueModMap.TryGetValue(residue, out var modIndexArr))
            {
                return modIndexArr;
            }
            return new int[0];
        }

        /// <summary>
        /// Get the composition of <paramref name="sequence"/> using this amino acid map (with static modification masses added)
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        /// <remarks>This method ignores N- and C-term specific modifications.</remarks>
        public Composition.Composition GetComposition(string sequence)
        {
            var composition = Composition.Composition.Zero;

            foreach (var residue in sequence)
            {
                var aaArr = GetAminoAcid(residue);
                if (aaArr == null)
                {
                    return null;
                }
                composition += aaArr.Composition;
            }

            return composition;
        }

        //public Composition.Composition GetComposition(AminoAcid nTerm, string sequence, AminoAcid cTerm)

        /// <summary>
        /// Get the parameters governing use of dynamic modifications
        /// </summary>
        /// <returns></returns>
        public ModificationParams GetModificationParams()
        {
            return _modificationParams;
        }

        /// <summary>
        /// Display this amino acid set on the console
        /// </summary>
        public void Display()
        {
            foreach (var location in AllSequenceLocations)
            {
                Console.WriteLine("===== " + location + " =====");
                var keys = _locationSpecificResidueMap[location].Keys.ToArray();
                Array.Sort(keys);
                foreach (var residue in keys)
                {
                    var aa = GetAminoAcid(residue, location);
                    Console.Write("{0}\t{1}\t{2}", residue, aa.Mass, aa.Composition);
                    foreach (var modIndex in GetModificationIndices(residue, location))
                    {
                        Console.Write("\t" + _modificationParams.GetModification(modIndex));
                    }
                    Console.WriteLine();
                }
            }
        }

        /// <summary>
        /// Get the set of unique compositions within this amino acid set
        /// </summary>
        /// <returns></returns>
        public Composition.Composition[] GetUniqueCompositions()
        {
            var uniqueCompositionList = new List<Composition.Composition>();
            var tempModTable = new List<Composition.Composition>[AminoAcid.StandardAminoAcidArr.Length];
            for (var i = 0; i < AminoAcid.StandardAminoAcidArr.Length; i++)
            {
                tempModTable[i] = new List<Composition.Composition>();
                uniqueCompositionList.Add(AminoAcid.StandardAminoAcidArr[i].Composition);
            }

            foreach (var location in AllSequenceLocations)
            {
                //var keys = _locationSpecificResidueMap[location].Keys.ToArray();
                for (var i = 0; i < AminoAcid.StandardAminoAcidArr.Length; i++)
                {
                    var residue = AminoAcid.StandardAminoAcidArr[i].Residue;
                    var aa = GetAminoAcid(residue, location);

                    foreach (var modIndex in GetModificationIndices(residue, location))
                    {
                        var mAa = new ModifiedAminoAcid(aa, _modificationParams.GetModification(modIndex));
                        if (tempModTable[i].Contains(mAa.Composition)) continue;
                        tempModTable[i].Add(mAa.Composition);
                        uniqueCompositionList.Add(mAa.Composition);
                    }
                    //Console.WriteLine();
                }
            }

            return uniqueCompositionList.ToArray();
        }

        /// <summary>
        /// Get the standard amino acid set, with no modifications
        /// </summary>
        /// <returns></returns>
        public static AminoAcidSet GetStandardAminoAcidSet()
        {
            return _standardAminoAcidSet ?? (_standardAminoAcidSet = new AminoAcidSet());
        }

        /// <summary>
        /// Get the standard amino acid set, with static Cysteine Alkylation (Carboamidomethyl)
        /// </summary>
        /// <returns></returns>
        public static AminoAcidSet GetStandardAminoAcidSetWithCarboamidomethylCys()
        {
            return _standardAminoAcidSetWithCarboamidomethylCys ??
                (_standardAminoAcidSetWithCarboamidomethylCys = new AminoAcidSet(Modification.Carbamidomethylation));
        }

        private static AminoAcidSet _standardAminoAcidSet;
        private static AminoAcidSet _standardAminoAcidSetWithCarboamidomethylCys;

        #region Private Members

        private readonly Dictionary<SequenceLocation, Dictionary<char, AminoAcid>> _locationSpecificResidueMap;
        private readonly Dictionary<SequenceLocation, Dictionary<char, int[]>> _locationSpecificResidueModMap;

        private readonly ModificationParams _modificationParams;

        #endregion

        private static readonly IEnumerable<SequenceLocation> AllSequenceLocations;
        private static readonly Dictionary<SequenceLocation, char> SequenceLocationToLocationResidue;
        private static readonly Dictionary<SequenceLocation, IEnumerable<SequenceLocation>> AffectedLocations;

        static AminoAcidSet()
        {
            AllSequenceLocations = System.Enum.GetValues(typeof(SequenceLocation)) as SequenceLocation[];
            SequenceLocationToLocationResidue = new Dictionary<SequenceLocation, char>
            {
                {SequenceLocation.PeptideNTerm, AminoAcid.PeptideNTerm.Residue},
                {SequenceLocation.PeptideCTerm, AminoAcid.PeptideCTerm.Residue},
                {SequenceLocation.ProteinNTerm, AminoAcid.ProteinNTerm.Residue},
                {SequenceLocation.ProteinCTerm, AminoAcid.ProteinCTerm.Residue}
            };

            AffectedLocations = new Dictionary<SequenceLocation, IEnumerable<SequenceLocation>>
            {
                {
                    SequenceLocation.Everywhere, new[]
                    {
                        SequenceLocation.Everywhere,
                        SequenceLocation.PeptideNTerm,
                        SequenceLocation.PeptideCTerm,
                        SequenceLocation.ProteinNTerm,
                        SequenceLocation.ProteinCTerm
                    }
                },
                {
                    SequenceLocation.PeptideNTerm, new[]
                    {
                        SequenceLocation.PeptideNTerm,
                        SequenceLocation.ProteinNTerm
                    }
                },
                {
                    SequenceLocation.PeptideCTerm, new[]
                    {
                        SequenceLocation.PeptideCTerm,
                        SequenceLocation.ProteinCTerm
                    }
                },
                {
                    SequenceLocation.ProteinNTerm, new[]
                    {
                        SequenceLocation.ProteinNTerm
                    }
                },
                {
                    SequenceLocation.ProteinCTerm, new[]
                    {
                        SequenceLocation.ProteinCTerm
                    }
                }
            };

        }
    }
}
