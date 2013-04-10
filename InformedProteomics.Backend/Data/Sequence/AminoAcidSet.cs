using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Enum;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class AminoAcidSet
    {
        #region Constructors

        // Generate an amino acid map with 20 standard amino acids
        public AminoAcidSet()
        {
            _residueMap.Add(AminoAcid.ProteinNTerm.Residue, AminoAcid.ProteinNTerm);
            _residueMap.Add(AminoAcid.ProteinCTerm.Residue, AminoAcid.ProteinCTerm);
            _residueMap.Add(AminoAcid.PeptideNTerm.Residue, AminoAcid.PeptideNTerm);
            _residueMap.Add(AminoAcid.PeptideCTerm.Residue, AminoAcid.PeptideCTerm);

            foreach (var aa in AminoAcid.StandardAminoAcidArr)
			{
				_residueMap.Add(aa.Residue, aa);
			}

        	_modificationParams = new ModificationParams();
            _residueModMap = new Dictionary<char, int[]>();
        }

        // Generate an amino acid map with Cys static modification
        public AminoAcidSet(Modification cysMod): this()
        {
            var cys = AminoAcid.GetStandardAminoAcid('C');
            var modCys = new ModifiedAminoAcid(cys, cysMod);
            _residueMap['C'] = modCys;
        }

        public AminoAcidSet(string modFilePath): this(new ModFileParser(modFilePath))
        {
        }

        private AminoAcidSet(ModFileParser modFileParser)
            : this(modFileParser.SearchModifications, modFileParser.MaxNumDynModsPerSequence)
        {
            
        }

        public AminoAcidSet(IEnumerable<SearchModification> searchModifications, int maxNumModsPerSequence): this()
        {
            var modifications = searchModifications as SearchModification[] ?? searchModifications.ToArray();
            // apply fixed modifications
            foreach (var searchModification in modifications)
            {
                if (searchModification.IsFixedModification)
                {
                    var targetResidue = searchModification.TargetResidue;
                    _residueMap[targetResidue] = new ModifiedAminoAcid(_residueMap[targetResidue], searchModification.Modification);
                }
            }

            // apply dynamic modifications
            var dynamicModifications = new HashSet<Modification>();
            var residueModMap = new Dictionary<char, List<Modification>>();
            foreach (var searchModification in modifications)
            {
                if (!searchModification.IsFixedModification)
                {
                    dynamicModifications.Add(searchModification.Modification);
                    var targetResidue = searchModification.TargetResidue;
                    List<Modification> modList;
                    if (residueModMap.TryGetValue(targetResidue, out modList))
                    {
                        modList.Add(searchModification.Modification);
                    }
                    else
                    {
                        residueModMap.Add(targetResidue, new List<Modification> {searchModification.Modification});
                    }
                }
            }

            var dynModArray = dynamicModifications.ToArray();
            _modificationParams = new ModificationParams(dynModArray, maxNumModsPerSequence);

            // initialize resideModMap
            _residueModMap = new Dictionary<char, int[]>();

            foreach (var entry in residueModMap)
            {
                var residue = entry.Key;
                var modList = entry.Value;

                var modIndexArr = new int[modList.Count];
                int index = -1;
                foreach (var mod in modList)
                {
                    int modIndex = Array.IndexOf(dynModArray, mod);
                    modIndexArr[++index] = modIndex;
                }
                _residueModMap.Add(residue, modIndexArr);
            }
        }

        #endregion

        #region Properties

        #endregion
        
        public AminoAcid GetAminoAcid(char residue)
        {
            AminoAcid aminoAcid;
            if (_residueMap.TryGetValue(residue, out aminoAcid))
            {
                return aminoAcid;
            }
            return null;
        }

        public AminoAcid GetAminoAcid(char residue, SequenceLocation location)
        {
            var residueMap = _locationSpecificResidueMap[location];
            AminoAcid aminoAcid;
            if (residueMap.TryGetValue(residue, out aminoAcid))
            {
                return aminoAcid;
            }
            return null;
        }

        public int[] GetModificationIndices(char residue)
        {
            int[] modIndexArr;
            if (_residueModMap.TryGetValue(residue, out modIndexArr))
            {
                return modIndexArr;
            }
            return new int[0];
        }

        public int[] GetModificationIndices(char residue, SequenceLocation location)
        {
            throw new System.NotImplementedException();
        }

        public Composition GetComposition(string sequence)
        {
            Composition composition = Composition.Zero;

            foreach (char residue in sequence)
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

        public ModificationParams GetModificationParams()
        {
            return _modificationParams;
        }

        #region Private Members

        private readonly Dictionary<char, AminoAcid> _residueMap = new Dictionary<char, AminoAcid>();
        private readonly Dictionary<SequenceLocation, Dictionary<char, AminoAcid>> _locationSpecificResidueMap = new Dictionary<SequenceLocation, Dictionary<char, AminoAcid>>();
        private readonly Dictionary<char, int[]> _residueModMap;

        private readonly ModificationParams _modificationParams;

        #endregion
    }

    public class ModifiedAminoAcid: AminoAcid
    {

        public ModifiedAminoAcid(AminoAcid aa, Modification modification) 
            : base(aa.Residue, aa.Name+"+"+modification.Name, aa.Composition+modification.Composition)
        {
            _modification = modification;
        }

        private readonly Modification _modification;

        public Modification Modification
        {
            get { return _modification; }
        }
    }
}
