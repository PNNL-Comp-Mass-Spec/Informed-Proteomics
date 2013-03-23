using System;
using System.Collections.Generic;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class AminoAcidSet
    {
        #region Constructors

        // Generate an amino acid map with 20 standard amino acids
        public AminoAcidSet()
        {
			foreach (var aa in AminoAcid.StandardAminoAcidArr)
			{
				_residueMap.Add(aa.Residue, new[] {aa});
			}

        	_numDynMods = 0;
        }

        // Generate an amino acid map with Cys static modification
        public AminoAcidSet(Modification cysMod)
        {
            foreach (var aa in AminoAcid.StandardAminoAcidArr)
            {
                AminoAcid cys = AminoAcid.GetStandardAminoAcid('C');
				if (aa != cys)
				{
					_residueMap.Add(aa.Residue, new[] {aa});
				}
				else
				{
					var modCys = new AminoAcid('C', cysMod.Name + " C", cys.Composition + cysMod.Composition);
					_residueMap.Add('C', new[] { modCys });
				}
            }
            _numDynMods = 0;
        }

        // TODO: Read modification information from modFileName
        public AminoAcidSet(string modFileName)
        {
            _numDynMods = 0;
        }


        #endregion

        #region Properties

        private int _maxNumDynModsPerPeptide = 2;

        public int MaxNumDynModsPerPeptide
        {
            get { return _maxNumDynModsPerPeptide; }
            set { _maxNumDynModsPerPeptide = value; }
        }

        #endregion

        public AminoAcid GetUnmodifiedAminoAcid(char residue)
        {
            AminoAcid[] aaArr = GetAminoAcids(residue);
			return aaArr == null ? null : aaArr[0];
        }
            
        public AminoAcid[] GetAminoAcids(char residue)
        {
            return _residueMap[residue];
        }

        public IEnumerable<Composition> GetCompositions(string sequence)
        {
            if (_numDynMods == 0)
            {
                Composition composition = Composition.Zero;

                foreach (char residue in sequence)
                {
                    var aaArr = GetAminoAcids(residue);
					if (aaArr == null)
					{
						return null;
					}
                	composition += aaArr[0].Composition;
                }

                return new[] {composition};
            }

            // TODO: not yet implemented
            throw new NotImplementedException("Dynamic mods have not been implemented.");
        }

        //public AminoAcid[] GetAminoAcids(AminoAcid aa, SequenceLocation location)
        //{
        //    throw new System.NotImplementedException();
        //}

        #region Private Members

        private readonly Dictionary<char, AminoAcid[]> _residueMap = new Dictionary<char, AminoAcid[]>();
        private readonly int _numDynMods;

        #endregion
    }
}
