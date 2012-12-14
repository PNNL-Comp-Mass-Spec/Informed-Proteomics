namespace InformedProteomics.Backend.Utils
{
    public class ModifiedAminoAcid : AminoAcid
    {
        public ModifiedAminoAcid(char residue, string name, Composition comp) : base(residue, name, comp)
        {
        }

        public int Modification
        {
            get
            {
                throw new System.NotImplementedException();
            }
            set
            {
            }
        }
    
        public Modification GetModification()
        {
            throw new System.NotImplementedException();
        }
    }
}
