using System.Linq;

namespace InformedProteomics.Backend.Data.Biology
{
    public class Enzyme
    {
        private Enzyme(string name, string residues, bool isNTerm, string description, string psiCvAccession)
        {
            _name = name;
            _residues = residues.ToCharArray();
            _isNTerm = isNTerm;
            _description = description;
            _psiCvAccession = psiCvAccession;
        }

        private readonly string _name;
        private readonly char[] _residues;
        private readonly bool _isNTerm;
        private readonly string _description;
        private readonly string _psiCvAccession;

        public string Name
        {
            get { return _name; }
        }

        public char[] Residues
        {
            get { return _residues; }
        }

        public bool IsNTerm
        {
            get { return _isNTerm; }
        }

        public string Description
        {
            get { return _description; }
        }

        public string PsiCvAccession
        {
            get { return _psiCvAccession; }
        }

        public bool IsCleavable(char residue)
        {
            if (_residues.Length == 0)
                return true;
            return _residues.Any(r => r == residue);
        }

        public static readonly Enzyme UnspecificCleavage = new Enzyme("UnspecificCleavage", "", false, "unspecific cleavage", "MS:1001956");
        public static readonly Enzyme Trypsin = new Enzyme("Trypsin", "KR", false, "Trypsin", "MS:1001251");
        public static readonly Enzyme Chymotrypsin = new Enzyme("Chymotrypsin", "FYWL", false, "Chymotrypsin", "MS:1001306");
        public static readonly Enzyme LysC = new Enzyme("LysC", "K", false, "Lys-C", "MS:1001309");
        public static readonly Enzyme LysN = new Enzyme("LysN", "K", true, "Lys-N", null);
        public static readonly Enzyme GluC = new Enzyme("GluC","E",false, "glutamyl endopeptidase", "MS:1001917");
        public static readonly Enzyme ArgC = new Enzyme("ArgC","R",false, "Arg-C", "MS:1001303");
        public static readonly Enzyme AspN = new Enzyme("AspN","D",true, "Asp-N", "MS:1001304");
        public static readonly Enzyme Alp = new Enzyme("aLP", "", false, "alphaLP", null);
        public static readonly Enzyme NoCleavage = new Enzyme("NoCleavage", "", false, "no cleavage", "MS:1001955");
    }
}
