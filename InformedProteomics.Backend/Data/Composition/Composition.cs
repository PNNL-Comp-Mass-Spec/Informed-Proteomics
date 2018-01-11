using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Composition
{
    /// <summary>
    /// Composition, consisting of elements C, H, N, O, S, or P, and optionally additional elements
    /// </summary>
    public class Composition: AbstractComposition
    {
        /// <summary>
        /// Empty composition
        /// </summary>
        public static readonly Composition Zero = new Composition(0, 0, 0, 0, 0);

        /// <summary>
        /// Water
        /// </summary>
        public static readonly Composition H2O = new Composition(0, 2, 0, 1, 0);

        /// <summary>
        /// Ammonia
        /// </summary>
        public static readonly Composition NH3 = new Composition(0, 3, 1, 0, 0);

        /// <summary>
        /// Amide
        /// </summary>
        public static readonly Composition NH2 = new Composition(0, 2, 1, 0, 0);

        /// <summary>
        /// Hydroxide
        /// </summary>
        public static readonly Composition OH = new Composition(0, 1, 0, 1, 0);

        /// <summary>
        /// Carbon Monoxide
        /// </summary>
        public static readonly Composition CO = new Composition(1, 0, 0, 1, 0);

        /// <summary>
        /// Hydrogen
        /// </summary>
        public static readonly Composition Hydrogen = new Composition(0, 1, 0, 0, 0);

        #region Constructors

        /// <summary>
        /// Constructor, with provided counts for C, H, N, O, S, P
        /// </summary>
        /// <param name="c"></param>
        /// <param name="h"></param>
        /// <param name="n"></param>
        /// <param name="o"></param>
        /// <param name="s"></param>
        /// <param name="p"></param>
        public Composition(int c, int h, int n, int o, int s, int p)
        {
            _c = (short)c;
            _h = (short)h;
            _n = (short)n;
            _o = (short)o;
            _s = (short)s;
            _p = (short)p;
        }

        /// <summary>
        /// Constructor, with provided counts for C, H, N, O, S
        /// </summary>
        /// <param name="c"></param>
        /// <param name="h"></param>
        /// <param name="n"></param>
        /// <param name="o"></param>
        /// <param name="s"></param>
        public Composition(int c, int h, int n, int o, int s)
            : this(c, h, n, o, s, 0)
        {
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        /// <param name="composition"></param>
        public Composition(Composition composition)
            : this(composition.C, composition.H,
                composition.N, composition.O, composition.S, composition.P)
        {
            if (composition._additionalElements != null)
            {
                _additionalElements = new Dictionary<Atom, short>(composition._additionalElements);
            }
        }

        /// <summary>
        /// Constructor, with provided counts for C, H, N, O, S, P, and a tuple of additional element and count
        /// </summary>
        /// <param name="c"></param>
        /// <param name="h"></param>
        /// <param name="n"></param>
        /// <param name="o"></param>
        /// <param name="s"></param>
        /// <param name="p"></param>
        /// <param name="additionalElement"></param>
        public Composition(int c, int h, int n, int o, int s, int p, Tuple<Atom, short> additionalElement)
            : this(c, h, n, o, s, p)
        {
            _additionalElements = new Dictionary<Atom, short> { { additionalElement.Item1, additionalElement.Item2 } };
        }

        /// <summary>
        /// Constructor, with provided counts for C, H, N, O, S, P, and tuples of additional elements and respective counts
        /// </summary>
        /// <param name="c"></param>
        /// <param name="h"></param>
        /// <param name="n"></param>
        /// <param name="o"></param>
        /// <param name="s"></param>
        /// <param name="p"></param>
        /// <param name="additionalElements"></param>
        public Composition(int c, int h, int n, int o, int s, int p, IEnumerable<Tuple<Atom, short>> additionalElements)
            : this(c, h, n, o, s, p)
        {
            _additionalElements = new Dictionary<Atom, short>();
            foreach (var element in additionalElements)
            {
                _additionalElements.Add(element.Item1, element.Item2);
            }
        }

        /// <summary>
        /// Empty constructor
        /// </summary>
        public Composition()
        {
            _additionalElements = new Dictionary<Atom, short>();
        }

        private Composition(int c, int h, int n, int o, int s, int p, Dictionary<Atom, short> additionalElements)
            : this(c, h, n, o, s, p)
        {
            _additionalElements = additionalElements;
        }
        #endregion

        #region Properties

        /// <summary>
        /// Count of Carbon atoms in the composition
        /// </summary>
        public short C
        {
            get => _c;
            set
            {
                _c = value;
                _mass = null;
                _nominalMass = null;
            }
        }

        /// <summary>
        /// Count of Hydrogen atoms in the composition
        /// </summary>
        public short H
        {
            get => _h;
            set
            {
                _h = value;
                _mass = null;
                _nominalMass = null;
            }
        }

        /// <summary>
        /// Count of Nitrogen atoms in the composition
        /// </summary>
        public short N
        {
            get => _n;
            set
            {
                _n = value;
                _mass = null;
                _nominalMass = null;
            }
        }

        /// <summary>
        /// Count of Oxygen atoms in the composition
        /// </summary>
        public short O
        {
            get => _o;
            set
            {
                _o = value;
                _mass = null;
                _nominalMass = null;
            }
        }

        /// <summary>
        /// Count of Sulfur atoms in the composition
        /// </summary>
        public short S
        {
            get => _s;
            set
            {
                _s = value;
                _mass = null;
                _nominalMass = null;
            }
        }

        /// <summary>
        /// Count of Phosphorus atoms in the composition
        /// </summary>
        public short P
        {
            get => _p;
            set
            {
                _p = value;
                _mass = null;
                _nominalMass = null;
            }
        }

        /// <summary>
        /// Map of additional elements in the composition and their respective counts
        /// </summary>
        public Dictionary<Atom, short> AdditionalElements
        {
            get => _additionalElements;
            set
            {
                _additionalElements = value;
                _mass = null;
                _nominalMass = null;
            }
        }

        /// <summary>
        /// Mass of the composition
        /// </summary>
        public override double Mass => (double)(_mass ?? (_mass = GetMonoIsotopicMass()));

        /// <summary>
        /// Nominal mass of the composition
        /// </summary>
        public override int NominalMass => (int)(_nominalMass ?? (_nominalMass = GetNominalMass()));

        #endregion

        #region Private members

        private Dictionary<Atom, short> _additionalElements;
        private short _c;
        private short _h;
        private short _n;
        private short _o;
        private short _s;
        private short _p;

        private double? _mass;
        private int? _nominalMass;

        #endregion

        /// <inheritdoc />
        public override int GetHashCode()
        {
            return Mass.GetHashCode();
        }

        /// <inheritdoc />
        public override bool Equals(object obj)
        {
            if (!(obj is Composition other))
                return false;

            if (_c != other.C || _h != other.H || _n != other.N || _o != other.O || _s != other.S || _p != other.P)
                return false;

            if (_additionalElements == null) return other._additionalElements == null;

            if (other._additionalElements == null) return false;

            // Both have additional elements
            foreach (var entry in _additionalElements)
            {
                if (!(other._additionalElements.TryGetValue(entry.Key, out var otherValue))) return false;
                if (entry.Value != otherValue) return false;
            }
            return true;
        }

        #region Operators

        /// <summary>
        /// Return a new composition that consists of this composition and <paramref name="c"/>
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        protected Composition AddComposition(Composition c)
        {
            var numC = _c + c._c;
            var numH = _h + c._h;
            var numN = _n + c._n;
            var numO = _o + c._o;
            var numS = _s + c._s;
            var numP = _p + c._p;

            if (_additionalElements == null && c._additionalElements == null)
                return new Composition(numC, numH, numN, numO, numS, numP);

            Dictionary<Atom, short> additionalElements = null;
            if (_additionalElements != null && c._additionalElements != null)
            {
                additionalElements = new Dictionary<Atom, short>(_additionalElements);
                foreach (var element in c._additionalElements)
                {
                    var atom = element.Key;

                    if (_additionalElements.TryGetValue(atom, out var numAtoms))
                    {
                        // atom was in _additionalElements
                        additionalElements[atom] = (short)(numAtoms + element.Value);
                    }
                    else
                    {
                        additionalElements[atom] = element.Value;
                    }
                }
            }
            else if (_additionalElements != null)
            {
                additionalElements = new Dictionary<Atom, short>(_additionalElements);
            }
            else if (c._additionalElements != null)
            {
                additionalElements = new Dictionary<Atom, short>(c._additionalElements);
            }

            var newComposition = new Composition(numC, numH, numN, numO, numS, numP, additionalElements);

            return newComposition;
        }

        /// <summary>
        /// Return a new composition that consists of this composition and <paramref name="c"/>
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public Composition Add(Composition c)
        {
            if (c is CompositionWithDeltaMass comWithDelta)
                return comWithDelta.Add(this);

            return AddComposition(c);
        }

        /// <summary>
        /// Overloaded '+' operator to add 2 compositions
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Composition operator +(Composition c1, Composition c2)
        {
            if (!(c1 is CompositionWithDeltaMass compWithDelta))
                return c1.Add(c2);

            return compWithDelta.Add(c2);
        }

        /// <summary>
        /// Return the inverse composition (all counts negative)
        /// </summary>
        /// <returns></returns>
        public Composition Negate()
        {
            if (_additionalElements == null)
                return new Composition(-_c, -_h, -_n, -_o, -_s, -_p);
            var additionalElements =
                _additionalElements.ToDictionary(element => element.Key, element => (short)(-element.Value));
            return new Composition(-_c, -_h, -_n, -_o, -_s, -_p, additionalElements);
        }

        /// <summary>
        /// Unary -
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Composition operator -(Composition c)
        {
            if (!(c is CompositionWithDeltaMass compWithDelta))
                return c.Negate();

            return compWithDelta.Negate();
        }

        /// <summary>
        /// Overloaded '-' operator to subtract one composition from another
        /// </summary>
        /// <param name="c1"></param>
        /// <param name="c2"></param>
        /// <returns></returns>
        public static Composition operator -(Composition c1, Composition c2)
        {
            return c1 + (-c2);
        }

        #endregion

        #region ToString and Parsing from String

        /// <inheritdoc />
        public override string ToString()
        {
            var basicCompositionStr = "C(" + C + ") H(" + H + ") N(" + N + ") O(" + O + ") S(" + S + (P != 0 ? ") P(" + P + ")" : ")");
            if (_additionalElements == null)
                return basicCompositionStr;

            var buf = new StringBuilder(basicCompositionStr);
            foreach (var element in _additionalElements)
            {
                buf.Append(" " + element.Key.Code + "(" + element.Value + ")");
            }
            return buf.ToString();
        }

        /// <summary>
        /// Return the composition as an empirical formula
        /// </summary>
        /// <returns></returns>
        public string ToPlainString()
        {
            var basicCompositionStr =
                (C == 0 ? "" : "C" + C)
                + (H == 0 ? "" : "H" + H)
                + (N == 0 ? "" : "N" + N)
                + (O == 0 ? "" : "O" + O)
                + (S == 0 ? "" : "S" + S)
                + (P == 0 ? "" : "P" + P);

            if (_additionalElements == null)
                return basicCompositionStr;

            var buf = new StringBuilder(basicCompositionStr);
            foreach (var element in _additionalElements)
            {
                buf.Append(element.Key.Code + element.Value);
            }
            return buf.ToString();
        }

        /// <summary>
        /// Parse a plain-string empirical formula, for example
        /// C2H3N1O1S
        /// </summary>
        /// <param name="plaincompositionStr"></param>
        /// <remarks>Empirical formula cannot have parentheses or spaces</remarks>
        /// <returns>Composition object, or null if a parse error</returns>
        public static Composition ParseFromPlainString(string plaincompositionStr)
        {
            if (!Regex.IsMatch(plaincompositionStr, @"^([A-Z][a-z]?-?\d*)+$")) return null;

            var unimodString = new StringBuilder();

            var matches = Regex.Matches(plaincompositionStr, @"[A-Z][a-z]?-?\d*");
            foreach (Match match in matches)
            {
                var element = match.Value;
                var atom = Regex.Match(element, @"[A-Z][a-z]?");
                var num = element.Substring(atom.Index + atom.Length);
                if (num.Length == 0) num = "1";
                if (unimodString.Length != 0) unimodString.Append(" ");
                unimodString.AppendFormat("{0}({1})", atom, num);
            }

            return Parse(unimodString.ToString());
        }

        /// <summary>
        /// Parse unimod-like composition string, for example
        /// H(117) C(77) N(17) O(26) S(2)
        /// </summary>
        /// <param name="compositionStr"></param>
        /// <remarks>Requires the use of parentheses for element counts.  Also requires whitespace (typically a space) between each element and count</remarks>
        /// <returns>Composition object, or null if a parse error</returns>
        public static Composition Parse(string compositionStr)
        {
            var c = 0;
            var h = 0;
            var n = 0;
            var o = 0;
            var s = 0;
            var p = 0;
            Dictionary<Atom, short> additionalElements = null;

            var token = compositionStr.Split();
            foreach (var e in token)
            {
                if (Regex.IsMatch(e, @"^\d*[a-zA-Z]+(\(-?\d+\))?$"))
                {
                    string element;
                    int num;
                    if (Regex.IsMatch(e, @"^\d*?[a-zA-Z]+$"))
                    {
                        element = e;
                        num = 1;
                    }
                    else
                    {
                        element = e.Substring(0, e.IndexOf('('));
                        num = int.Parse(e.Substring(e.IndexOf('(') + 1, e.LastIndexOf(')') - e.IndexOf('(') - 1));
                    }
                    if (element.Equals("C")) c += num;
                    else if (element.Equals("H")) h += num;
                    else if (element.Equals("N")) n += num;
                    else if (element.Equals("O")) o += num;
                    else if (element.Equals("S")) s += num;
                    else if (element.Equals("P")) p += num;
                    else
                    {
                        var atom = Atom.Get(element);
                        if (atom == null) return null;
                        if (additionalElements == null)
                            additionalElements = new Dictionary<Atom, short>();

                        if (additionalElements.TryGetValue(atom, out var currentAtomCount))
                            additionalElements[atom] = (short)(currentAtomCount + num);
                        else
                            additionalElements.Add(atom, (short)num);
                    }
                }
                else // illegal string
                {
                    return null;
                }
            }
            return new Composition(c, h, n, o, s, p, additionalElements);
        }

        #endregion

        #region Private Methods
        /// <summary>
        /// Gets the mono-isotopic mass
        /// </summary>
        private double GetMonoIsotopicMass()
        {
            var mass = _c * MassC + _h * MassH + _n * MassN + _o * MassO + _s * MassS + _p * MassP;
            if (_additionalElements != null) mass += _additionalElements.Sum(entry => entry.Key.Mass * entry.Value);
            return mass;
        }

        /// <summary>
        /// Gets the mono-isotopic nominal mass
        /// </summary>
        private int GetNominalMass()
        {
            var nominalMass = _c * NominalMassC + _h * NominalMassH + _n * NominalMassN + _o * NominalMassO + _s * NominalMassS + _p * NominalMassP;
            if (_additionalElements != null) nominalMass += _additionalElements.Sum(entry => entry.Key.NominalMass * entry.Value);
           return nominalMass;
        }

        #endregion

        #region Masses of Atoms

        private static readonly double MassC = Atom.Get("C").Mass;
        private static readonly double MassH = Atom.Get("H").Mass;
        private static readonly double MassN = Atom.Get("N").Mass;
        private static readonly double MassO = Atom.Get("O").Mass;
        private static readonly double MassS = Atom.Get("S").Mass;
        private static readonly double MassP = Atom.Get("P").Mass;
        //private static readonly double MassIsotope = Atom.Get("13C").AveragineMass - MassC;

        private static readonly int NominalMassC = Atom.Get("C").NominalMass;
        private static readonly int NominalMassH = Atom.Get("H").NominalMass;
        private static readonly int NominalMassN = Atom.Get("N").NominalMass;
        private static readonly int NominalMassO = Atom.Get("O").NominalMass;
        private static readonly int NominalMassS = Atom.Get("S").NominalMass;
        private static readonly int NominalMassP = Atom.Get("P").NominalMass;
        #endregion
    }
}
