using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class Composition : IMolecule
    {
        public static readonly Composition Zero = new Composition(0, 0, 0, 0, 0);
		public static readonly Composition H2O = new Composition(0, 2, 0, 1, 0);
        public static readonly Composition NH3 = new Composition(0, 3, 1, 0, 0);
        public static readonly Composition NH2 = new Composition(0, 2, 1, 0, 0);
        public static readonly Composition OH = new Composition(0, 1, 0, 1, 0);
        public static readonly Composition CO = new Composition(1, 0, 0, 1, 0);
        public static readonly Composition Hydrogen = new Composition(0, 1, 0, 0, 0);

        public Composition(int c, int h, int n, int o, int s)
        {
            _c = (short)c;
            _h = (short)h;
            _n = (short)n;
            _o = (short)o;
            _s = (short)s;
            _additionalElements = null;
        }

        public Composition(int c, int h, int n, int o, int s, int p)
            : this(c, h, n, o, s, new Tuple<Atom, short>(AtomP, (short)p))
        {
        }
        
        public Composition(Composition composition): this(composition.C, composition.H, 
            composition.N, composition.O, composition.S)
        {
            if (composition._additionalElements != null)
            {
                _additionalElements = new Dictionary<Atom, short>(composition._additionalElements);
            }
        }

        public Composition(int c, int h, int n, int o, int s, Tuple<Atom, short> additionalElement)
            : this(c, h, n, o, s)
        {
            _additionalElements = new Dictionary<Atom, short> {{additionalElement.Item1, additionalElement.Item2}};
        }

        public Composition(int c, int h, int n, int o, int s, IEnumerable<Tuple<Atom, short>> additionalElements)
            : this(c, h, n, o, s)
        {
            _additionalElements = new Dictionary<Atom, short>();
            foreach (var element in additionalElements)
            {
                _additionalElements.Add(element.Item1, element.Item2);
            }
        }

        private Composition(int c, int h, int n, int o, int s, Dictionary<Atom, short> additionalElements)
            : this(c, h, n, o, s)
        {
            _additionalElements = additionalElements;
        }

        #region Properties

        public short C
        {
            get { return _c; }
        }

        public short H
        {
            get { return _h; }
        }

        public short N
        {
            get { return _n; }
        }

        public short O
        {
            get { return _o; }
        }

        public short S
        {
            get { return _s; }
        }

        #endregion

        #region Private members

        private readonly Dictionary<Atom, short> _additionalElements;
        private readonly short _c;
        private readonly short _h;
        private readonly short _n;
        private readonly short _o;
        private readonly short _s;

        #endregion


        /// <summary>
        /// Gets the mono-isotopic mass
        /// </summary>
        public double GetMass()
        {
            double mass = _c*MassC + _h*MassH + _n*MassN + _o*MassO + _s*MassS;
            if (_additionalElements != null)
                mass += _additionalElements.Sum(entry => entry.Key.Mass*entry.Value);
            return mass;
        }

        /// <summary>
        /// Gets the m/z of ith isotope
        /// </summary>
        /// <param name="isotopeIndex">isotope index. 0 means mono-isotope, 1 means 2nd isotope, etc.</param>
        /// <returns></returns>
        public double GetIsotopeMass(int isotopeIndex)
        {
            return GetMass() + isotopeIndex * MassIsotope;
        }

        /// <summary>
        /// Gets the mono-isotopic nominal mass
        /// </summary>
        public int GetNominalMass()
        {
            int nominalMass = _c*NominalMassC + _h*NominalMassH + _n*NominalMassN + _o*NominalMassO + _s*NominalMassS;
            if (_additionalElements != null)
                nominalMass += _additionalElements.Sum(entry => entry.Key.NominalMass * entry.Value);
            return nominalMass;
        }

        public Composition GetComposition()
        {
            return this;
        }

        public override int GetHashCode()
        {
            var hashCode = _c * 0x01000000 + _h * 0x00010000 + _n * 0x00000400 + _o * 0x00000010 + _s;
            if (_additionalElements == null)
                return hashCode;
            return hashCode + _additionalElements.Sum(element => element.Key.GetHashCode() * element.Value);
        }

        public override bool Equals(object obj)
        {
            var other = obj as Composition;
            if (other == null)
                return false;
            if (_c != other.C || _h != other.H || _n != other.N || _o != other.O || _s != other.S)
                return false;

            if (_additionalElements != null)
            {
                if (other._additionalElements == null)
                {
                    return false;
                }
                foreach (var entry in _additionalElements)
                {
                    short num;
                    if (_additionalElements.TryGetValue(entry.Key, out num))
                    {
                        if (entry.Value != num)
                            return false;
                    }
                    else
                    {
                        return false;
                    }
                }
                return true;
            }
            return other._additionalElements == null;
        }

        public float[] GetIsotopomerEnvelop()
        {
            return GetApproximatedIsotopomerEnvelop();
        }

        public float[] GetApproximatedIsotopomerEnvelop()
        {
            var mean = _c*(1 - ProbC12) + _h*(1 - ProbH1) + _n*(1 - ProbN14) + _o*(1 - ProbO16) + _s*(1 - ProbS32);
          
            var dist = new float[(int)mean + 4];
            var exp = Math.Exp(-mean);
            for (var i = 0; i < dist.Length; i++)
            {
                dist[i] = (float) (Math.Pow(mean, i)*exp/MathNet.Numerics.SpecialFunctions.Factorial(i));

            }
            float max = dist.Concat(new float[] {0}).Max();
            for (var i = 0; i < dist.Length; i++)
            {
                dist[i] = dist[i] / max;
            }
            return dist;
        }

        public static Composition operator +(Composition c1, Composition c2)
        {
            if (c1 == Composition.Zero)
                return c2;
            if (c2 == Composition.Zero)
                return c1;
            int numC = c1._c + c2._c;
            int numH = c1._h + c2._h;
            int numN = c1._n + c2._n;
            int numO = c1._o + c2._o;
            int numS = c1._s + c2._s;

            if(c1._additionalElements == null && c2._additionalElements == null)
                return new Composition(numC, numH, numN, numO, numS);

            Dictionary<Atom, short> additionalElements = null;
            if (c1._additionalElements != null && c2._additionalElements != null)
            {
                additionalElements = new Dictionary<Atom, short>(c1._additionalElements);
                foreach (var element in c2._additionalElements)
                {
                    var atom = element.Key;
                    short numAtoms;
                    if (c1._additionalElements.TryGetValue(atom, out numAtoms))
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
            else if (c1._additionalElements != null)
            {
                additionalElements = new Dictionary<Atom, short>(c1._additionalElements);
            }
            else if(c2._additionalElements != null)
            {
                additionalElements = new Dictionary<Atom, short>(c2._additionalElements);
            }
            return new Composition(numC, numH, numN, numO, numS, additionalElements);
        }

        /// <summary>
        /// Unary -
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static Composition operator -(Composition c)
        {
            if(c._additionalElements == null)
                return new Composition(-c._c, -c._h, -c._n, -c._o, -c._s);
            var additionalElements = 
                c._additionalElements.ToDictionary(element => element.Key, element => (short) (-element.Value));
            return new Composition(-c._c, -c._h, -c._n, -c._o, -c._s, additionalElements);
        }

        public static Composition operator -(Composition c1, Composition c2)
        {
            return c1 + (-c2);
        }

		public override string ToString()
		{
		    string basicCompositionStr = "C(" + C + ") H(" + H + ") N(" + N + ") O(" + O + ") S(" + S + ")";
            if(_additionalElements == null)
                return basicCompositionStr;

            var buf = new StringBuilder(basicCompositionStr);
		    foreach (var element in _additionalElements)
		    {
		        buf.Append(" " + element.Key.Code + "(" + element.Value + ")");
		    }
		    return buf.ToString();
		}

        #region Masses of Atoms

        private static readonly double MassC = Atom.Get("C").Mass;
        private static readonly double MassH = Atom.Get("H").Mass;
        private static readonly double MassN = Atom.Get("N").Mass;
        private static readonly double MassO = Atom.Get("O").Mass;
        private static readonly double MassS = Atom.Get("S").Mass;
        private static readonly double MassIsotope = Atom.Get("13C").Mass - MassC;

        private static readonly Atom AtomP = Atom.Get("P");

        private static readonly int NominalMassC = Atom.Get("C").NominalMass;
        private static readonly int NominalMassH = Atom.Get("H").NominalMass;
        private static readonly int NominalMassN = Atom.Get("N").NominalMass;
        private static readonly int NominalMassO = Atom.Get("O").NominalMass;
        private static readonly int NominalMassS = Atom.Get("S").NominalMass;

        #endregion

        #region Isotope probabilities of Atoms // added by kyowon

        private const double ProbC12 = .9890;
        private const double ProbH1 = .99985;
        private const double ProbN14 = .99636;
        private const double ProbO16 = .99762;
        private const double ProbS32 = .95029;

        #endregion
    }
}
