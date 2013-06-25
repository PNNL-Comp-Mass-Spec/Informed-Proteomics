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

        public Composition(Composition composition)
            : this(composition.C, composition.H,
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
            _additionalElements = new Dictionary<Atom, short> { { additionalElement.Item1, additionalElement.Item2 } };
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
            double mass = _c * MassC + _h * MassH + _n * MassN + _o * MassO + _s * MassS;
            if (_additionalElements != null)
                mass += _additionalElements.Sum(entry => entry.Key.Mass * entry.Value);
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
            int nominalMass = _c * NominalMassC + _h * NominalMassH + _n * NominalMassN + _o * NominalMassO + _s * NominalMassS;
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
            return GetApproximatedIsotopomerEnvelop(10);
        }

        private float GetIsotopeProbability(int[] number, double[] means, double[] exps)
        {
            var prob = 1.0f;
            for (var i = 1; i <= Math.Min(ProbC.Length - 1, number.Length); i++)
            {
                var mean = means[i];
                var exp = exps[i];
                if (number[i - 1] == 0) prob *= (float)exp;
                else prob *= (float)(Math.Pow(mean, number[i - 1]) * exp / MathNet.Numerics.SpecialFunctions.Factorial(number[i - 1]));
            }
            return prob;
        }

        private static int[][][] GetPossibleIsotopeCombinations(int max) // called just once. 
        {
            var comb = new List<int[]>[max + 1];
            var maxIsotopeNumberInElement = ProbC.Length - 1;
            comb[0] = new List<int[]> { (new int[maxIsotopeNumberInElement]) };

            for (var n = 1; n <= max; n++)
            {
                comb[n] = new List<int[]>();
                for (var j = 1; j <=maxIsotopeNumberInElement; j++)
                {
                    var index = n - j;
                    if (index < 0) continue;
                    foreach (var t in comb[index])
                    {
                        var add = new int[maxIsotopeNumberInElement];
                        add[j - 1]++;
                        for (var k = 0; k < t.Length; k++)
                            add[k] += t[k];
                        var toAdd = comb[n].Select(v => !v.Where((t1, p) => t1 != add[p]).Any()).All(equal => !equal);
                        if(toAdd) comb[n].Add(add);
                    }
                }
            }
            var possibleIsotopeCombinations = new int[max][][];
            for (var i = 0; i < possibleIsotopeCombinations.Length; i++)
            {
                possibleIsotopeCombinations[i] = new int[comb[i].Count][];
                var j = 0;
                foreach (var t in comb[i])
                {
                    possibleIsotopeCombinations[i][j++] = t;
                }
            }
            return possibleIsotopeCombinations;
        }

        public float[] GetApproximatedIsotopomerEnvelop()
        {
            return GetApproximatedIsotopomerEnvelop(7);
        }

        public float[] GetApproximatedIsotopomerEnvelop(int maxIsotope)
        {
            var dist = new float[Math.Min(maxIsotope, PossibleIsotopeCombinations.Length)];
            var means = new double[PossibleIsotopeCombinations[0][0].Length + 1];
            var exps = new double[means.Length];
            for (var i = 0; i < means.Length; i++) // precalculate means and thier exps
            {
                means[i] = _c * ProbC[i] + _h * ProbH[i] + _n * ProbN[i] + _o * ProbO[i] + _s * ProbS[i];
                exps[i] = Math.Exp(means[i]);
            }
            for (var i = 0; i < dist.Length; i++)
            {
                foreach (var isopeCombinations in PossibleIsotopeCombinations[i])
                {
                    dist[i] += GetIsotopeProbability(isopeCombinations, means, exps);
                }
            }
            var max = dist.Concat(new float[] { 0 }).Max();
            var maxPassed = false;
            var cutIndex = dist.Length - 1;
            for (var i = 0; i < dist.Length; i++)
            {
                dist[i] = dist[i] / max;
                if (dist[i] >= 1) maxPassed = true;
                if (!maxPassed || dist[i] >= 0.01) continue;
                cutIndex = i;
                break;
            }

            var truncatedDist = new float[cutIndex + 1];
            for (var i = 0; i < truncatedDist.Length; i++)
                truncatedDist[i] = dist[i];

            return truncatedDist;
        }

        public int GetMostAbundantIsotopeZeroBasedIndex()
        {
            var index = 0;
            var isotopeEnvelope = GetApproximatedIsotopomerEnvelop();
            var max = 0.0f;
            for (var i = 0; i < isotopeEnvelope.Length; i++)
            {
                if (!(max < isotopeEnvelope[i])) continue;
                max = isotopeEnvelope[i];
                index = i;
                if (max >= 1) break;
            }
            return index;
        }

        public static Composition operator +(Composition c1, Composition c2)
        {
            // ReSharper disable PossibleUnintendedReferenceComparison
            if (c1 == Zero)
                return c2;
            if (c2 == Zero)
                return c1;
            // ReSharper restore PossibleUnintendedReferenceComparison
            int numC = c1._c + c2._c;
            int numH = c1._h + c2._h;
            int numN = c1._n + c2._n;
            int numO = c1._o + c2._o;
            int numS = c1._s + c2._s;

            if (c1._additionalElements == null && c2._additionalElements == null)
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
            else if (c2._additionalElements != null)
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
            if (c._additionalElements == null)
                return new Composition(-c._c, -c._h, -c._n, -c._o, -c._s);
            var additionalElements =
                c._additionalElements.ToDictionary(element => element.Key, element => (short)(-element.Value));
            return new Composition(-c._c, -c._h, -c._n, -c._o, -c._s, additionalElements);
        }

        public static Composition operator -(Composition c1, Composition c2)
        {
            return c1 + (-c2);
        }

        public override string ToString()
        {
            string basicCompositionStr = "C(" + C + ") H(" + H + ") N(" + N + ") O(" + O + ") S(" + S + ")";
            if (_additionalElements == null)
                return basicCompositionStr;

            var buf = new StringBuilder(basicCompositionStr);
            foreach (var element in _additionalElements)
            {
                buf.Append(" " + element.Key.Code + "(" + element.Value + ")");
            }
            return buf.ToString();
        }

        public string ToPlainString()
        {
            return
                (C == 0 ? "" : "C" + C)
                + (H == 0 ? "" : "H" + H)
                + (N == 0 ? "" : "N" + N)
                + (O == 0 ? "" : "O" + O)
                + (S == 0 ? "" : "S" + S);
        }

        //added by kyowon, incomplete for additional Elementss
        public static Composition Parse(string s)
        {
            var t = s.Split(' ');
            if (t.Length < 5) return null;
            var basicComposition = new int[5];
            for (var i = 0; i < basicComposition.Length; i++)
            {
                basicComposition[i] = int.Parse(t[i].Substring(t[i].IndexOf('(') + 1, t[i].IndexOf(')') - t[i].IndexOf('(') - 1));
            }
            if (t.Length == basicComposition.Length)
                return new Composition(basicComposition[0], basicComposition[1], basicComposition[2], basicComposition[3], basicComposition[4]);
            throw new System.NotImplementedException();
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

        #region Isotope probabilities of Atoms // added by kyowon - the numbers of elements should be the same. Now it is 4.

        private static readonly double[] ProbC = new[] { .9893, 0.0107, 0, 0 };
        private static readonly double[] ProbH = new[] { .999885, .000115, 0, 0 };
        private static readonly double[] ProbN = new[] { 0.99632, 0.00368, 0, 0};
        private static readonly double[] ProbO = new[] { 0.99757, 0.00038, 0.00205, 0};
        private static readonly double[] ProbS = new[] { 0.9493, 0.0076, 0.0429, 0.0002 };
        #endregion

        #region Possible combinations of isotopes // added by kyowon

        private static readonly int[][][] PossibleIsotopeCombinations = GetPossibleIsotopeCombinations(100);

        #endregion
    }
}
