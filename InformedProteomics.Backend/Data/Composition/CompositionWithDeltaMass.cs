using System;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Composition
{
    public class CompositionWithDeltaMass: Composition
    {
        public CompositionWithDeltaMass(double mass, int nominalMass) : base(0, 0, 0, 0, 0)
        {
            _deltaMass = mass;
            _deltaNominalMass = nominalMass;
        }

        public CompositionWithDeltaMass(double mass)
            : this(mass, (int)Math.Round(mass * Constants.RescalingConstant))
        {
        }

        public override double Mass
        {
            get { return base.Mass + _deltaMass; }
        }

        public override int NominalMass
        {
            get { return base.NominalMass + _deltaNominalMass; }
        }

        public static CompositionWithDeltaMass operator +(CompositionWithDeltaMass c1, CompositionWithDeltaMass c2)
        {
            var comp = (Composition)c1 + (Composition)c2;
            return new CompositionWithDeltaMass(comp, c1._deltaMass+c2._deltaMass, c1._deltaNominalMass+c2._deltaNominalMass);
        }

        public static CompositionWithDeltaMass operator +(Composition c1, CompositionWithDeltaMass c2)
        {
            var comp = c1 + (Composition)c2;
            return new CompositionWithDeltaMass(comp, c2._deltaMass, c2._deltaNominalMass);
        }

        public static CompositionWithDeltaMass operator +(CompositionWithDeltaMass c1, Composition c2)
        {
            return c2 + c1;
        }

        /// <summary>
        /// Unary -
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public static CompositionWithDeltaMass operator -(CompositionWithDeltaMass c)
        {
            var comp = -((Composition) c);
            return new CompositionWithDeltaMass(comp, -c._deltaMass, -c._deltaNominalMass);
        }

        public static CompositionWithDeltaMass operator -(CompositionWithDeltaMass c1, CompositionWithDeltaMass c2)
        {
            return c1 + (-c2);
        }

        public static CompositionWithDeltaMass operator -(Composition c1, CompositionWithDeltaMass c2)
        {
            return c1 + (-c2);
        }

        public static CompositionWithDeltaMass operator -(CompositionWithDeltaMass c1, Composition c2)
        {
            return c1 + (-c2);
        }

        public override string ToString()
        {
            return base.ToString() + " " + string.Format("{0:N3}", _deltaMass);
        }

        protected bool Equals(CompositionWithDeltaMass other)
        {
            return base.Equals(other) && _deltaMass.Equals(other._deltaMass);
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != GetType()) return false;
            return Equals((CompositionWithDeltaMass)obj);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                return (base.GetHashCode() * 397) ^ _deltaMass.GetHashCode();
            }
        }

        private CompositionWithDeltaMass(Composition composition, double deltaMass, int deltaNominalMass): base(composition)
        {
            _deltaMass = deltaMass;
            _deltaNominalMass = deltaNominalMass;
        }

        private readonly double _deltaMass;
        private readonly int _deltaNominalMass;

    }
}
