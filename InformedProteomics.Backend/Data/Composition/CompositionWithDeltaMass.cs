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

        public new Composition Add(Composition c)
        {
            var comWithDelta = c as CompositionWithDeltaMass;
            return comWithDelta == null ? 
                new CompositionWithDeltaMass(AddComposition(c), _deltaMass, _deltaNominalMass)
                : new CompositionWithDeltaMass(AddComposition(c), _deltaMass + comWithDelta._deltaMass, _deltaNominalMass + comWithDelta._deltaNominalMass);
        }

        public new Composition Negate()
        {
            return new CompositionWithDeltaMass(base.Negate(), -_deltaMass, -_deltaNominalMass);
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

        private CompositionWithDeltaMass(Composition composition, double deltaMass, int deltaNominalMass)
            : base(composition)
        {
            _deltaMass = deltaMass;
            _deltaNominalMass = deltaNominalMass;
        }

        private readonly double _deltaMass;
        private readonly int _deltaNominalMass;

    }
}
