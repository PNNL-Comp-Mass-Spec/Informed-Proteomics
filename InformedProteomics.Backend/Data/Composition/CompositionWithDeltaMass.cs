using System;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Composition
{
    /// <summary>
    /// A composition with a delta mass
    /// </summary>
    public class CompositionWithDeltaMass : Composition
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="mass"></param>
        /// <param name="nominalMass"></param>
        public CompositionWithDeltaMass(double mass, int nominalMass) : base(0, 0, 0, 0, 0)
        {
            _deltaMass = mass;
            _deltaNominalMass = nominalMass;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="mass"></param>
        public CompositionWithDeltaMass(double mass)
            : this(mass, (int)Math.Round(mass * Constants.RescalingConstant))
        {
        }

        /// <summary>
        /// Monoisotopic mass
        /// </summary>
        public override double Mass => base.Mass + _deltaMass;

        /// <summary>
        /// Nominal mass
        /// </summary>
        public override int NominalMass => base.NominalMass + _deltaNominalMass;

        /// <summary>
        /// Return a new composition that consists of this composition added to <paramref name="c"/>
        /// </summary>
        /// <param name="c"></param>
        /// <returns></returns>
        public new Composition Add(Composition c)
        {
            if (!(c is CompositionWithDeltaMass comWithDelta))
            {
                return new CompositionWithDeltaMass(AddComposition(c), _deltaMass, _deltaNominalMass);
            }

            return new CompositionWithDeltaMass(AddComposition(c), _deltaMass + comWithDelta._deltaMass, _deltaNominalMass + comWithDelta._deltaNominalMass);
        }

        /// <summary>
        /// Return the negation of this composition
        /// </summary>
        /// <returns></returns>
        public new Composition Negate()
        {
            return new CompositionWithDeltaMass(base.Negate(), -_deltaMass, -_deltaNominalMass);
        }

        /// <inheritdoc />
        public override string ToString()
        {
            return base.ToString() + " " + string.Format("{0:N3}", _deltaMass);
        }

        /// <summary>
        /// Check 2 CompositionWithDeltaMasses for equality
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        protected bool Equals(CompositionWithDeltaMass other)
        {
            return base.Equals(other) && _deltaMass.Equals(other._deltaMass);
        }

        /// <inheritdoc />
        public override bool Equals(object obj)
        {
            if (obj == null)
            {
                return false;
            }

            if (ReferenceEquals(this, obj))
            {
                return true;
            }

            if (obj.GetType() != GetType())
            {
                return false;
            }

            return Equals((CompositionWithDeltaMass)obj);
        }

        /// <inheritdoc />
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
