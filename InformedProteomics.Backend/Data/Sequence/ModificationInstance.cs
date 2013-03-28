using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Sequence
{
    public class ModificationInstance
    {
        public ModificationInstance(IEnumerable<Modification> modifications)
        {
            Modifications = modifications;
            Composition = Composition.Zero;
            foreach(var modification in Modifications)
            {
                Composition += modification.Composition;
            }
        }

        private Composition _composition;
        private double _mass;

        public IEnumerable<Modification> Modifications { get; private set; }
        public Composition Composition 
        { 
            get { return _composition; }
            private set 
            { 
                _composition = value;
                _mass = value.GetMass();
            }
        }

        public int GetNumModifications()
        {
            return Modifications.Count();
        }

        public double GetMass()
        {
            return _mass;
        }
    }
}
