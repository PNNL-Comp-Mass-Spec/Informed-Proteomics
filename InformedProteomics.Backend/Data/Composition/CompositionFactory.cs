using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;

namespace InformedProteomics.Backend.Data.Composition
{
    public abstract class CompositionFactory
    {
        public Composition GetComposition(int c, int h, int n, int o, int s, int p)
        {
            return null;
        }

        public Composition GetComposition(int c, int h, int n, int o, int s)
        {
            return GetComposition(c, h, n, o, s, 0);
        }

        public Composition GetComposition(Data.Composition.Composition composition)
        {
            return null;
        }

        public Composition GetComposition(int c, int h, int n, int o, int s, int p, IEnumerable<Tuple<Atom, short>> additionalElements)
        {
            return null;
        }
    }
}

