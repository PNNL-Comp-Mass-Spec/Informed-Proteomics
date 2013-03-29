﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class NeutralLoss
    {
        public NeutralLoss(string name, Composition composition)
        {
            Name = name;
            Composition = composition;
        }

        public string Name { get; private set; }
        public Composition Composition { get; private set; }

        public static readonly NeutralLoss NoLoss = new NeutralLoss("", Composition.Zero);
        public static readonly NeutralLoss H2O = new NeutralLoss("-H2O", Composition.H2O);
        public static readonly NeutralLoss NH3 = new NeutralLoss("-NH3", Composition.NH3);

        public static readonly IEnumerable<NeutralLoss> CommonNeutralLosses = new List<NeutralLoss> { NoLoss, H2O, NH3 };
    }
}