using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public enum ActivationMethod : byte
    {
        CID,
        ETD,
        HCD,
        ECD,
        PQD,
        Unknown
    }
}
