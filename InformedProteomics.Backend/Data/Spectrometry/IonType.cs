using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class IonType
    {
        public class CtermIonType // temporarily defined by kyowon
        {
            
        }

        public class NtermIonType
        {
            
        }

        public class PrecursorIonType
        {
            
        }

        public int charge { get; private set; }
        
        public double GetMz(Composition cutComposition)
        {
            throw new System.NotImplementedException();
        }
    }
}
