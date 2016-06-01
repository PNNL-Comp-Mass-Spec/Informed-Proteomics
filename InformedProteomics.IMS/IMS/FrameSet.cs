using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.IMS
{
    public class FrameSet
    {
        public int FrameID
        {
            get
            {
                throw new System.NotImplementedException();
            }
            set
            {
            }
        }

        // Usually the first frame
        public Frame GetMS1Frame()
        {
            throw new System.NotImplementedException();
        }

        public Frame[] GetMS2Frames()
        {
            throw new System.NotImplementedException();
        }

        public Frame GetFrame(int frameNum)
        {
            throw new System.NotImplementedException();
        }

        public List<IMSPrecursor> GetPrecursors(Ion ion)
        {
            throw new System.NotImplementedException();
        }
    }
}
