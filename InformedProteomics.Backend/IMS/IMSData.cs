using System.Collections.Generic;
using InformedProteomics.Backend.Data.Science;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.IMS
{
    public class IMSData
    {
        public int NumFrameSets
        {
            get
            {
                throw new System.NotImplementedException();
            }
            set
            {
            }
        }

        public int NumFramesPerFrameSet
        {
            get
            {
                throw new System.NotImplementedException();
            }
            set
            {
            }
        }
    
        // MS1 only
        public XIC GetXIC(Ion ion)
        {
            throw new System.NotImplementedException();
        }

        public void GetFrameSet(int frameSetID)
        {
            throw new System.NotImplementedException();
        }

        public List<Frame> GetMS1Frames()
        {
            throw new System.NotImplementedException();
        }

        public List<FrameSet> GetFrameSetList()
        {
            throw new System.NotImplementedException();
        }
    }
}
