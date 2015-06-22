using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;


namespace InformedProteomics.Backend.MassFeature
{
    public interface INodeComparer<T>
    {
        bool SameCluster(T node1, T node2);
    }
}
