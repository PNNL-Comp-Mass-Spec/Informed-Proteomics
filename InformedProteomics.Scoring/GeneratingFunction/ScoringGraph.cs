using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Scoring.GeneratingFunction
{
    public class ScoringGraph : IScoringGraph
    {
        public ScoringGraph(IList<DeconvolutedPeak> peaks, IList<int> peakScore, MzComparerWithBinning comparer)
        {
            for(var i = 0; i < peaks.Count; i++)
            {
                var massBinIndex = comparer.GetBinNumber(peaks[i].Mass);
                var massPeakScore = peakScore[i];

                // implement here
            }
        }

        public ScoringGraph(ProductSpectrum spectrum)
        {
            // implment later
        }
        
        public int GetNodeScore(int nodeIndex)
        {
            throw new NotImplementedException();
        }

        public IEnumerable<ScoringGraphEdge> GetEdges(int nodeIndex)
        {
            throw new NotImplementedException();
        }

        public int GetNumNodes()
        {
            throw new NotImplementedException();
        }
    }
}
