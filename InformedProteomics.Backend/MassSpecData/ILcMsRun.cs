using System.Collections.Generic;

namespace InformedProteomics.Backend.MassSpecData
{
    public interface ILcMsRun: IChromatogramExtractor, ISpectrumExtractor
    {
        int MinLcScan { get; }
        int MaxLcScan { get; }

        double GetElutionTime(int scanNum);
        int GetMsLevel(int scanNum);
        int GetPrevScanNum(int scanNum, int msLevel);
        int GetNextScanNum(int scanNum, int msLevel);
        IList<int> GetScanNumbers(int msLevel);
    }
}
