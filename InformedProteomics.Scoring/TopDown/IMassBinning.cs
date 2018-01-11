namespace InformedProteomics.Scoring.TopDown
{
    public interface IMassBinning
    {
        int GetBinNumber(double mass);

        double GetMass(int binNumber);
        double GetMassStart(int binNumber);
        double GetMassEnd(int binNumber);

        double MaxMass { get; }
        double MinMass { get; }
        int NumberOfBins { get; }

        bool Filtered { get; }
    }
}
