namespace InformedProteomics.Backend.Utils
{
    public class GeneratingFunctionForSkyline
    {
        //public GeneratingFunctionForSkyline(
        //    double minPrecursorIonMass, 
        //    double maxPrecursorIonMass, 
        //    double[] minProductIonMasses, 
        //    double[] maxProductIonMasses
        //    )
        //{
        //    _peptideNominalMass = peptideNominalMass;
        //    _fragmentMasses = fragmentMasses;
        //    _hist = new double[peptideNominalMass + 1, _fragmentMasses.Length + 1];
        //    _aminoAcidIntMasses = AminoAcidMasses.Select(GetBinNumber).ToArray();

        //}

        //private void ComputeGeneratingFunction()
        //{
        //    _hist[0, 0] = 1.0;
        //    var score = new int[_peptideNominalMass];
        //    foreach (var fragmentMass in _fragmentMasses) score[fragmentMass] = 1;
        //    for (var nominalMass = 1; nominalMass <= _peptideNominalMass; nominalMass++)
        //    {
        //        var curScore = score[nominalMass];
        //        foreach (var aminoAcidMass in AminoAcidMasses)
        //        {
        //            var prevMass = nominalMass - aminoAcidMass;
        //            if(prevMass < 0) continue;
        //            for (var s = 0; s <= _fragmentMasses.Length; s++)
        //            {
                        
        //            }
        //        }
        //    }
        //}

        //public static int GetBinNumber(double doubleValue)
        //{
        //    return (int) Math.Round(doubleValue*RescalingConstant);
        //}

        //private readonly int[] _aminoAcidIntMasses;

        //private readonly int _min
        //private readonly int _peptidenIntMass;
        //private readonly int[] _fragmentMasses;
        //private readonly int _numTransitions;
        //private readonly double[,] _hist;

        //private static readonly double[] AminoAcidMasses;
        //private static readonly double[] AminoAcidProbabilities;
        //private static readonly double[] PrefixIonOffsets
        //private const double RescalingConstant = 274.335215;

        //static GeneratingFunctionForSkyline()
        //{
        //    AminoAcidMasses = new[]
        //    {
        //        57.021464, 71.037114, 87.032029, 97.052764, 99.068414,
        //        101.04768, 103.00919, 113.08406, 113.08406, 114.04293,
        //        115.02694, 128.05858, 128.09496, 129.04259, 131.04048,
        //        137.05891, 147.06841, 156.10111, 163.06333, 186.07931
        //    };
        //    AminoAcidProbabilities = new double[20];
        //    for (var i = 0; i < AminoAcidProbabilities.Length; i++) AminoAcidProbabilities[i] = 0.05;
        //}
    }
}
