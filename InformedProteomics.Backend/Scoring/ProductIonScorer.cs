using System;
using InformedProteomics.Backend.Data.Results;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Scoring
{
    class ProductIonScorer
    {
        public Dictionary<int, Dictionary<IonType, double[]>> IonXICs { get; private set; }
        public Dictionary<int, Dictionary<IonType, double>[]> Spectra { get; private set; }
        public float Score { get; private set; }
        public Sequence Peptide { get; private set; }
        public DatabaseSubTargetResult PrecursorResultRep { get; private set; }
        
        public ProductIonScorer(DatabaseMultipleSubTargetResult matchedResult, Sequence peptide)
        {
            PrecursorResultRep = matchedResult.PrecursorResultRep;
            GetFragmentResults(matchedResult);
            Peptide = peptide;
            Score = GetScore();
        }

        private void GetFragmentResults(DatabaseMultipleSubTargetResult matchedResult)
        {
            IonXICs = new Dictionary<int, Dictionary<IonType, double[]>>();
            Spectra = new Dictionary<int, Dictionary<IonType, double>[]>();

            foreach (var fragmentTargetResult in matchedResult.FragmentResultList)
            {
                var fragment = fragmentTargetResult.DatabaseFragmentTarget.Fragment;
                var residueNumber = fragment.ResidueNumber;
                //Console.WriteLine(residueNumber + "\t");
                var ion = new IonType(fragment.IonSymbol, fragment.ChargeState);
                //Console.WriteLine(ion);
                if (!IonXICs.ContainsKey(residueNumber))
                {
                    IonXICs.Add(residueNumber, new Dictionary<IonType, double[]>());
                    Spectra.Add(residueNumber, new Dictionary<IonType, double>[fragmentTargetResult.XYData.Yvalues.Length]);
                    for (var i = 0; i < Spectra[residueNumber].Length;i++ )
                        Spectra[residueNumber][i] = new Dictionary<IonType, double>();
                }

                if (!IonXICs[residueNumber].ContainsKey(ion)){
                    IonXICs[residueNumber].Add(ion, fragmentTargetResult.XYData.Yvalues);
                }
                var spectraPerFragment = Spectra[residueNumber];
                
                /*Console.WriteLine(residueNumber + "\t" + PrecursorResultRep.XYData.Xvalues.Length + "\t" + fragmentTargetResult.XYData.Yvalues.Length + "\t" + fragmentTargetResult.XYData.Xvalues.Length + "\t" + spectraPerFragment.Length);
                foreach (var t in PrecursorResultRep.XYData.Xvalues)
                    Console.Write(t + " ");
                Console.WriteLine();
                foreach (var t in fragmentTargetResult.XYData.Xvalues)
                    Console.Write(t+" ");
                Console.WriteLine();
                foreach (var t in fragmentTargetResult.XYData.Yvalues)
                    Console.Write(t + " ");
                Console.WriteLine();
                Console.WriteLine();
                */
                for (var i = 0; i < fragmentTargetResult.XYData.Yvalues.Length;i++)
                {
                    //Console.WriteLine(spectraPerFragment.Length + "\t" + i);
                    if (!spectraPerFragment[i].ContainsKey(ion))
                    {
                        spectraPerFragment[i].Add(ion, fragmentTargetResult.XYData.Yvalues[i]);
                    }
                }
            }
        }


        private float GetScore()
        {
            var score = 0f;//
            for (var residueNumber = 1; residueNumber < Peptide.Count; residueNumber++)
            {
                float subScore;
                if (IonXICs.ContainsKey(residueNumber))
                {
                    var specScorer = new FragmentSpectrumScorer(Spectra[residueNumber]); // check if +1 or -1..


                    
                    //foreach(var q in specScorer.UsedIons) Console.Write(q.ToString() + " ");
                   // Console.WriteLine();
                   // foreach(var q in IonXICs[residueNumber].Keys) Console.Write(q.ToString()+" ");
                   // Console.WriteLine();
                    if (specScorer.UsedIons.Count == 0)
                    {
                        subScore = -3;
                    }
                    else
                    {
                        var fragXICScorer = new FragmentXICScorer(IonXICs[residueNumber], specScorer.UsedIons,
                                                                  PrecursorResultRep.XYData.Yvalues);
                        subScore = fragXICScorer.Score + specScorer.Score;

                        Console.WriteLine("Frag Score : " + residueNumber + "\t" + specScorer.Score + "\t" + fragXICScorer.RawScore + "\t" + fragXICScorer.Score + "\t" + specScorer.UsedIons.Count + "\t" + IonXICs[residueNumber].Keys.Count);

                    }
                    //Console.WriteLine("Frag Score : " + residueNumber + "\t" + fragXICScorer.Score + "\t" + specScorer.Score + "\t" + specScorer.UsedIons.Count);
                }
                else
                {
                    subScore = -3;
                }

                score += subScore;
                
            }
            Console.WriteLine();
            return score;                        
        }

    }
}
