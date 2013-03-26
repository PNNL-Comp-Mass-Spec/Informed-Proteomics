using System;
using InformedProteomics.Backend.Data.Results;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Sequence;

namespace InformedProteomics.Backend.Scoring
{
    class ProductIonScorer
    {
        public FragmentXICSet[] FragmentXICSets { get; private set; } // residue
        public FragmentSpectrum[] FragmentSpectra { get; private set; } // residue
        public float Score { get; private set; }
        public Sequence Peptide { get; private set; }
        public DatabaseSubTargetResult PrecursorResultRep { get; private set; }
        private HashSet<int> _matchedResidues;
        private int _apexIndex;
        private static readonly Dictionary<FragmentSpectrumParameter, float> NoXICScore = new Dictionary<FragmentSpectrumParameter, float>();
        
        public ProductIonScorer(DatabaseMultipleSubTargetResult matchedResult, Sequence peptide)
        {
            PrecursorResultRep = matchedResult.PrecursorResultRep;
            Peptide = peptide;
            GetFragmentResults(matchedResult);
            Score = GetScore();
        }
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

        private void GetFragmentResults(DatabaseMultipleSubTargetResult matchedResult)
        {
            FragmentXICSets = new FragmentXICSet[Peptide.Count];
            FragmentSpectra = new FragmentSpectrum[Peptide.Count];
            _matchedResidues = new HashSet<int>();
            for (var i = 0; i < FragmentXICSets.Length; i++)
            {
                FragmentXICSets[i] = new FragmentXICSet();
            }

            for (var i = 0; i < FragmentSpectra.Length; i++)
            {
                FragmentSpectra[i] = new FragmentSpectrum();
            }
            var precursorIon = IonType.GetPrecursorIon(PrecursorResultRep.DatabaseSubTarget.ChargeState);

            var apexIntensity = 0.0;
            for (var i = 0; i < PrecursorResultRep.XYData.Yvalues.Length; i++)
            {
                if (!(apexIntensity < PrecursorResultRep.XYData.Yvalues[i])) continue;
                apexIntensity = PrecursorResultRep.XYData.Yvalues[i];
                _apexIndex = i;
            }
            
            foreach (var fragmentTargetResult in matchedResult.FragmentResultList)
            {

                var fragment = fragmentTargetResult.DatabaseFragmentTarget.Fragment;
                var residueNumber = fragment.ResidueNumber;
                var ion = new IonType(fragment.IonSymbol, fragment.ChargeState); //TODO isotope
                //  fragmentTargetResult.PeakQualityData.IsotopicProfile.
                if (residueNumber < 0 || FragmentXICSets.Length - 1 < residueNumber) continue;
                FragmentXICSets[residueNumber].AddIonXIC(ion, fragmentTargetResult.XYData.Yvalues);
                FragmentXICSets[residueNumber].AddIonXIC(precursorIon, PrecursorResultRep.XYData.Yvalues);
                FragmentSpectra[residueNumber].AddIntensity(ion, fragmentTargetResult.XYData.Yvalues[_apexIndex]);
                FragmentSpectra[residueNumber].AddIntensity(precursorIon, PrecursorResultRep.XYData.Yvalues[_apexIndex]);

                _matchedResidues.Add(residueNumber);
            }
        }



        private float GetScore() 
        {
            var score = 0f;
            for (var residueNumber = 1; residueNumber < Peptide.Count; residueNumber++)
            {
                float subScore;
                var fragmentXICSet = FragmentXICSets[residueNumber];
                var fragmentSpectrum = FragmentSpectra[residueNumber];
                var par = new FragmentSpectrumParameter(Peptide, residueNumber);

                if (fragmentXICSet.Count == 0 && NoXICScore.ContainsKey(par))
                {
                    subScore = NoXICScore[par];
                }
                else
                {
                    var specScorer = new FragmentSpectrumScorer(fragmentSpectrum, par);
                    var fragXICScorer = new FragmentXICSetScorer(fragmentXICSet, _apexIndex, par);
                    subScore = fragXICScorer.Score + specScorer.Score;
                    if (fragmentXICSet.Count == 0) NoXICScore[par] = subScore;
                    Console.WriteLine("Frag Score : " + residueNumber + "\t" + specScorer.Score + "\t" + fragXICScorer.RawScore + "\t" + 
                                      fragXICScorer.Score + "\t" + "\t" + fragmentXICSet.Count);
                }
                score += subScore; 
            }
            Console.WriteLine("Score : " + score);
            return score;                        
        }

    }
}
