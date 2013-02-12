using System.Collections.Generic;

namespace InformedProteomics.Backend.Scoring
{
    class FragmentSpectrumScorer
    {
        public FragmentSpectrum[] FragmentSpectra { get; private set; }
        public FragmentSpectrumParameter Parameter { get; private set; }
        public float Score { get; private set; }
        public List<IonType> UsedIons { get; private set; }
        public int BestSpectrumIndex { get; private set; }

        public FragmentSpectrumScorer(FragmentSpectrum[] fragmentSpectra, FragmentSpectrumParameter par)
        {
            FragmentSpectra = fragmentSpectra;
            Parameter = par;
            Score = GetScore();
        }

        private float GetScore()
        {
            var score = -50f;
            UsedIons = new List<IonType>();
            for (var i = 0; i < FragmentSpectra.Length; i++)
            {
                var t = FragmentSpectra[i];
                var ions = new List<IonType>();
                var subScore = 0f;
                for (var c = 1; c <= 5; c++)
                {
                    subScore += GetScoreOfOneSpectrum(t, c, ions, Parameter);
                }
                if (subScore > score && ions.Count!=0)
                {
                    score = subScore;
                    UsedIons = ions;
                    BestSpectrumIndex = i;
                }
            }
            //Console.WriteLine("***" + UsedIons.Count);
            if (score < -49f) score = -.1f;
            return score;
        }

        private static float GetScoreOfOneSpectrum(FragmentSpectrum spectrum, int charge, List<IonType> usedIons, FragmentSpectrumParameter par)
        {
            IonType bestSuffixIonType = null, bestPrefixIonType = null;
            IonType bestSuffixIonTypeWN = null, bestPrefixIonTypeWN = null;
                
            var suffixScore = -500f;
            var prefixScore = -500f;
            var prefixSuffixScore = -500f;
            var suffixNeutralLossScore = -500f;
            var prefixNeutralLossScore = -500f; // update later.... 

            foreach (var ion in spectrum.Keys)// get best suffix, prefix iontypes
            {
                if (ion.Charge != charge) continue;
                if (ion.NeutralLoss.Length != 0) continue;

                var s = SubScoreFactory.GetProductIonSpectrumScore(spectrum, ion, par);
                //Console.Write(s+"\t");
                if (ion.IsPrefix)
                {
                    if (s > prefixScore)
                    {
                        prefixScore = s;
                        bestPrefixIonType = ion;
                    }
                }
                else
                {
                    if (s > suffixScore)
                    {
                        suffixScore = s;
                        bestSuffixIonType = ion;
                    }
                }
            }


            foreach (var ion in spectrum.Keys)// get best suffix, prefix scores w/ neutral losses
            {
                if (ion.Charge != charge) continue;
                if (ion.NeutralLoss.Length == 0) continue;
                
                if (bestPrefixIonType == null && ion.IsPrefix)
                {
                    var s = SubScoreFactory.GetProductIonSpectrumScore(spectrum, ion, par);
                    if (s > prefixNeutralLossScore)
                    {
                        prefixNeutralLossScore = s;
                        bestPrefixIonTypeWN = ion;
                    }
                }

                if (bestPrefixIonType!=null && ion.Type.Equals(bestPrefixIonType.Type))
                {
                    var s = SubScoreFactory.GetProductIonSpectrumScore(spectrum, bestPrefixIonType, ion, par);
                    if (s > prefixNeutralLossScore)
                    {
                        prefixNeutralLossScore = s;
                        bestPrefixIonTypeWN = ion;
                    }
                }

                if (bestSuffixIonType == null && !ion.IsPrefix)
                {
                    var s = SubScoreFactory.GetProductIonSpectrumScore(spectrum, ion, par);
                    if (s > suffixNeutralLossScore)
                    {
                        suffixNeutralLossScore = s;
                        bestSuffixIonTypeWN = ion;
                    }
                }

                if (bestSuffixIonType != null && ion.Type.Equals(bestSuffixIonType.Type))
                {
                    var s = SubScoreFactory.GetProductIonSpectrumScore(spectrum, bestSuffixIonType, ion, par);
                    if (s > suffixNeutralLossScore)
                    {
                        suffixNeutralLossScore = s;
                        bestSuffixIonTypeWN = ion;
                    }
                }
            }

            if (bestPrefixIonType != null && bestSuffixIonType != null)
            {
                prefixSuffixScore = SubScoreFactory.GetProductIonSpectrumScore(spectrum, bestSuffixIonType, bestPrefixIonType, par);
            }

            if (spectrum.GetIntensity(bestPrefixIonType)>0) usedIons.Add(bestPrefixIonType);
            if (spectrum.GetIntensity(bestSuffixIonType)>0) usedIons.Add(bestSuffixIonType);
            if (spectrum.GetIntensity(bestPrefixIonTypeWN)>0) usedIons.Add(bestPrefixIonTypeWN);
            if (spectrum.GetIntensity(bestSuffixIonTypeWN)>0) usedIons.Add(bestSuffixIonTypeWN);

            
           // Console.WriteLine("**");
           // foreach(var io in usedIons) Console.Write(io.ToString() + " ");
            
           // Console.WriteLine();
            if (suffixScore < -49f) suffixScore = -.2f / charge; // not updated .. opt later
            if (prefixScore < -49f) prefixScore = -.2f / charge; // not updated .. opt later
            if (suffixNeutralLossScore < -49f) suffixNeutralLossScore = -.2f / charge; // not updated .. opt later
            if (prefixNeutralLossScore < -49f) prefixNeutralLossScore = -.2f / charge; // not updated .. opt later
            if (prefixSuffixScore < -49f) prefixSuffixScore = -.2f / charge; // not updated .. opt later

            return suffixScore + prefixScore + suffixNeutralLossScore + prefixNeutralLossScore + prefixSuffixScore;
        }
    }
}
