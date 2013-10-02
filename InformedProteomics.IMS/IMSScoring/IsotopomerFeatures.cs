using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.IMS.IMS;

namespace InformedProteomics.IMS.IMSScoring
{
    public class IsotopomerFeatures : List<Feature>
    {
        public int BinningMultiplyFactor { get; private set; }
        public float[] TheoreticalIsotopomerEnvelope { get; private set; }
        public int MaxIntensityIndex { get; private set; } // max index for _theoreticalIsotopomerEnvelope, not for this
        private IsotopomerFeatures(ImsDataCached imsData, Ion ion, Feature precursorFeature, bool isPrecurosr)
        {
          //  BinningMultiplyFactor = 6;
            if (ion.Charge >= 3) BinningMultiplyFactor = 1;
                //else if (ion.Charge == 1) BinningMultiplyFactor = 6;
                //else if (ion.Charge == 2) BinningMultiplyFactor = 2;
            else BinningMultiplyFactor = 2;

             var te = ion.Composition.GetIsotopomerEnvelop();
             TheoreticalIsotopomerEnvelope = new float[(te.Length + FeatureNode.OffsetFromMonoIsotope) * BinningMultiplyFactor];
             for (var i = FeatureNode.OffsetFromMonoIsotope; i < te.Length + FeatureNode.OffsetFromMonoIsotope; i++)
            {
                TheoreticalIsotopomerEnvelope[i * BinningMultiplyFactor] = te[i - FeatureNode.OffsetFromMonoIsotope];
            }

            MaxIntensityIndex = GetMaximumIndex(TheoreticalIsotopomerEnvelope);
            //
            var start = - FeatureNode.OffsetFromMonoIsotope*BinningMultiplyFactor;
            for (var i = start; i < TheoreticalIsotopomerEnvelope.Length + start; i++)
            {
                var mz = ion.GetIsotopeMz((float)i / BinningMultiplyFactor);
                if (isPrecurosr)
                {
                    if (mz > imsData.MaxPrecursorMz || mz < imsData.MinPrecursorMz) Add(null);
                    else Add(imsData.GetPrecursorFeature(mz, precursorFeature));
                }
                else
                {
                    if (mz > imsData.MaxFragmentMz || mz < imsData.MinFragmentMz) Add(null);
                    else Add(imsData.GetFramentFeature(mz, precursorFeature));
                    //var tt = imsData.GetFramentFeature(mz, precursorFeature);
                   // if(mz == 361.198813765) 
                    //    Console.WriteLine(mz + " " + (tt == null ? "NULL" : tt.GetHighestPoint().ToString()));
                }
                
            }
        }

        static public IsotopomerFeatures GetFramentIsotopomerFeatures(ImsDataCached imsData, Composition cutComposition, IonType fragmentIonClassBase, Feature precursorFeature)
        {
            return new IsotopomerFeatures(imsData, fragmentIonClassBase.GetIon(cutComposition), precursorFeature, false);
        }

        static public IsotopomerFeatures GetPrecursorIsotopomerFeatures(ImsDataCached imsData, Ion precursorIon, Feature precursorFeature)
        {
            return new IsotopomerFeatures(imsData, precursorIon, precursorFeature, true);
        }

        public Feature GetMostAbundantFeature()
        {
            var index = 0;
            var max = 0.0;
            for (var i = 0; i < Count;i++ )
            {
                var feature = this[i];
                if (feature != null && max < feature.IntensityMax)
                {
                    max = feature.IntensityMax;
                    index = i;
                }
            }
            return this[index];
        }

        public Feature GetTheoreticallyMostAbundantFeature()
        {
            return this[MaxIntensityIndex];
        }

        static private int GetMaximumIndex(float[] e)
        {
            var maxIndex = 0;
            var max = float.NegativeInfinity;
            for (var i = 0; i < e.Length; i++)
            {
                if (!(e[i] > max)) continue;
                max = e[i];
                maxIndex = i;
            }
            return maxIndex;
        }
    }
}
