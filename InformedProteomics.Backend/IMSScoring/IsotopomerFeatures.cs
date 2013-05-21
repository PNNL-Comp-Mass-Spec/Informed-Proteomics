using System;
using System.Collections.Generic;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.IMS;
using Feature = InformedProteomics.Backend.IMS.Feature;

namespace InformedProteomics.Backend.IMSScoring
{
    public class IsotopomerFeatures : List<Feature>
    {
        private readonly float[] _theoreticalIsotopomerEnvelope;
        private readonly int _maxIntensityIndex;// max index for _theoreticalIsotopomerEnvelope, not for this
        public int SupportSize { get; private set; }
        private IsotopomerFeatures(ImsDataCached imsData, Ion ion, Feature precursorFeature, bool isPrecurosr)
        {
            _theoreticalIsotopomerEnvelope = ion.Composition.GetApproximatedIsotopomerEnvelop(FeatureNode.NumMaxSupport);
            _maxIntensityIndex = GetMaximumIndex(_theoreticalIsotopomerEnvelope);
            SupportSize = _theoreticalIsotopomerEnvelope.Length + FeatureNode.NumMinusIsotope;

            for (var i = -FeatureNode.NumMinusIsotope; i < _theoreticalIsotopomerEnvelope.Length - FeatureNode.NumMinusIsotope; i++)
            {
                var mz = ion.GetIsotopeMz(i + _maxIntensityIndex);
                if (isPrecurosr)
                {
                    if (mz > imsData.MaxPrecursorMz || mz < imsData.MinPrecursorMz) Add(null);
                    else Add(imsData.GetPrecursorFeature(mz, precursorFeature));
                }
                else
                {
                    if (mz > imsData.MaxFragmentMz || mz < imsData.MinFragmentMz) Add(null);
                    else Add(imsData.GetFramentFeature(mz, precursorFeature));
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

        public Feature GetNthFeatureFromTheoreticallyMostIntenseFeature(int n)
        {
            return this[Math.Min(Count - 1, n + _maxIntensityIndex + FeatureNode.NumMinusIsotope)]; // this[_maxIntensityIndex] corresponds to the max intensity istope - FeatureNode.NumMinusIsotope isotope
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

        public float GetTheoreticalIntensityOfNthFeature(int n)
        {
            return n + _maxIntensityIndex < 0 ? 0f : _theoreticalIsotopomerEnvelope[Math.Min(_theoreticalIsotopomerEnvelope.Length - 1, n + _maxIntensityIndex)];
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
