using System.Collections.Generic;
using InformedProteomics.Backend.MassSpecData;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsPeakCluster : LcMsFeature
    {
        public LcMsPeakCluster(LcMsRun run, ObservedIsotopeEnvelope observedEnvelope)
            : base(observedEnvelope.MonoMass, observedEnvelope.Charge, observedEnvelope.RepresentativePeak.Mz, observedEnvelope.ScanNum, observedEnvelope.Abundance)
        {
            Envelopes = new List<ObservedIsotopeEnvelope>();
            Run = run;
            AddObservedEnvelope(observedEnvelope);
        }

        public double[] EnvelopeDistanceScoreAcrossTime { get; internal set; }
        public double[] EnvelopeDistanceScoreAcrossCharge { get; internal set; }

        public void SetAbundance(double abu)
        {
            Abundance = abu;
        }

        public void AddObservedEnvelope(ObservedIsotopeEnvelope envelope)
        {
            if (MaxScanNum < 0 || envelope.ScanNum > MaxScanNum)
            {
                MaxScanNum = envelope.ScanNum;
            }
            if (MinScanNum < 0 || envelope.ScanNum < MinScanNum)
            {
                MinScanNum = envelope.ScanNum;
            }
            if (MaxCharge < 0 || envelope.Charge > MaxCharge)
            {
                MaxCharge = envelope.Charge;
            }
            if (MinCharge < 0 || envelope.Charge < MinCharge)
            {
                MinCharge = envelope.Charge;
            }
            Envelopes.Add(envelope);
        }

        public void Clear()
        {
            MaxScanNum = -1;
            MinScanNum = -1;
            MinCharge = -1;
            MaxCharge = -1;
            Envelopes.Clear();
        }
        
        public void ExpandElutionRange()
        {
            // considering DDA instrument
            if (Run.MaxMsLevel > 1 && NetLength < 0.01)
            {
                var ms1ScanNums = Run.GetMs1ScanVector();
                var ms1ScanNumToIndex = Run.GetMs1ScanNumToIndex();

                var minCol = ms1ScanNumToIndex[MinScanNum];
                var maxCol = ms1ScanNumToIndex[MaxScanNum];

                for (var i = minCol - 1; i >= 0; i--)
                {
                    if (MinScanNum - ms1ScanNums[i] > (minCol - i))
                    {
                        MinScanNum = ms1ScanNums[i];
                        break;
                    }
                }
                for (var i = maxCol + 1; i < ms1ScanNums.Length; i++)
                {
                    if (ms1ScanNums[i] - MaxScanNum > (i - maxCol))
                    {
                        MaxScanNum = ms1ScanNums[i];
                        break;
                    }
                }
            }
        }
        
        public List<ObservedIsotopeEnvelope> Envelopes { get; private set; }
    }
}
