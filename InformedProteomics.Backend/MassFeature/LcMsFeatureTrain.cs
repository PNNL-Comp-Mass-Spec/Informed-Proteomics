using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.Backend.MassFeature
{
    public class LcMsFeatureTrain
    {
        public static ICollection<IdentifiedProtein> CollectTrainSet(string pbfFilePath, string idFilePath, double netTh = 0.005)
        {
            Modification.RegisterAndGetModification(Modification.Cysteinyl.Name, Modification.Cysteinyl.Composition);
            
            var tsvReader = new TsvFileParser(idFilePath);
            var run = PbfLcMsRun.GetLcMsRun(pbfFilePath);
            var targetSet = new Dictionary<string, IdentifiedProtein>();

            var scanNums = new int[tsvReader.NumData];
            for (var i = 0; i < tsvReader.NumData; i++) scanNums[i] = int.Parse(tsvReader.GetData("Scan")[i]);
            var index = Enumerable.Range(0, scanNums.Length).ToArray();

            Array.Sort(scanNums, index);

            //for (var i = 0; i < tsvReader.NumData; i++)
            foreach (var i in index)
            {
                var qv = double.Parse(tsvReader.GetData("QValue")[i]);
                if (qv > 0.01) continue;

                var scan = int.Parse(tsvReader.GetData("Scan")[i]);
                var charge = int.Parse(tsvReader.GetData("Charge")[i]);
                var mass = double.Parse(tsvReader.GetData("Mass")[i]);
                var seq = tsvReader.GetData("Sequence")[i];
                var mod = tsvReader.GetData("Modifications")[i];
                var comp = tsvReader.GetData("Composition")[i];

                var spec = run.GetSpectrum(scan) as ProductSpectrum;
                //var seqKeyPair = SetModifications(seq, mod);

                var identifiedProtein = new IdentifiedProtein(mass, charge, scan, seq, mod, 0, i)
                {
                    Composition = comp,
                };

                var good = IsGoodTarget(spec, identifiedProtein.GetSequence());
                var row = tsvReader.GetRows()[i];

                if (!good) continue;

                IdentifiedProtein ip;

                if (targetSet.TryGetValue(comp, out ip))
                {
                    var minNet = run.GetElutionTime(ip.MinScanNum) / run.GetElutionTime(run.MaxLcScan);
                    var maxNet = run.GetElutionTime(ip.MaxScanNum) / run.GetElutionTime(run.MaxLcScan);
                    var net = run.GetElutionTime(scan) / run.GetElutionTime(run.MaxLcScan);
                    var netDiff = (net >= minNet && net <= maxNet) ? 0d : Math.Min(Math.Abs(minNet - net), Math.Abs(maxNet - net));

                    if (netDiff > netTh)
                    {
                        //Console.WriteLine("{0}\t{1}", scan, seq);
                    }
                    else
                    {
                        ip.MaxScanNum = Math.Max(ip.MaxScanNum, scan);
                        ip.MinScanNum = Math.Min(ip.MinScanNum, scan);
                        ip.MaxCharge = Math.Max(ip.MaxCharge, charge);
                        ip.MinCharge = Math.Min(ip.MinCharge, charge);
                        ip.Rows.Add(row);
                    }
                }
                else
                {
                    identifiedProtein.Rows.Add(row);
                    targetSet.Add(comp, identifiedProtein);
                }
            }

            return targetSet.Values;
        }

        private static bool IsGoodTarget(ProductSpectrum ms2Spec, Sequence sequence)
        {
            //var counter = new MatchedPeakCounter(ms2Spec, new Tolerance(10), 1, 10);
            var BaseIonTypesCID = new[] { BaseIonType.B, BaseIonType.Y };
            var BaseIonTypesETD = new[] { BaseIonType.C, BaseIonType.Z };
            var tolerance = new Tolerance(10);
            var minCharge = 1;
            var maxCharge = 10;
            var baseIonTypes = ms2Spec.ActivationMethod != ActivationMethod.ETD ? BaseIonTypesCID : BaseIonTypesETD;
            var cleavages = sequence.GetInternalCleavages();
            const double RelativeIsotopeIntensityThreshold = 0.8d;


            var nTheoreticalIonPeaks = 0;
            var nObservedIonPeaks = 0;
            var nObservedPrefixIonPeaks = 0;
            var nObservedsuffixIonPeaks = 0;

            int index = 0; // cleavage index 
            foreach (var c in cleavages)
            {
                foreach (var baseIonType in baseIonTypes)
                {
                    var fragmentComposition = baseIonType.IsPrefix
                                  ? c.PrefixComposition + baseIonType.OffsetComposition
                                  : c.SuffixComposition + baseIonType.OffsetComposition;

                    for (var charge = minCharge; charge <= maxCharge; charge++)
                    {
                        var ion = new Ion(fragmentComposition, charge);
                        int baseIsotopePeakIndex;
                        int nIsotopes;
                        int nMatchedIsotopes;
                        //if (FindIon(ion, tolerance, RelativeIsotopeIntensityThreshold, out baseIsotopePeakIndex, out nIsotopes, out nMatchedIsotopes))
                        if (ms2Spec.ContainsIon(ion, tolerance, RelativeIsotopeIntensityThreshold))
                        {
                            if (baseIonType.IsPrefix) nObservedPrefixIonPeaks++;
                            else nObservedsuffixIonPeaks++;
                            nObservedIonPeaks++;
                        }
                        //_nObservedIonPeaks += nMatchedIsotopes;
                        //_nTheoreticalIonPeaks += nIsotopes;
                        nTheoreticalIonPeaks++;
                    }
                }
                index++;
            }

            if ((double)nObservedPrefixIonPeaks / nObservedIonPeaks > 0.9 ||
                (double)nObservedsuffixIonPeaks / nObservedIonPeaks > 0.9) return false;

            return true;
        }

        


    }
}
