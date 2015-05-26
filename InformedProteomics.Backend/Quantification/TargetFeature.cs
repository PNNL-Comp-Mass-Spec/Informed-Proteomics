using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Utils;


namespace InformedProteomics.Backend.Quantification
{
    public class TargetFeature
    {
        public TargetFeature(double mass, int charge, int scan)
        {
            ScanNum = scan;
            Mass = mass;
            Charge = charge;
            
            MinCharge = charge;
            MaxCharge = charge;
            MinScanNum = scan;
            MaxScanNum = scan;
            _quantifiedMs1Feature = null;
        }        
        
        public TargetFeature(int datasetId, int id, double mass, int charge, int scan) : this(mass, charge, scan)
        {
            DataSetId = datasetId;
            Id = id;
        }
        
        public int ScanNum { get; set; }
        public int Charge { get; set; }
        public double Mass { get; set; }
        public int DataSetId { get; set; }
        public int Id { get; set; }
        public string Sequence { get; set; }
        public string Modifications { get; set; }
        public string SequenceText { get; set; }

        public string ProteinName { get; set; }
        public string Composition { get; set; }

        public int MinCharge { get; set; }
        public int MaxCharge { get; set; }
        public int ChargeLength { get { return (MaxCharge == 0) ? 0 : MaxCharge - MinCharge + 1; } }

        public int MinScanNum { get; set; }
        public int MaxScanNum { get; set; }
        public int ScanLength { get { return (MaxScanNum == 0) ? 0 : MaxScanNum - MinScanNum + 1; } }

        public Ms1Feature LinkedMs1Feature { get { return _quantifiedMs1Feature; }}

        public List<string> Rows = new List<string>();
        
        public void SetMs1Feature(Ms1Feature ms1Feature)
        {
            _quantifiedMs1Feature = ms1Feature;
        }

        public double QuantifiedAbundance
        {
            get { return _quantifiedMs1Feature == null ? 0d : _quantifiedMs1Feature.Abundance; }
        }

        public static ICollection<TargetFeature> CollectUniqueTargets(string pbfFilePath, string idFilePath, double netTh = 0.005)
        {
            var tsvReader = new TsvFileParser(idFilePath);
            var run = PbfLcMsRun.GetLcMsRun(pbfFilePath);
            var targetSet = new Dictionary<string, TargetFeature>();

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
                var seqKeyPair = SetModifications(seq, mod);
                var good = IsGoodTarget(spec, seqKeyPair.Item1);
                var row = tsvReader.GetRows()[i];

                if (!good) continue;
                
                TargetFeature tf;

                if (targetSet.TryGetValue(comp, out tf))
                {
                    var minNet = run.GetElutionTime(tf.MinScanNum)/run.GetElutionTime(run.MaxLcScan);
                    var maxNet = run.GetElutionTime(tf.MaxScanNum)/run.GetElutionTime(run.MaxLcScan);
                    var net = run.GetElutionTime(scan) / run.GetElutionTime(run.MaxLcScan);
                    var netDiff = (net >= minNet && net <= maxNet) ? 0d : Math.Min(Math.Abs(minNet - net), Math.Abs(maxNet - net));
                    
                    if (netDiff > netTh)
                    {
                        //Console.WriteLine("{0}\t{1}", scan, seq);
                    }
                    else
                    {
                        tf.MaxScanNum = Math.Max(tf.MaxScanNum, scan);
                        tf.MinScanNum = Math.Min(tf.MinScanNum, scan);
                        tf.MaxCharge = Math.Max(tf.MaxCharge, charge);
                        tf.MinCharge = Math.Min(tf.MinCharge, charge);
                        tf.Rows.Add(row);
                    }
                }
                else
                {
                    tf = new TargetFeature(mass, charge, scan)
                    {
                        MinScanNum = scan,
                        MaxScanNum = scan,
                        MinCharge = charge,
                        MaxCharge = charge,
                        Sequence = seq,
                        SequenceText = seqKeyPair.Item2,
                        Modifications = mod,
                        Composition = comp,
                    };
                    tf.Rows.Add(row);

                    targetSet.Add(comp, tf);
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

            if ((double) nObservedPrefixIonPeaks/nObservedIonPeaks > 0.9 ||
                (double) nObservedsuffixIonPeaks/nObservedIonPeaks > 0.9) return false;

            return true;
        }

        private static Tuple<Sequence, string> SetModifications(string cleanSequence, string modifications)
        {
            // Build Sequence AminoAcid list
            var sequence = new Sequence(cleanSequence, new AminoAcidSet());
            var sequenceText = cleanSequence;
            var parsedModifications = ParseModifications(modifications);

            // Add modifications to sequence
            parsedModifications.Sort(new CompareModByHighestPosition());   // sort in reverse order for insertion
            foreach (var mod in parsedModifications)
            {
                var pos = mod.Item1;
                if (pos > 0)
                {
                    pos--;
                }

                var modLabel = string.Format("[{0}]", mod.Item2.Name);
                sequenceText = sequenceText.Insert(mod.Item1, modLabel);
                var aa = sequence[pos];
                var modaa = new ModifiedAminoAcid(aa, mod.Item2);
                sequence[pos] = modaa;
            }

            return new Tuple<Sequence, string>(new Sequence(sequence), sequenceText);
        }

        private static List<Tuple<int, Modification>> ParseModifications(string modifications)
        {
            var mods = modifications.Split(',');
            var parsedMods = new List<Tuple<int, Modification>>();
            if (mods.Length < 1 || mods[0] == string.Empty)
            {
                return parsedMods;
            }

            foreach (var modParts in mods.Select(mod => mod.Split(' ')))
            {
                if (modParts.Length < 0)
                {
                    throw new FormatException("Unknown Modification");
                }

                var modName = modParts[0];
                var modPos = Convert.ToInt32(modParts[1]);
                var modification = Modification.Get(modName);
                parsedMods.Add(new Tuple<int, Modification>(modPos, modification));
                if (modification == null)
                {
                    throw new Exception(string.Format("Found an unrecognized modification: {0}", modName));
                }
            }

            return parsedMods;
        }

        private Ms1Feature _quantifiedMs1Feature;


        internal class CompareModByHighestPosition : IComparer<Tuple<int, Modification>>
        {
            public int Compare(Tuple<int, Modification> x, Tuple<int, Modification> y)
            {
                return y.Item1.CompareTo(x.Item1);
            }
        }     
   
    }

}
