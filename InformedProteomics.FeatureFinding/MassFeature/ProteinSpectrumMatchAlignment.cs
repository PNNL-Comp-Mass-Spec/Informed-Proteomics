using System.Collections.Generic;
using System.Linq;

namespace InformedProteomics.FeatureFinding.MassFeature
{
    public class ProteinSpectrumMatchAlignment
    {
        public List<ProteinSpectrumMatchSet> GroupingByPrsm(int dataid, IEnumerable<ProteinSpectrumMatch> matches, INodeComparer<ProteinSpectrumMatch> prsmComparer)
        {
            var prsmSet = new NodeSet<ProteinSpectrumMatch>(){};
            prsmSet.AddRange(matches);
            var groupList = prsmSet.ConnnectedComponents(prsmComparer);
            return groupList.Select(@group => new ProteinSpectrumMatchSet(dataid, @group)).ToList();
        }

        public ProteinSpectrumMatchSet[][] GroupAcrossRuns(List<ProteinSpectrumMatchSet>[] prsmGroup, INodeComparer<ProteinSpectrumMatchSet> prsmGroupComparer)
        {
            var nDataset = prsmGroup.Length;
            var prsmSet = new NodeSet<ProteinSpectrumMatchSet>() { };

            for (var i = 0; i < nDataset; i++)
            {
                var groupedPrsms = prsmGroup[i];
                if (groupedPrsms == null) continue;
                prsmSet.AddRange(groupedPrsms);
            }

            var alignedPrsms = prsmSet.ConnnectedComponents(prsmGroupComparer);
            var alignedResult = new ProteinSpectrumMatchSet[alignedPrsms.Count][];
            for (var i = 0; i < alignedResult.Length; i++) alignedResult[i] = new ProteinSpectrumMatchSet[nDataset];

            for(var i = 0; i < alignedPrsms.Count; i++)
            {
                foreach (var set in alignedPrsms[i])
                {
                    if (alignedResult[i][set.DataId] != null)
                    {
                        alignedResult[i][set.DataId].Merge(set);
                        //Console.WriteLine("[{4}] {0}-{1}...{2}-{3}", set.MinScanNum, set.MaxScanNum, alignedResult[i][set.DataId].MinScanNum, alignedResult[i][set.DataId].MaxScanNum, set.DataId);
                    }
                    else
                    {
                        alignedResult[i][set.DataId] = set;
                    }
                }
            }
            return alignedResult;
        }

        /*
        public LcMsFeature[][] GetLcMsFeatures(ProteinSpectrumMatcheSet[][] alignedPrsmSet, LcMsRun[] runs)
        {
            var nDataset = runs.Length;

            var minMaxRegions = new List<Tuple<double, int, double, double>>();
            for (var j = 0; j < alignedPrsmSet.Length; j++)
            {
                minMaxRegions.Add(GetMinMaxScanRegion(runs, alignedPrsmSet[j]));
            }

            for (var i = 0; i < nDataset; i++)
            {
                for (var j = 0; j < alignedPrsmSet.Length; j++)
                {
                    var region = minMaxRegions[j];

                    //scanregions
                }
            }

            Console.WriteLine(alignedPrsmSet.Length);
        }

        private Tuple<double, int, double, double> GetMinMaxScanRegion(LcMsRun[] runs, ProteinSpectrumMatcheSet[] matches)
        {
            var minElutionTime = double.MaxValue;
            var maxElutionTime = 0d;
            var massList = new List<double>();
            var chargeList = new List<int>();

            for (var i = 0; i < matches.Length; i++)
            {
                if (matches[i] == null) continue;
                var match = matches[i];

                minElutionTime = Math.Min(minElutionTime, runs[i].GetElutionTime(match.MinScanNum));
                maxElutionTime = Math.Max(maxElutionTime, runs[i].GetElutionTime(match.MaxScanNum));

                massList.AddRange(match.Select(mt=>mt.Mass));
                chargeList.AddRange(match.Select(mt => mt.Charge));
            }
            massList.Sort();
            var mass = massList[(int)(massList.Count * 0.5)];
            var charge = chargeList[(int)(chargeList.Count * 0.5)];

            return new Tuple<double, int, double, double>(mass, charge, minElutionTime, maxElutionTime);
        }
        */
    }
}
