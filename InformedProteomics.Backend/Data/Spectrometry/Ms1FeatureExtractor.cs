using System;
using System.CodeDom.Compiler;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Windows.Forms;
using InformedProteomics.Backend.Data.Biology;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Utils;
using QuickGraph;

namespace InformedProteomics.Backend.Data.Spectrometry
{
    public class Ms1Feature
    {
        public double MonoIosotopeMass { get; private set; }

        //scan, monomass, mz, charge, abundance, score
        internal Ms1Feature(ChargeLcScanCluster cluster, double monoMass, int[] ms1ScanNums, int minCharge)
        {
            _cluster = cluster;
            MonoIosotopeMass = monoMass;
            _minCharge = minCharge;
            _ms1ScanNums = ms1ScanNums;
        }
        
        public double Abundance
        {
            get { return _cluster.TotalIntensity; }
        }

        public IEnumerable<Tuple<int, double>> GetScanNumIntensityPair()
        {
            var ret = new List<Tuple<int, double>>();

            for (var col = _cluster.MinCol; col <= _cluster.MaxCol; col++)
            {
                ret.Add(new Tuple<int, double>(_ms1ScanNums[col], _cluster.IntensityAlongCol[col - _cluster.MinCol]));
            }
            return ret;
        }

        public double Score
        {
            get { return _cluster.Score; }
        }

        public int Charge
        {
            get
            {
                var cell = _cluster.GetHighestIntensityCell();
                return cell.Row + _minCharge;
            }
        }

        public int ScanNum
        {
            get
            {
                var cell = _cluster.GetHighestIntensityCell();
                return _ms1ScanNums[cell.Col];
            }
        }

        public double Mz
        {
            get
            {
                var isoEnv = Averagine.GetIsotopomerEnvelope(MonoIosotopeMass);
                return Ion.GetIsotopeMz(MonoIosotopeMass, Charge, isoEnv.MostAbundantIsotopeIndex);                
            }
        }

        public int MaxScanNum
        {
            get { return _ms1ScanNums[_cluster.MaxCol]; }
        }

        public int MinScanNum
        {
            get { return _ms1ScanNums[_cluster.MinCol]; }
        }

        public int MinCharge
        {
            get { return _minCharge + _cluster.MinRow; }
        }

        public int MaxCharge
        {
            get { return _minCharge + _cluster.MaxRow; }
        }


        private readonly ChargeLcScanCluster _cluster;
        private readonly int _minCharge;
        private readonly int[] _ms1ScanNums;
    }
    
    public class Ms1FeatureExtractor
    {
        private readonly ChargeLcScanMatrix _csm;
        private readonly double _minSearchMass;
        private readonly double _maxSearchMass;
        private readonly int _minSearchCharge;
        private readonly int _maxSearchCharge;
        private readonly int _minSearchMassBin;
        private readonly int _maxSearchMassBin;
        private readonly MzComparerWithBinning _comparer;
        private readonly LcMsRun _run;


        private readonly List<ChargeLcScanCluster>[] _massBinToClusterMap;

        public Ms1FeatureExtractor(LcMsRun run, double searchMinMass = 3000, double searchMaxMass = 50000, int searchMinCharge = 2, int searchMaxCharge = 50)
        {
            _minSearchMass = searchMinMass;
            _maxSearchMass = searchMaxMass;
            _minSearchCharge = searchMinCharge;
            _maxSearchCharge = searchMaxCharge;
            _run = run;
            _csm = new ChargeLcScanMatrix(run, 27, _minSearchCharge, _maxSearchCharge);
            _comparer = _csm.GetMzComparerWithBinning();

            _minSearchMassBin = _comparer.GetBinNumber(_minSearchMass);
            _maxSearchMassBin = _comparer.GetBinNumber(_maxSearchMass);
            _massBinToClusterMap = new List<ChargeLcScanCluster>[_maxSearchMassBin - _minSearchMassBin + 1];

            for(var i = 0; i < _massBinToClusterMap.Length; i++)
                _massBinToClusterMap[i] = new List<ChargeLcScanCluster>();
        }

        public void Run()
        {
            //var ms1ScanNumbers = _run.GetMs1ScanVector();

            for (var binNum = _minSearchMassBin; binNum <= _maxSearchMassBin; binNum++)
            {
                if ((binNum - _minSearchMassBin)%1000 == 0)
                {
                    Console.WriteLine("Processed {0} mass bins for total {1} bins", binNum - _minSearchMassBin, _maxSearchMassBin - _minSearchMassBin + 1);
                }
                
                var monoMass = _comparer.GetMzAverage(binNum);
                //_massBinToClusterMap[binNum - _minSearchMassBin].AddRange(_csm.GetMs1Features(monoMass));


                var neighborMassBins = new List<int>();

                for (var i = -2; i <= 2; i++)
                {
                    if (i == 0) continue;
                    var neighborIsotopebinNum = _comparer.GetBinNumber(monoMass + i);
                    if (neighborIsotopebinNum >= _minSearchMassBin && neighborIsotopebinNum <= _maxSearchMassBin) neighborMassBins.Add(neighborIsotopebinNum);
                }
                
                foreach (var cluster in _csm.GetMs1Features(monoMass))
                {
                    // Before adding it, Check if there are existing neighbor that can be merged
                    // search mass range = +-2 [Da]
                    cluster.Active = true;
                    var foundNeighbor = false;
                    
                    foreach (var neighborBin in neighborMassBins)
                    {
                        foreach (var neighborCluster in _massBinToClusterMap[neighborBin - _minSearchMassBin])
                        {
                            if (neighborCluster.Overlaps(cluster))
                            {
                                if (neighborCluster.Active)
                                {
                                    if (neighborCluster.Score+neighborCluster.Score2 > cluster.Score+cluster.Score2)
                                    {
                                        cluster.Active = false;
                                    }
                                    else
                                    {
                                        neighborCluster.Active = false;
                                        cluster.Active = true;
                                    }
                                }
                                else
                                {
                                    cluster.Active = false;
                                }

                                foundNeighbor = true;
                                break;
                            }
                        }

                        if (foundNeighbor) break;
                    }

                    _massBinToClusterMap[binNum - _minSearchMassBin].Add(cluster);
                    //scan, monomass, mz, charge, abundance, score
                    //Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", scanNum, monoMass, mostAbundantIsotopeMz, charge, cluster.TotalIntensity, cluster.GetScore());
                }
            }
        }

        public IEnumerable<Ms1Feature> GetMs1Features()
        {
            var comparer = _csm.GetMzComparerWithBinning();
            var ms1ScanNums = _run.GetMs1ScanVector();
            var features = new List<Ms1Feature>();

            for (var i = 0; i < _massBinToClusterMap.Length; i++)
            {
                var monoMass = comparer.GetMzAverage(i + _minSearchMassBin);

                foreach (var cluster in _massBinToClusterMap[i])
                {
                    if (cluster.Active == false) continue;
                    features.Add(new Ms1Feature(cluster, monoMass, ms1ScanNums, _minSearchCharge));
                }

            }
            return features.OrderBy(x => x.ScanNum);
        }
        
    }
}
