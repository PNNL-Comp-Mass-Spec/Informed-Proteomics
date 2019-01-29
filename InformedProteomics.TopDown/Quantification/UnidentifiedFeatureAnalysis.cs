using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using InformedProteomics.Backend.Data.Spectrometry;
using InformedProteomics.Backend.Database;
using InformedProteomics.Backend.MassSpecData;
using InformedProteomics.TopDown.SequenceTag;

namespace InformedProteomics.TopDown.Quantification
{
    public class UnidentifiedFeatureAnalysis
    {
        /*
         * This is a class that Jung Kap asked me to complete to analyze unidentified features from a cross tab file
         * This class checks to see for each feature how many spectrum contained it, from the spectrum how many tags were generated, and then
         * how many of those tags exist in our data base.
         *
         */
        public UnidentifiedFeatureAnalysis(string[] rawFiles, string crossTabFile, string databaseFile)
        {
            _rawFiles = rawFiles;
            _crossTabFile = crossTabFile;
            _databaseFile = databaseFile;
            _searchableDB = new SearchableDatabase(new FastaDatabase(_databaseFile));
            _filteredFeatures = new List<Tuple<int,double, double, double, double, double>>();
        }

        /*
         * This filters the features to get the unidentified features that have a p-val > .01
         * */
        public void FilterFeatures(string outputFolder)
        {
            var filteredCount = 0;
            string featureLine;
            var tsv = new StringBuilder();
            var featureFile = new StreamReader(_crossTabFile);
            featureFile.ReadLine();
            while (((featureLine = featureFile.ReadLine()) != null))
            {
                var featureElements = featureLine.Split('\t');
                var sequence = featureElements[1];
                var pVal = double.Parse(featureElements[9]);
                var logFc = double.Parse(featureElements[10]);
                if (!(string.IsNullOrEmpty(sequence) && pVal > .01)) continue;
                var mass = featureElements[12];
                var minElution = featureElements[13];
                var maxElution = featureElements[14];
                filteredCount++;
                tsv.Append(string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n", filteredCount,mass, minElution, maxElution, pVal,logFc));
            }

            _filteredFile = outputFolder + "filteredFeatures.tsv";
            File.WriteAllText(_filteredFile, tsv.ToString());

            _spectrumMatchesMatrix = new int[filteredCount][];
            _tagsGeneratedMatrix = new int[filteredCount][];
            _dataBaseHitMatrix = new int[filteredCount][];
            _identifiedFeatures = new Tuple<int, double, int>[filteredCount][];
            for (var i = 0; i < _identifiedFeatures.Length; i++)
            {
                _identifiedFeatures[i] = new Tuple<int, double, int>[_rawFiles.Length];
            }

            Console.WriteLine("# of features left after filter: {0}", filteredCount);
        }

        /**
         * This prints two files:
         *
         * 1) a general file that looks at each feature and gives the spectrum match count, tag count and database hit count
         * 2) a more specific file that matches a feature to a scan # and then the mz and charge value being looked at.
         *
         * Format is specific to data Jung gave me
         */
        public void PrintAnalysisToFile(string outputFolder)
        {
            var outputFile = outputFolder + "analysis.tsv";
            var tsv = new StringBuilder();
            const string headerString = "32A\t32B\t32C\t32D\t32E\t32F\t32G\t33A\t33B\t33C\t33D\t33E\t33F\t33G";
            tsv.Append(string.Format("id\tmass\tminElutionTime\tmaxElutionTime\tpVal\tfold-change\t{0}\t{1}\t{2}\n",headerString,headerString,headerString));
            for (var i = 0; i < _filteredFeatures.Count; i++)
            {
                var line = string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t", _filteredFeatures[i].Item1,_filteredFeatures[i].Item2, _filteredFeatures[i].Item3, _filteredFeatures[i].Item4, _filteredFeatures[i].Item5,_filteredFeatures[i].Item6);
                _spectrumMatchesMatrix[i].ToList().ForEach(x => line += x.ToString() + "\t");
                _tagsGeneratedMatrix[i].ToList().ForEach(x => line += x.ToString() + "\t");
                _dataBaseHitMatrix[i].ToList().ForEach(x => line += x.ToString() + "\t");
                line += "\n";
                tsv.Append(line);
            }
            File.WriteAllText(outputFile, tsv.ToString());

            const string head = "Scan#\tMz\tCharge";
            tsv.Clear();
            tsv.Append("feature\tMonoMass\t32A\t\t\t32B\t\t\t32C\t\t\t32D\t\t\t32E\t\t\t32F\t\t\t32G\t\t\t33A\t\t\t33B\t\t\t33C\t\t\t33D\t\t\t33E\t\t\t33F\t\t\t33G\n");
            tsv.Append(string.Format("\t\t\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n",head,head,head,head,head,head,head,head,head,head,head,head,head,head));
            for (var i = 0; i < _identifiedFeatures.Length; i++)
            {
                var line = string.Format("{0}\t{1}\t", i, _filteredFeatures[i].Item2);
                _identifiedFeatures[i].ToList().ForEach(x => line += x.Item1 + "\t" + x.Item2 +"\t" + x.Item3 + "\t");
                line += "\n";
                tsv.Append(line);
            }
            File.WriteAllText(outputFolder + "featureDetails.tsv",tsv.ToString());
        }

        public void DoAnalysis()
        {
            InitializeMatrix(_spectrumMatchesMatrix);
            InitializeMatrix(_tagsGeneratedMatrix);
            InitializeMatrix(_dataBaseHitMatrix);
            GetFilteredFeatures(_filteredFile);

            for (var i = 0; i < _rawFiles.Length; i++)
            {
                Console.WriteLine("Processing File {0}...............", i);
                var run = PbfLcMsRun.GetLcMsRun(_rawFiles[i]);
                var ms2List = run.GetScanNumbers(2);
                Console.WriteLine("# of scans {0}",ms2List.Count);
                for (var j = 0; j < _filteredFeatures.Count; j++)
                {
                    var matchedSpecList = GetMatchedSpectra(run, ms2List, _filteredFeatures[j],i);
                    _spectrumMatchesMatrix[j][i] = matchedSpecList.Count;

                    var tags = GetTags(matchedSpecList);
                    _tagsGeneratedMatrix[j][i] = tags.Count;

                    var hitCount = TagsInDatabase(tags);
                    _dataBaseHitMatrix[j][i] = hitCount;
                }
            }
        }

        private List<ProductSpectrum> GetMatchedSpectra(ILcMsRun run, IList<int> ms2List ,Tuple<int,double, double, double, double,double> feature, int fileIndex)
        {
            var spectrumList = new List<ProductSpectrum>();
            var mass = feature.Item2;
            var minElution = feature.Item3;
            var maxElution = feature.Item4;
            var featureId = feature.Item1 - 1;
            var det = new Tuple<int, double, int>(-1, -1.0, -1);
            foreach (var scanNum in ms2List)
            {
                var scanElutionTime = run.GetElutionTime(scanNum);
                if (scanElutionTime < minElution || scanElutionTime > maxElution) continue;

                if (!(run.GetSpectrum(scanNum) is ProductSpectrum spectrum))
                {
                    Console.WriteLine("Unable to retrieve the spectrum for scan " + scanNum);
                    continue;
                }

                var window = spectrum.IsolationWindow;
                var minMz = window.MinMz - .5;
                var maxMz = window.MaxMz + .5;
                var mzTable = GetFeatureMassTable(mass);

                for (var j = 0; j < mzTable.Length; j++)
                {
                    var mz = mzTable[j];
                    if (mz < minMz || mz > maxMz) continue;
                    spectrumList.Add(spectrum);
                    if(!(det.Item1 > -1)) det = new Tuple<int, double, int>(scanNum,mz,j+2);
                    break;
                }
                _identifiedFeatures[featureId][fileIndex] = det;
            }
            return spectrumList;
        }

        private List<SequenceTag.SequenceTag> GetTags(IReadOnlyCollection<ProductSpectrum> spectra)
        {
            var tagDict = new Dictionary<string,SequenceTag.SequenceTag>();
            if (spectra.Count == 0) return tagDict.Values.ToList();
            foreach (var spectrum in spectra)
            {
                var tagFinder = new SequenceTagFinder(spectrum, new Tolerance(10), 4);
                var tags = tagFinder.GetAllSequenceTagString();
                foreach (var t in tags)
                {
                    if (tagDict.ContainsKey(t.Sequence)) continue;
                    tagDict.Add(t.Sequence,t);
                }
            }
            return tagDict.Values.ToList();
        }

        private int TagsInDatabase(IReadOnlyCollection<SequenceTag.SequenceTag> tags)
        {
            var hits = 0;
            if (tags.Count == 0) return hits;
            foreach (var tag in tags)
            {
                var index = _searchableDB.Search(tag.Sequence);
                if (index >= 0) hits++;
            }
            return hits;
        }

        private void InitializeMatrix(IList<int[]> matrix)
        {
            for (var i = 0; i < matrix.Length; i++)
            {
                matrix[i] = new int[_rawFiles.Length];
            }
        }

        private void GetFilteredFeatures(string filterFile)
        {
            var filterReader = new StreamReader(filterFile);
            string filterLine;
            while (((filterLine = filterReader.ReadLine()) != null))
            {
                var filteredElements = filterLine.Split('\t');
                var id = int.Parse(filteredElements[0]);
                var mass = double.Parse(filteredElements[1]);
                var minElution = double.Parse(filteredElements[2]);
                var maxElution = double.Parse(filteredElements[3]);
                var pVal = double.Parse(filteredElements[4]);
                var fVal = double.Parse(filteredElements[5]);
                _filteredFeatures.Add(new Tuple<int, double, double,double,double,double>(id, mass, minElution, maxElution,pVal,fVal));
            }
        }

        private double[] GetFeatureMassTable(double mass, int chargeLimit = 60)
        {
            return Enumerable.Range(2, chargeLimit - 1).Select(x => ((mass + x) / x)).ToArray();
        }

        private readonly string[] _rawFiles;
        private readonly string _crossTabFile;
        private string _filteredFile;
        private readonly SearchableDatabase _searchableDB;
        private readonly List<Tuple<int,double, double, double,double,double>> _filteredFeatures;
        private Tuple<int, double, int>[][] _identifiedFeatures;
        private int[][] _spectrumMatchesMatrix;
        private int[][] _tagsGeneratedMatrix;
        private int[][] _dataBaseHitMatrix;
    }
}
