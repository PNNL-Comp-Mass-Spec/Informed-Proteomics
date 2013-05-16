using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.IMSScoring
{
    public class SubScoreFactory
    {
        // IMPORTANT : if ion1 == ion2 for isotope it means relaton between ion1 and precursor!
        // Now isotope score for precursor is not trained (mgf does not contain precursor profile), and node score of precursor is zero should be fixed..

        private const double LowScore = -3; // TODO how to??
        //TODO : ims, lc corr score should be dependent on ratio scores..

        private readonly Dictionary<GroupParameter, List<IonType>> _ionTypeDictionary;

        private Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _ratioScoreDictionary; // trained
        private readonly Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _ratioTargetProbDictionary; // trained
        private readonly Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _ratioDecoyProbDictionary; // trained

        private Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _lcCorrScoreDictionary; 
        private readonly Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _lcCorrTargetProbDictionary;
        private readonly Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _lcCorrDecoyProbDictionary;

        private Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _imsCorrScoreDictionary;
        private readonly Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _imsCorrTargetProbDictionary;
        private readonly Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _imsCorrDecoyProbDictionary;

        private readonly Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>> _isotopeLcCorrScoreDictionary;
        private readonly Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>> _isotopeImsCorrScoreDictionary;
        private Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>> _isotopeIntensityCorrScoreDictionary; // trained

        private readonly Dictionary<GroupParameter, Dictionary<int, double>> _precursorIsotopeLcCorrScoreDictionary;
        private readonly Dictionary<GroupParameter, Dictionary<int, double>> _precursorIsotopeImsCorrScoreDictionary;
        private readonly Dictionary<GroupParameter, Dictionary<int, double>> _precursorIsotopeIntensityCorrScoreDictionary;

        private readonly Dictionary<GroupParameter, double> _noIonScore; // trained

        public SubScoreFactory(string filePath)
        {
            _ionTypeDictionary = new Dictionary<GroupParameter, List<IonType>>();
            _ratioTargetProbDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();
            _ratioDecoyProbDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();
            _ratioScoreDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();

            _lcCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();
            _lcCorrTargetProbDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();
            _lcCorrDecoyProbDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();

            _imsCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();
            _imsCorrTargetProbDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();
            _imsCorrDecoyProbDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();

            _isotopeLcCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>>();
            _isotopeImsCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>>();
            _isotopeIntensityCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>>();

            _precursorIsotopeLcCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();
            _precursorIsotopeImsCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();
            _precursorIsotopeIntensityCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();

            _noIonScore = new Dictionary<GroupParameter, double>();

            Read(filePath);
        }


        private void Read(string filePath)
        {
            var ratioProbDictionary = _ratioTargetProbDictionary;
            var isotopeIntensityCorrTargetProbDictionary = new Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>>();
            var isotopeIntensityCorrDecoyProbDictionary = new Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>>();
            var isotopeIntensityCorrProbDictionary = isotopeIntensityCorrTargetProbDictionary;
            
            var noIonTargetProbDictionary = new Dictionary<GroupParameter, double>();
            var noIonDecoyProbDictionary = new Dictionary<GroupParameter, double>();
            var noIonProbDictionary = noIonTargetProbDictionary;

            var stremaReader = new StreamReader(filePath);
            string s;
            var mode = -1;
            GroupParameter groupParameter = null;
            List<int> ionTypeIndices = null;
            while ((s = stremaReader.ReadLine()) != null)
            {
                if (s.StartsWith("###DECOY"))
                {
                    ratioProbDictionary = _ratioDecoyProbDictionary;
                    isotopeIntensityCorrProbDictionary = isotopeIntensityCorrDecoyProbDictionary;
                    noIonProbDictionary = noIonDecoyProbDictionary;
                    continue;
                }
                if (s.StartsWith("##IONTYPES"))
                {
                    mode = 1;
                    continue;
                }
                if (s.StartsWith("##ISOTOPE"))
                {
                    mode = 2;
                    continue;
                }
                if (s.StartsWith("##RATIO"))
                {
                    mode = 3;
                    continue;
                }
                if (s.StartsWith("##NOION"))
                {
                    mode = 4;
                    continue;
                }

                var token = s.Split('\t');
                if (s.StartsWith("#G"))
                {
                    groupParameter = GroupParameter.Parse(token[1]);
                    continue;
                }

                if (s.StartsWith("#I"))
                {
                    if (mode == 1)
                    {
                        var ionType = IonType.Parse(token[1]);
                        if (groupParameter == null || ionType == null) continue;
                        if (!_ionTypeDictionary.ContainsKey(groupParameter))
                            _ionTypeDictionary[groupParameter] = new List<IonType>();
                        _ionTypeDictionary[groupParameter].Add(ionType);
                    }
                    else
                    {
                        ionTypeIndices = new List<int>();
                        for (var i = 1; i < token.Length; i++)
                            ionTypeIndices.Add(int.Parse(token[i]));
                    }
                    continue;
                }
                if (!s.StartsWith("#") && token.Length > 0)
                {
                    if (mode == 2)
                    {
                        if (groupParameter == null || ionTypeIndices == null || ionTypeIndices.Count == 0) continue;
                        if(!isotopeIntensityCorrProbDictionary.ContainsKey(groupParameter))
                            isotopeIntensityCorrProbDictionary[groupParameter] = new Dictionary<int, Dictionary<int, double>>();
                        var si = isotopeIntensityCorrProbDictionary[groupParameter];
                        if (!si.ContainsKey(ionTypeIndices[0])) si[ionTypeIndices[0]] = new Dictionary<int, double>();
                        var ssi = si[ionTypeIndices[0]];
                        foreach(var scoreStr in token)
                        {
                            var st = scoreStr.Split(',');
                            if(st.Length ==2)
                                ssi[int.Parse(st[0])] = double.Parse(st[1]);
                        }
                    }else if (mode == 3)
                    {
                        if (groupParameter == null || ionTypeIndices == null || ionTypeIndices.Count < 2) continue;
                        if (!ratioProbDictionary.ContainsKey(groupParameter))
                            ratioProbDictionary[groupParameter] = new Dictionary<Tuple<int, int>, Dictionary<int, double>>();
                        var sr = ratioProbDictionary[groupParameter];
                        var key = new Tuple<int, int>(ionTypeIndices[0], ionTypeIndices[1]);
                        if (!sr.ContainsKey(key)) sr[key] = new Dictionary<int, double>();
                        var ssr = sr[key];
                        foreach (var scoreStr in token)
                        {
                            var st = scoreStr.Split(',');
                            if(st.Length == 2)
                                ssr[int.Parse(st[0])] = double.Parse(st[1]);
                        }
                    }else if (mode == 4)
                    {
                        if (groupParameter == null) continue;
                        noIonProbDictionary[groupParameter] = double.Parse(s);
                    }
                } 
            }
            
            foreach (var k in noIonTargetProbDictionary.Keys)
            {
                _noIonScore[k] = GetLogLRScore(noIonTargetProbDictionary[k], noIonDecoyProbDictionary[k]);
            }

            _isotopeIntensityCorrScoreDictionary = GetScoreDictionaryFromProbDictionaries(isotopeIntensityCorrTargetProbDictionary, isotopeIntensityCorrDecoyProbDictionary);
            GetScoreDictionariesFromProbDictionaries();
        }

        internal double GetNoIonScore(GroupParameter parameter)
        {
            return _noIonScore[parameter];
        }

        internal List<IonType> GetIonTypes(GroupParameter parameter)
        {
            if(_ionTypeDictionary.ContainsKey(parameter))
                return _ionTypeDictionary[parameter];
            return new List<IonType>();
        }

        internal double GetRatioScore(IonType ionType1, IonType ionType2, int ratio, GroupParameter parameter)
        {
            return GetScoreFromDictionary(_ratioScoreDictionary, ionType1, ionType2, ratio, parameter);
        }

        internal double GetLcCorrelationScore(IonType ionType1, IonType ionType2, double lcCorr, GroupParameter parameter)
        {
            if (_lcCorrScoreDictionary.Count == 0)
                return CorrToInt(lcCorr) / 5.0 * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_lcCorrScoreDictionary, ionType1, ionType2, CorrToInt(lcCorr), parameter);
        }

        internal double GetImsCorrelationScore(IonType ionType1, IonType ionType2, double imsCorr, GroupParameter parameter)
        {
            if (_imsCorrScoreDictionary.Count == 0)
                return CorrToInt(imsCorr) / 5.0 * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_imsCorrScoreDictionary, ionType1, ionType2, CorrToInt(imsCorr), parameter);
        }

        internal double GetRatioScore(IonType ionType, int ratio, GroupParameter parameter)  // between precursor and frag ion
        {
            return GetRatioScore(ionType, ionType, ratio, parameter);
        }

        internal double GetLcCorrelationScore(IonType ionType, double lcCorr, GroupParameter parameter)  // between precursor and frag ion
        {
            return GetLcCorrelationScore(ionType, ionType, lcCorr, parameter);
        }

        internal double GetImsCorrelationScore(IonType ionType, double imsCorr, GroupParameter parameter) // between precursor and frag ion
        {
            return GetImsCorrelationScore(ionType, ionType, imsCorr, parameter);
        }

        internal double GetIsotopeIntensityCorrelationScore(IonType ionType, double isotopeCorr, GroupParameter groupParameter)
        {
            if (_isotopeIntensityCorrScoreDictionary.Count == 0)
                return CorrToInt(isotopeCorr) / 5.0 * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_isotopeIntensityCorrScoreDictionary, _ionTypeDictionary[groupParameter].IndexOf(ionType), CorrToInt(isotopeCorr), groupParameter);
        }

        internal double GetIsotopeLcCorrelationScore(IonType ionType, double lcCorr, GroupParameter parameter)
        {
            if (_isotopeLcCorrScoreDictionary.Count == 0)
                return CorrToInt(lcCorr) / 5.0 * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_isotopeLcCorrScoreDictionary, _ionTypeDictionary[parameter].IndexOf(ionType), CorrToInt(lcCorr), parameter);
        }

        internal double GetIsotopeImsCorrelationScore(IonType ionType, double imsCorr, GroupParameter groupParameter)
        {
            if (_isotopeImsCorrScoreDictionary.Count == 0)
                return CorrToInt(imsCorr) / 5.0 * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_isotopeImsCorrScoreDictionary, _ionTypeDictionary[groupParameter].IndexOf(ionType), CorrToInt(imsCorr), groupParameter);
        }

        internal double GetIsotopeIntensityCorrelationScore(double isotopeCorr, GroupParameter groupParameter)
        {
            if (_precursorIsotopeIntensityCorrScoreDictionary.Count == 0)
                return CorrToInt(isotopeCorr) / 5.0 * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_precursorIsotopeIntensityCorrScoreDictionary, CorrToInt(isotopeCorr), groupParameter);
        }

        internal double GetIsotopeLcCorrelationScore(double lcCorr, GroupParameter parameter)
        {
            if (_precursorIsotopeLcCorrScoreDictionary.Count == 0)
                return CorrToInt(lcCorr) / 5.0 * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_precursorIsotopeLcCorrScoreDictionary, CorrToInt(lcCorr), parameter);
        }

        internal double GetIsotopeImsCorrelationScore(double imsCorr, GroupParameter parameter)
        {
            if (_precursorIsotopeImsCorrScoreDictionary.Count == 0)
                return CorrToInt(imsCorr)/5.0 * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_precursorIsotopeImsCorrScoreDictionary, CorrToInt(imsCorr), parameter);
        }

        internal double GetKLDivergence(IonType ionType, IonType ionType1, int ratio, double lcCorrelation, double imsCorrelation, GroupParameter groupParameter)
        {
            double prob1 = 1, prob2 = 1;
            var si = _ionTypeDictionary[groupParameter];
            var ionTypeIndices = new Tuple<int, int> ( si.IndexOf(ionType), si.IndexOf(ionType1) );
            if(_ratioTargetProbDictionary.ContainsKey(groupParameter))
            {
                prob1 = prob1 * _ratioTargetProbDictionary[groupParameter][ionTypeIndices][ratio];
                prob2 = prob2 * _ratioDecoyProbDictionary[groupParameter][ionTypeIndices][ratio];
            }
            if (_lcCorrTargetProbDictionary.ContainsKey(groupParameter))
            {
                prob1 = prob1 * _lcCorrTargetProbDictionary[groupParameter][ionTypeIndices][CorrToInt(lcCorrelation)];
                prob2 = prob2 * _lcCorrDecoyProbDictionary[groupParameter][ionTypeIndices][CorrToInt(lcCorrelation)];    
            }
            if(_imsCorrTargetProbDictionary.ContainsKey(groupParameter))
            {
                prob1 = prob1 * _imsCorrTargetProbDictionary[groupParameter][ionTypeIndices][CorrToInt(imsCorrelation)];
                prob2 = prob2 * _imsCorrDecoyProbDictionary[groupParameter][ionTypeIndices][CorrToInt(imsCorrelation)];   
            }
            return GetKLDivergence(prob1, prob2);
        }

        internal double GetKLDivergence(IonType ionType, int ratio, double lcCorrelation, double imsCorrelation, GroupParameter groupParameter)
        {
            return GetKLDivergence(ionType, ionType, ratio, lcCorrelation, imsCorrelation, groupParameter);
        }

        private double GetKLDivergence(double prob1, double prob2)
        {
            return prob1*Math.Log(prob1/prob2) + (1 - prob1)*Math.Log((1 - prob1)/(1 - prob2));
        }

        public static int CorrToInt(double corr) // the higher the better from 1 to 5
        {
            var m = 1.0;
            var score = 0;
            for (; score < 5; score++)
            {
                if (1 - m > corr)
                    break;
                m = m * 0.4;
            }
            return score;
        }

        public static int[] GetAllCorrIntScore()
        {
            return new []{1,2,3,4,5};
        }

        private double GetScoreFromDictionary(Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> scoreDictionary, IonType ionType1, IonType ionType2, int rawScore, GroupParameter groupParameter)
        {
            if (!scoreDictionary.ContainsKey(groupParameter)) return LowScore;
            var l1 = scoreDictionary[groupParameter];
            var si = _ionTypeDictionary[groupParameter];
            var key = new Tuple<int, int> ( si.IndexOf(ionType1), si.IndexOf(ionType2) );
            if (!l1.ContainsKey(key)) return LowScore;
            var l2 = l1[key];
            return !l2.ContainsKey(rawScore) ? LowScore : l2[rawScore];
        }

        private double GetScoreFromDictionary(Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>> scoreDictionary, int ionTypeIndex, int rawScore, GroupParameter parameter)
        {
            if (!scoreDictionary.ContainsKey(parameter)) return LowScore;
            var l1 = scoreDictionary[parameter];
            if (!l1.ContainsKey(ionTypeIndex)) return LowScore;
            var l2 = l1[ionTypeIndex];
            return !l2.ContainsKey(rawScore) ? LowScore : l2[rawScore];
        }

        private double GetScoreFromDictionary(Dictionary<GroupParameter, Dictionary<int, double>> scoreDictionary, int rawScore, GroupParameter parameter)
        {
            var key = parameter.GetPrecursorGroupParameter();
            if (!scoreDictionary.ContainsKey(key)) return LowScore;
            var l1 = scoreDictionary[key];
            return !l1.ContainsKey(rawScore) ? LowScore : l1[rawScore];
        }

        private Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> GetScoreDictionaryFromProbDictionaries(Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> target, Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> decoy)
        {
            var scoreDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();

            foreach (var groupParameter in target.Keys)
            {
                var st = target[groupParameter];
                var sd = decoy[groupParameter];
                scoreDictionary[groupParameter] = new Dictionary<Tuple<int, int>, Dictionary<int, double>>();
                foreach (var ionTypes in st.Keys)
                {
                    var sst = st[ionTypes];
                    var ssd = sd[ionTypes];
                    scoreDictionary[groupParameter][ionTypes] = new Dictionary<int, double>();
                    foreach (var score in sst.Keys)
                    {
                        scoreDictionary[groupParameter][ionTypes][score] = GetLogLRScore(sst[score], ssd[score]);
                    }
                }
            }
            return scoreDictionary;
        }

        private Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>> GetScoreDictionaryFromProbDictionaries(Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>> target, Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>> decoy)
        {
            var scoreDictionary = new Dictionary<GroupParameter, Dictionary<int, Dictionary<int, double>>>();

            foreach (var groupParameter in target.Keys)
            {
                var st = target[groupParameter];
                var sd = decoy[groupParameter];
                scoreDictionary[groupParameter] = new Dictionary<int, Dictionary<int, double>>();
                foreach (var ionTypes in st.Keys)
                {
                    var sst = st[ionTypes];
                    var ssd = sd[ionTypes];
                    scoreDictionary[groupParameter][ionTypes] = new Dictionary<int, double>();
                    foreach (var score in sst.Keys)
                    {
                        scoreDictionary[groupParameter][ionTypes][score] = GetLogLRScore(sst[score], ssd[score]);
                    }
                }
            }
            return scoreDictionary;
        }

        private void GetScoreDictionariesFromProbDictionaries()
        {
            _ratioScoreDictionary = GetScoreDictionaryFromProbDictionaries(_ratioTargetProbDictionary, _ratioDecoyProbDictionary);
            _lcCorrScoreDictionary = GetScoreDictionaryFromProbDictionaries(_lcCorrTargetProbDictionary, _lcCorrDecoyProbDictionary);
            _imsCorrScoreDictionary = GetScoreDictionaryFromProbDictionaries(_imsCorrTargetProbDictionary, _imsCorrDecoyProbDictionary);
        }

        private double GetLogLRScore(double targetProb, double decoyProb)
        {
            return Math.Log(targetProb/decoyProb);
        }
    }
}
