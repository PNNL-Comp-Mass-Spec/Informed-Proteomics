using System;
using System.Collections.Generic;
using System.IO;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.IMS.IMSScoring
{
    public class SubScoreFactory
    {
        // IMPORTANT : if ion1 == ion2 for isotope it means relaton between ion1 and precursor!
        // Now isotope score for precursor is not trained (mgf does not contain precursor profile), and node score of precursor is zero should be fixed..

        private const int MaxCorrIntScore = 5;
        private const double CorrIntScoreFactorForIsotopeScore = 0.75; // the lower the more stringent
        private const double CorrIntScoreFactorForIsotopeScoreForPrecursor = 0.6;
        private const double CorrIntScoreFactorForLcImsScore = 0.5;
        private const double LowScoreForIsotopeScore = -3;
        private const double HighScoreForIsotopeScore = 3;
        private const double LowScoreForLcImsScore = -.5;
        private const double HighScoreForLcImsScore = 1.5;


        //TODO : ims, lc corr score should be dependent on ratio scores..

        private readonly Dictionary<GroupParameter, List<IonType>> _ionTypeDictionary;

        private Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _ratioScoreDictionary; // trained
        private readonly Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _ratioTargetProbDictionary; // trained
        private readonly Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> _ratioDecoyProbDictionary; // trained

        private Dictionary<GroupParameter, Dictionary<int, double>> _lcCorrScoreDictionary; 
        private readonly Dictionary<GroupParameter, Dictionary<int, double>> _lcCorrTargetProbDictionary;
        private readonly Dictionary<GroupParameter, Dictionary<int, double>> _lcCorrDecoyProbDictionary;

        private Dictionary<GroupParameter, Dictionary<int, double>> _imsCorrScoreDictionary;
        private readonly Dictionary<GroupParameter, Dictionary<int, double>> _imsCorrTargetProbDictionary;
        private readonly Dictionary<GroupParameter, Dictionary<int, double>> _imsCorrDecoyProbDictionary;

        private readonly Dictionary<GroupParameter, Dictionary<int, double>> _isotopeIntensityCorrScoreDictionary;

        private readonly Dictionary<GroupParameter, double> _noIonScore; // trained

        private readonly Dictionary<GroupParameter, Dictionary<IonType, int>> _ionTypeIndexDictionary; 

        public SubScoreFactory(string filePath)
        {
            _ionTypeDictionary = new Dictionary<GroupParameter, List<IonType>>();
            _ratioTargetProbDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();
            _ratioDecoyProbDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();
            _ratioScoreDictionary = new Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>>();

            _lcCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();
            _lcCorrTargetProbDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();
            _lcCorrDecoyProbDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();

            _imsCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();
            _imsCorrTargetProbDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();
            _imsCorrDecoyProbDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();

            _isotopeIntensityCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();

            _noIonScore = new Dictionary<GroupParameter, double>();

            
            Read(filePath);

            _ionTypeIndexDictionary = new Dictionary<GroupParameter, Dictionary<IonType, int>>();

            foreach (var gp in _ionTypeDictionary.Keys)
            {
                var ions = _ionTypeDictionary[gp];
                _ionTypeIndexDictionary[gp] = new Dictionary<IonType, int>();
                for (var i = 0; i < ions.Count; i++)
                {
                    _ionTypeIndexDictionary[gp][ions[i]] = i;
                }
            }

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
            return GetScoreFromDictionary(_ratioScoreDictionary, ionType1, ionType2, ratio, parameter, LowScoreForIsotopeScore);
        }

        internal double GetLcCorrelationScore(double lcCorr, GroupParameter parameter)
        {
            if (_lcCorrScoreDictionary.Count == 0)
                return CorrToIntForLcImsScore(lcCorr) / (double)MaxCorrIntScore * (HighScoreForLcImsScore - LowScoreForLcImsScore) + LowScoreForLcImsScore;
            return GetScoreFromDictionary(_lcCorrScoreDictionary, CorrToIntForLcImsScore(lcCorr), parameter, LowScoreForLcImsScore);
        }

        internal double GetImsCorrelationScore(double imsCorr, GroupParameter parameter)
        {
            if (_imsCorrScoreDictionary.Count == 0)
                return CorrToIntForLcImsScore(imsCorr) / (double)MaxCorrIntScore * (HighScoreForLcImsScore - LowScoreForLcImsScore) + LowScoreForLcImsScore;
            return GetScoreFromDictionary(_imsCorrScoreDictionary, CorrToIntForLcImsScore(imsCorr), parameter, LowScoreForLcImsScore);
        }

        internal double GetRatioScore(IonType ionType, int ratio, GroupParameter parameter)  // between precursor and frag ion
        {
            return GetRatioScore(ionType, ionType, ratio, parameter);
        }

        internal double GetIsotopeIntensityCorrelationScore(double isotopeCorr, GroupParameter groupParameter)
        {
            if (_isotopeIntensityCorrScoreDictionary.Count == 0)
            {
                var ret =  CorrToIntForIsotopeScore(isotopeCorr)/(double) MaxCorrIntScore*
                       (HighScoreForIsotopeScore - LowScoreForIsotopeScore) + LowScoreForIsotopeScore;
                return ret;
            }
            return GetScoreFromDictionary(_isotopeIntensityCorrScoreDictionary, CorrToIntForIsotopeScore(isotopeCorr), groupParameter, LowScoreForIsotopeScore);
        }

        internal double GetIsotopeIntensityCorrelationScoreForPrecursor(double isotopeCorr, GroupParameter groupParameter)
        {
            if (_isotopeIntensityCorrScoreDictionary.Count == 0)
            {
                var ret = CorrToIntForIsotopeScoreForPrecursor(isotopeCorr) / (double)MaxCorrIntScore *
                       (HighScoreForIsotopeScore - LowScoreForIsotopeScore) + LowScoreForIsotopeScore;
                return ret;
            }
            return GetScoreFromDictionary(_isotopeIntensityCorrScoreDictionary, CorrToIntForIsotopeScoreForPrecursor(isotopeCorr), groupParameter, LowScoreForIsotopeScore);
        }

        internal double GetKLDivergence(IonType ionType, IonType ionType1, int ratio, GroupParameter groupParameter)
        {
            double prob1 = 1, prob2 = 1;
            var si = _ionTypeIndexDictionary[groupParameter];
            var ionTypeIndices = new Tuple<int, int> ( si[ionType], si[ionType1]); 
            if(_ratioTargetProbDictionary.ContainsKey(groupParameter))
            {
                prob1 = prob1 * _ratioTargetProbDictionary[groupParameter][ionTypeIndices][ratio];
                prob2 = prob2 * _ratioDecoyProbDictionary[groupParameter][ionTypeIndices][ratio];
            }
           
            return GetKLDivergence(prob1, prob2);
        }

        internal double GetKLDivergence(IonType ionType, int ratio, GroupParameter groupParameter)
        {
            return GetKLDivergence(ionType, ionType, ratio, groupParameter);
        }

        private double GetKLDivergence(double prob1, double prob2)
        {
            return prob1*Math.Log(prob1/prob2) + (1 - prob1)*Math.Log((1 - prob1)/(1 - prob2));
        }

        public static int CorrToIntForIsotopeScore(double corr) // the higher the better from 1 to 5
        {
            var m = 1.0;
            var score = 0;
            for (; score < MaxCorrIntScore; score++)
            {
                if (1 - m > corr)
                    break;
                m = m * CorrIntScoreFactorForIsotopeScore; // 
            }
            return score;
        }

        public static int CorrToIntForIsotopeScoreForPrecursor(double corr) // the higher the better from 1 to 5
        {
            var m = 1.0;
            var score = 0;
            for (; score < MaxCorrIntScore; score++)
            {
                if (1 - m > corr)
                    break;
                m = m * CorrIntScoreFactorForIsotopeScoreForPrecursor; // 
            }
            return score;
        }

        public static int CorrToIntForLcImsScore(double corr) // the higher the better from 1 to 5
        {
            var m = 1.0;
            var score = 0;
            for (; score < MaxCorrIntScore; score++)
            {
                if (1 - m > corr)
                    break;
                m = m * CorrIntScoreFactorForLcImsScore; // 
            }
            return score;
        }

        public static int[] GetAllCorrIntScore()
        {
            var s = new int[MaxCorrIntScore+1];
            for (var i = 0; i <= MaxCorrIntScore; i++)
                s[i] = i;
            return s;
        }

        private double GetScoreFromDictionary(Dictionary<GroupParameter, Dictionary<Tuple<int, int>, Dictionary<int, double>>> scoreDictionary, IonType ionType1, IonType ionType2, int rawScore, GroupParameter groupParameter, double lowScore)
        {
            if (!scoreDictionary.ContainsKey(groupParameter)) return lowScore;
            var l1 = scoreDictionary[groupParameter];
            var si = _ionTypeIndexDictionary[groupParameter];
            var key = new Tuple<int, int> ( si[ionType1], si[ionType2] );
            if (!l1.ContainsKey(key)) return lowScore;
            var l2 = l1[key];
            return !l2.ContainsKey(rawScore) ? lowScore : l2[rawScore];
        }

        private double GetScoreFromDictionary(Dictionary<GroupParameter, Dictionary<int, double>> scoreDictionary, int rawScore, GroupParameter parameter, double lowScore)
        {
            var key = parameter.GetPrecursorGroupParameter();
            if (!scoreDictionary.ContainsKey(key)) return lowScore;
            var l1 = scoreDictionary[key];
            return !l1.ContainsKey(rawScore) ? lowScore : l1[rawScore];
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

        private Dictionary<GroupParameter, Dictionary<int, double>> GetScoreDictionaryFromProbDictionaries(Dictionary<GroupParameter, Dictionary<int, double>> target, Dictionary<GroupParameter, Dictionary<int, double>> decoy)
        {
            var scoreDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();

            foreach (var groupParameter in target.Keys)
            {
                var st = target[groupParameter];
                var sd = decoy[groupParameter];
                scoreDictionary[groupParameter] = new Dictionary<int, double>();
                foreach (var score in st.Keys)
                {
                    scoreDictionary[groupParameter][score] = GetLogLRScore(st[score], sd[score]);
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
