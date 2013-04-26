using System.Collections.Generic;
using InformedProteomics.Backend.Data.Spectrometry;

namespace InformedProteomics.Backend.IMSScoring
{
    public static class SubScoreFactory
    {
        private const double LowScore = -5;

        static private Dictionary<GroupParameter, List<IonType>> _ionTypeDictionary;

        static private Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>> _ratioScoreDictionary;
        static private Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>> _ratioTargetProbDictionary;
        static private Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>> _ratioDecoyProbDictionary;

        static private Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>> _lcCorrScoreDictionary;
        static private Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>> _lcCorrTargetProbDictionary;
        static private Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>> _lcCorrDecoyProbDictionary;

        static private Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>> _imsCorrScoreDictionary;
        static private Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>> _imsCorrTargetProbDictionary;
        static private Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>> _imsCorrDecoyProbDictionary;

        private static Dictionary<GroupParameter, Dictionary<IonType, Dictionary<int, double>>> _isotopeLCCorrScoreDictionary;
        private static Dictionary<GroupParameter, Dictionary<IonType, Dictionary<int, double>>> _isotopeIMSCorrScoreDictionary;
        private static Dictionary<GroupParameter, Dictionary<IonType, Dictionary<int, double>>> _isotopeIntensityCorrScoreDictionary;

        private static Dictionary<GroupParameter, Dictionary<int, double>> _precursorIsotopeLCCorrScoreDictionary;
        private static Dictionary<GroupParameter, Dictionary<int, double>> _precursorIsotopeIMSCorrScoreDictionary;
        private static Dictionary<GroupParameter, Dictionary<int, double>> _precursorIsotopeIntensityCorrScoreDictionary; 

        private static void ClearDictionaries()
        {
            _ionTypeDictionary = new Dictionary<GroupParameter, List<IonType>>();
            _ratioTargetProbDictionary = new Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>>();
            _ratioDecoyProbDictionary = new Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>>();
            _ratioScoreDictionary = new Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>>();

            _lcCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>>();
            _lcCorrTargetProbDictionary = new Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>>();
            _lcCorrDecoyProbDictionary = new Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>>();

            _imsCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>>();
            _imsCorrTargetProbDictionary = new Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>>();
            _imsCorrDecoyProbDictionary = new Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>>();

            _isotopeLCCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<IonType, Dictionary<int, double>>>();
            _isotopeIMSCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<IonType, Dictionary<int, double>>>();
            _isotopeIntensityCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<IonType, Dictionary<int, double>>>();

            _precursorIsotopeLCCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();
            _precursorIsotopeIMSCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();
            _precursorIsotopeIntensityCorrScoreDictionary = new Dictionary<GroupParameter, Dictionary<int, double>>();
        }

        internal static void Read(string fileName)//TODO implement ..!!
        {
            ClearDictionaries();
        }

        internal static List<IonType> GetIonTypes(GroupParameter parameter)
        {
            if(_ionTypeDictionary.ContainsKey(parameter))
                return _ionTypeDictionary[parameter];
            return new List<IonType>();
        }

        internal static double GetRatioScore(IonType ionType1, IonType ionType2, int ratio, GroupParameter parameter)
        {
            return GetScoreFromDictionary(_ratioScoreDictionary, ionType1, ionType2, ratio, parameter);
        }
 
        internal static double GetLCCorrelationScore(IonType ionType1, IonType ionType2, double lcCorr, GroupParameter parameter)
        {
            if (_lcCorrScoreDictionary.Count == 0)
                return lcCorr*(-2*LowScore) + LowScore;
            return GetScoreFromDictionary(_lcCorrScoreDictionary, ionType1, ionType2, CorrToInt(lcCorr), parameter);
        }

        internal static double GetIMSCorrelationScore(IonType ionType1, IonType ionType2, double imsCorr, GroupParameter parameter)
        {
            if (_imsCorrScoreDictionary.Count == 0)
                return imsCorr * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_imsCorrScoreDictionary, ionType1, ionType2, CorrToInt(imsCorr), parameter);
        }

        internal static double GetRatioScore(IonType ionType, int ratio, GroupParameter parameter)  // between precursor and frag ion
        {
            return GetRatioScore(ionType, ionType, ratio, parameter);
        }

        internal static double GetLCCorrelationScore(IonType ionType, double lcCorr, GroupParameter parameter)  // between precursor and frag ion
        {
            return GetLCCorrelationScore(ionType, ionType, lcCorr, parameter);
        }

        internal static double GetIMSCorrelationScore(IonType ionType, double imsCorr, GroupParameter parameter) // between precursor and frag ion
        {
            return GetIMSCorrelationScore(ionType, ionType, imsCorr, parameter);
        }

        internal static double GetIsotopeIntensityCorrelationScore(IonType ionType, double isotopeCorr, GroupParameter parameter)
        {
            return GetScoreFromDictionary(_isotopeIntensityCorrScoreDictionary, ionType, CorrToInt(isotopeCorr), parameter);
        }

        internal static double GetIsotopeLCCorrelationScore(IonType ionType, double lcCorr, GroupParameter parameter)
        {
            if (_isotopeLCCorrScoreDictionary.Count == 0)
                return lcCorr * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_isotopeLCCorrScoreDictionary, ionType, CorrToInt(lcCorr), parameter);
        }

        internal static double GetIsotopeIMSCorrelationScore(IonType ionType, double imsCorr, GroupParameter parameter)
        {
            if (_isotopeIMSCorrScoreDictionary.Count == 0)
                return imsCorr * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_isotopeIMSCorrScoreDictionary, ionType, CorrToInt(imsCorr), parameter);
        }

        internal static double GetIsotopeIntensityCorrelationScore(double isotopeCorr, GroupParameter parameter)
        {
            return GetScoreFromDictionary(_precursorIsotopeIntensityCorrScoreDictionary, CorrToInt(isotopeCorr), parameter);
        }

        internal static double GetIsotopeLCCorrelationScore(double lcCorr, GroupParameter parameter)
        {
            if (_precursorIsotopeLCCorrScoreDictionary.Count == 0)
                return lcCorr * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_precursorIsotopeLCCorrScoreDictionary, CorrToInt(lcCorr), parameter);
        }

        internal static double GetIsotopeIMSCorrelationScore(double imsCorr, GroupParameter parameter)
        {
            if (_precursorIsotopeIMSCorrScoreDictionary.Count == 0)
                return imsCorr * (-2 * LowScore) + LowScore;
            return GetScoreFromDictionary(_precursorIsotopeIMSCorrScoreDictionary, CorrToInt(imsCorr), parameter);
        }

        internal static double GetKLDivergence(IonType ionType, IonType ionType1, int ratio, double lcCorrelation, double imsCorrelation, GroupParameter parameter)
        {
            throw new System.NotImplementedException();
        }

        internal static double GetKLDivergence(IonType ionType, int ratio, double lcCorr, double imsCorr, GroupParameter parameter)
        {
            throw new System.NotImplementedException();
        }

        private static int CorrToInt(double corr)
        {
            return (int) corr*10;
        }

        private static double GetScoreFromDictionary(Dictionary<GroupParameter, Dictionary<List<IonType>, Dictionary<int, double>>> scoreDictionary, IonType ionType1, IonType ionType2, int rawScore, GroupParameter parameter)
        {
            if (!scoreDictionary.ContainsKey(parameter)) return LowScore;
            var l1 = scoreDictionary[parameter];
            var key = new List<IonType> { ionType1, ionType2 };
            if (!l1.ContainsKey(key)) return LowScore;
            var l2 = l1[key];
            return !l2.ContainsKey(rawScore) ? LowScore : l2[rawScore];
        }

        private static double GetScoreFromDictionary(Dictionary<GroupParameter, Dictionary<IonType, Dictionary<int, double>>> scoreDictionary, IonType ionType, int rawScore, GroupParameter parameter)
        {
            if (!scoreDictionary.ContainsKey(parameter)) return LowScore;
            var l1 = scoreDictionary[parameter];
            if (!l1.ContainsKey(ionType)) return LowScore;
            var l2 = l1[ionType];
            return !l2.ContainsKey(rawScore) ? LowScore : l2[rawScore];
        }

        private static double GetScoreFromDictionary(Dictionary<GroupParameter, Dictionary<int, double>> scoreDictionary, int rawScore, GroupParameter parameter)
        {
            var key = parameter.GetPrecursorGroupParameter();
            if (!scoreDictionary.ContainsKey(key)) return LowScore;
            var l1 = scoreDictionary[key];
            return !l1.ContainsKey(rawScore) ? LowScore : l1[rawScore];
        }
    }
}
