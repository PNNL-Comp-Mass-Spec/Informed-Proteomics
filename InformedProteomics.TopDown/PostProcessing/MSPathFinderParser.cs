using System;
using System.Collections.Generic;
using System.Linq;
using InformedProteomics.Backend.Data.Composition;
using InformedProteomics.Backend.Utils;

namespace InformedProteomics.TopDown.PostProcessing
{
    public class MsPathFinderParser
    {
        public MsPathFinderParser(string fileName)
        {
            _idList = new List<MsPathFinderId>();
            _scanNumToPrSm = new Dictionary<int, MsPathFinderId>();
            if (!Parse(fileName)) _scanNumToPrSm = null;
        }

        public MsPathFinderId GetId(int scanNum)
        {
            MsPathFinderId id;
            if (_scanNumToPrSm != null && _scanNumToPrSm.TryGetValue(scanNum, out id)) return id;
            return null;
        }

        public IList<MsPathFinderId> GetIdWithQValuesNoLargerThan(double qValueThreshold)
        {
            return _idList.TakeWhile(id => id.QValue <= qValueThreshold).ToList();
        }

        public IList<MsPathFinderId> GetIdList()
        {
            return _idList;
        }

        private readonly IList<MsPathFinderId> _idList;
        private readonly Dictionary<int, MsPathFinderId> _scanNumToPrSm;
        private bool Parse(string fileName)
        {
            var parser = new TsvFileParser(fileName);
            var scan = parser.GetData("Scan").Select(s => Convert.ToInt32(s)).ToArray();
            var pre = parser.GetData("Pre").Where(s => s.Length == 1).Select(p => p[0]).ToArray();
            if (pre.Length != parser.NumData) return false;
            var sequence = parser.GetData("Sequence").ToArray();
            var post = parser.GetData("Post").Where(s => s.Length == 1).Select(p => p[0]).ToArray();
            if (post.Length != parser.NumData) return false;
            var mod = parser.GetData("Modifications").ToArray();
            var composition = parser.GetData("Composition").Select(Composition.Parse).ToArray();
            var proteinName = parser.GetData("ProteinName").ToArray();
            var proteinDesc = parser.GetData("ProteinDesc").ToArray();
            var proteinLength = parser.GetData("ProteinLength").Select(s => Convert.ToInt32(s)).ToArray();
            var start = parser.GetData("Start").Select(s => Convert.ToInt32(s)).ToArray();
            var end = parser.GetData("End").Select(s => Convert.ToInt32(s)).ToArray();
            var charge = parser.GetData("Charge").Select(s => Convert.ToInt32(s)).ToArray();
            var mostAbundantIsotopeMz = parser.GetData("MostAbundantIsotopeMz").Select(Convert.ToDouble).ToArray();
            var mass = parser.GetData("Mass").Select(Convert.ToDouble).ToArray();
            var numMatchedFragment = parser.GetData("#MatchedFragments").Select(s => Convert.ToInt32(s)).ToArray();
            var qValue = parser.GetData("QValue").Select(Convert.ToDouble).ToArray();
            var pepQValue = parser.GetData("PepQValue").Select(Convert.ToDouble).ToArray();

            for (var i = 0; i < parser.NumData; i++)
            {
                var id = new MsPathFinderId(scan[i], pre[i], sequence[i], post[i], mod[i],
                    composition[i], proteinName[i], proteinDesc[i], proteinLength[i],
                    start[i], end[i], charge[i], mostAbundantIsotopeMz[i], mass[i],
                    numMatchedFragment[i], qValue[i], pepQValue[i])
                ;
                _idList.Add(id);

                if(!_scanNumToPrSm.ContainsKey(scan[i])) _scanNumToPrSm.Add(scan[i], id);
            }
            return true;
        }
    }
}
