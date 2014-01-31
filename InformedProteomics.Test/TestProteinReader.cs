using System.IO;
using System.Linq;
using InformedProteomics.Backend.Data.Sequence;
using InformedProteomics.Backend.Utils;
using NUnit.Framework;

namespace InformedProteomics.Test
{
    [TestFixture]
    class TestProteinReader
    {
        private string Trim(string prot)
        {
            int start = prot.IndexOf('.') + 1;
            int length = prot.LastIndexOf('.') - start;
            return prot.Substring(start, length);
        }
        [Test]
        public void Reader()
        {
            const string fileName = @"\\protoapps\UserData\Sangtae\ForChris\E_coli_iscU_60_mock_MSAlign_ResultTable.txt";
            var fileParser = new TsvFileParser(fileName);
            var data = fileParser.GetData("Peptide");

            var aset = new AminoAcidSet();

            var writeFile = new StreamWriter(fileName+".out");
            writeFile.WriteLine("Peptide\tMass\n");

            foreach (var i in data.Where(x => !x.Contains('(') && !x.Contains(')') && !x.Contains('[') && !x.Contains(']')))
            {
                var trimmed = Trim(i);
                var composition = aset.GetComposition(trimmed);
                if (composition != null)
                    writeFile.WriteLine(trimmed + '\t' + composition.GetMass());
            }
            writeFile.Close();
        }
    }
}
